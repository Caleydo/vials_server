import os
from flask import json
from intervaltree import IntervalTree
import math
from misopy.parse_gene import parseGene
from misopy.sashimi_plot.plot_utils.plot_gene import readsToWiggle_pysam
import numpy as np
import pysam
import scipy
from scipy.signal import resample
import sqlite3 as sqlite
from vials_server.mapping_service import MappingService
import vials_server.vials_server_helper as v_helper
__author__ = 'Hendrik Strobelt'


vials_db_name = 'all_miso_summaries.sqlite'
miso_sample_info_file = 'samples.json'

class MisoHandler:
    def __init__(self, project_info):
        self.name = 'MisoHandler'
        self.data_type = 'miso'
        self.project_info = project_info

    def enrich_gene_names(self, info_object):

        # detect ID type
        con = sqlite.connect(os.path.join(self.project_info['dir'], vials_db_name))
        cur = con.cursor()
        cur.execute("SELECT name FROM event_names LIMIT 0,100")
        samples = [item for sublist in cur.fetchall() for item in sublist]  # flatten result, faster than map(lambda x: x[0],l)
        con.close()
        id_type = MappingService.guess_id_type(samples)
        info_object['id_type_guessed'] = id_type

        print('mapping...')
        #add genename
        if not (id_type=='genename'):
            mapping_handler = MappingService.get_handler(id_type,'genename')

            dbFile, dbName, fromField, toField, description = mapping_handler.get_sql_fields()
            con = sqlite.connect(os.path.join(self.project_info['dir'], vials_db_name))
            cur = con.cursor()

            cur.execute("ATTACH ? as mapper", (dbFile,))

            cur.execute("DROP TABLE IF EXISTS event_names_enriched ")
            cur.execute("CREATE TABLE event_names_enriched (name text, alt_name text, description text)")

            question_string = "INSERT INTO event_names_enriched (name, alt_name, description) " \
                              "SELECT b.name, a."+toField+", a."+description+\
                              " FROM event_names as b LEFT JOIN mapper."+dbName+ " as a ON b.name==a."+fromField

            cur.execute("CREATE INDEX all_names ON event_names_enriched(name,alt_name)")

            cur.execute(question_string)
            con.commit()
            con.close()
            info_object['gene_name_mapped'] = 'event_names_enriched'
        else:
            info_object['gene_name_mapped'] = 'event_names'

    def get_genes_in_project(self):
        if 'gene_name_mapped' not in self.project_info['info']:
            self.enrich_gene_names(self.project_info['info'])
            v_helper.update_project_file(self.project_info)

        con = sqlite.connect(os.path.join(self.project_info['dir'], vials_db_name))
        cur = con.cursor()

        if self.project_info['info']['gene_name_mapped'] == 'event_names':
            cur.execute("SELECT name FROM event_names")
            res = map(lambda x: [x[0], None, None], cur.fetchall())
        else:
            cur.execute("SELECT name, alt_name, description FROM event_names_enriched")
            res = cur.fetchall()
        # res = [item for sublist in cur.fetchall() for item in sublist]  # flatten result, faster than map(lambda x: x[0],l)


        con.close()

        return res

    def setify(self, seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x['id'] in seen or seen_add(x['id']))]

    def get_genes_in_project_filtered(self, name_filter, exact_match):
        if 'gene_name_mapped' not in self.project_info['info']:
            self.enrich_gene_names(self.project_info['info'])
            v_helper.update_project_file(self.project_info)

        con = sqlite.connect(os.path.join(self.project_info['dir'], vials_db_name))
        cur = con.cursor()

        like_string_a = name_filter+'%'
        like_string_b = '%' + name_filter + '%'

        if exact_match:
            like_string_a = name_filter
            like_string_b = name_filter

        if self.project_info['info']['gene_name_mapped'] == 'event_names':
            cur.execute("SELECT name FROM event_names WHERE alt_name LIKE ? ORDER BY name LIMIT 0,100",
                        (like_string_a,))
            fetched = map(lambda x: [x[0], None, None], cur.fetchall())
        else:
            cur.execute("SELECT name, alt_name, description FROM event_names_enriched WHERE alt_name LIKE ? ORDER BY alt_name LIMIT 0,100",
                        (like_string_a,))

            fetched = cur.fetchall()

            if len(fetched) < 100:
                cur.execute("SELECT name, alt_name, description FROM event_names_enriched WHERE name LIKE ? ORDER BY name LIMIT 0,?",
                            (like_string_b, 100-len(fetched),))
                fetched = fetched + cur.fetchall()

        # res = [item for sublist in cur.fetchall() for item in sublist]  # flatten result, faster than map(lambda x: x[0],l)
        res = map(lambda x: {'id': x[0], 'name': x[1], 'desc': x[2]}, fetched)

        # remove all duplicates from the results (e.g. if name == id)
        res = self.setify(res)

        con.close()

        return res

    def downsample(self,x, size):
        if x.size > size:
            ds_factor = math.floor(float(x.size)/ size)
            fill_size = ds_factor*size - x.size
            if fill_size > 0:
                x = np.append(x, np.zeros(fill_size)*np.NaN)
            else:
                x = np.resize(x, ds_factor*size)
            return scipy.nanmean(x.reshape(-1, ds_factor), axis=1)
        else:
            return x

    def get_wiggles_and_jxns(self, bam_file, gene_meta, downsample_size):

        print "BAM:", bam_file
        bamdata = pysam.Samfile(bam_file, 'rb')
        sample_reads = bamdata.fetch('chr'+str(gene_meta['chromID']), gene_meta['start'], gene_meta['end'])
        wiggle, sample_jxns = readsToWiggle_pysam(sample_reads, gene_meta['start'], gene_meta['end'])
        return sample_jxns, self.downsample(wiggle, downsample_size).tolist()

    def get_samples_and_measures(self, gene_id, gene_meta):
        con = sqlite.connect(os.path.join(self.project_info['dir'], vials_db_name))
        con.row_factory = sqlite.Row
        cur = con.cursor()

        cur.execute("SELECT * FROM miso_summaries WHERE event_name = ?",(gene_id,))
        all_samples = cur.fetchall()
        # res = [item for sublist in cur.fetchall() for item in sublist]  # flatten result, faster than map(lambda x: x[0],l)

        con.close()

        samples_info = {}
        iso_measures = []

        # create a mapping to use short names and reduce bandwidth requirements
        isoform_long_to_short_names = {}
        if 'isoforms' in gene_meta:
            for short_name, info in gene_meta['isoforms'].iteritems():
                isoform_long_to_short_names[info['exonNames']] = short_name

        with open(os.path.join(self.project_info['dir'], miso_sample_info_file)) as samples_file:
            sample_meta = json.load(samples_file)

        jxns = []
        wiggle_list =[]
        for sample in all_samples:

            sample_id = sample['sample']
            samples_info[sample_id] = {
                'id': sample_id,
                'type': 'miso',
                'meta': {}
            }

            isoforms = sample['isoforms'].split(',')
            isoforms = map(lambda iso: iso.strip("'"), isoforms)
            isoforms = map(lambda iso: isoform_long_to_short_names[iso] if (iso in isoform_long_to_short_names) else iso, isoforms)

            measures = sample['miso_posterior_mean'].split(',')

            if len(isoforms)==len(measures):
                for index, isoform in enumerate(isoforms):
                    iso_measures.append({
                        'id': isoform,
                        'sample': sample_id,
                        'weight': measures[index]
                    })

            # cache the complictaed results
            cache_jw_db_path = os.path.join(self.project_info['dir'], '_cache')
            if not os.path.isdir(cache_jw_db_path):
                os.makedirs(cache_jw_db_path)
            con = sqlite.connect(os.path.join(cache_jw_db_path, 'jxn_wiggle.sqlite'))
            con.row_factory = sqlite.Row
            cur = con.cursor()
            cur.execute('CREATE TABLE IF NOT EXISTS jxn_wiggle(id TEXT, jxn TEXT, wiggles BLOB)')



            # try to find the bam file
            if sample_id in sample_meta:
                bam_file = sample_meta[sample_id]['bam_file']
                if bam_file and not bam_file == '':
                    print gene_id+'__'+sample_id

                    cur.execute('SELECT * FROM jxn_wiggle WHERE id=?',(gene_id+'__'+sample_id,))
                    cache_hit = cur.fetchone()

                    if cache_hit:
                        jxn = json.loads(cache_hit[1])
                        wiggles = map(float, cache_hit[2].split(','))

                    else:
                        jxn, wiggles = self.get_wiggles_and_jxns(
                            os.path.join(self.project_info['info']['bam_root_dir'], bam_file),
                            gene_meta,
                            2000)
                        cur.execute('INSERT INTO jxn_wiggle VALUES (?,?,?)',
                                    (gene_id+'__'+sample_id, json.dumps(jxn), ','.join(map(str, wiggles)),))
                        con.commit()
                    wiggle_list.append({'sample':sample_id, 'data': wiggles})
                    jxns.append({'sample': sample_id, 'data': jxn})

            con.close()






        return samples_info, iso_measures, jxns, wiggle_list