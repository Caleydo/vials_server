import os
from flask import json
from intervaltree import IntervalTree
from misopy.parse_gene import parseGene
from scipy.signal import resample
import sqlite3 as sqlite
from vials_server.mapping_service import MappingService
import vials_server.vials_server_helper as v_helper
__author__ = 'Hendrik Strobelt'


vials_db_name = 'all_miso_summaries.sqlite'


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

    def get_genes_in_project_filtered(self, name_filter):
        if 'gene_name_mapped' not in self.project_info['info']:
            self.enrich_gene_names(self.project_info['info'])
            v_helper.update_project_file(self.project_info)

        con = sqlite.connect(os.path.join(self.project_info['dir'], vials_db_name))
        cur = con.cursor()

        if self.project_info['info']['gene_name_mapped'] == 'event_names':
            cur.execute("SELECT name FROM event_names WHERE alt_name LIKE ? LIMIT 0,100",
                        (name_filter+'%',))
            fetched = map(lambda x: [x[0], None, None], cur.fetchall())
        else:
            cur.execute("SELECT name, alt_name, description FROM event_names_enriched WHERE alt_name LIKE ? LIMIT 0,100",
                        (name_filter+'%',))

            fetched = cur.fetchall()

            if len(fetched) < 100:
                cur.execute("SELECT name, alt_name, description FROM event_names_enriched WHERE name LIKE ? LIMIT 0,?",
                            ('%'+name_filter+'%', 100-len(fetched),))
                fetched = fetched + cur.fetchall()



        # res = [item for sublist in cur.fetchall() for item in sublist]  # flatten result, faster than map(lambda x: x[0],l)
        res = map(lambda x: {'id': x[0], 'name': x[1], 'desc': x[2]}, fetched)

        con.close()

        return res