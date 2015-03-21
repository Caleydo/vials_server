import os
from flask import json
from intervaltree import IntervalTree
from misopy.parse_gene import parseGene

__author__ = 'hen'
class BodyMapHandler:

    def __init__(self, data_root):
        self.name = "BodyMapHandler"
        self.summary_file_in_project = lambda project, sample: os.path.join(project, "summary", "%s.miso_summary.json" % sample)
        self.data_type = "BodyMap"
        self.gene_info_in_ref_genome = lambda x: os.path.join(x, "genes_to_filenames.json")
        self.ref_genomes_root_dir = os.path.join(data_root, "ref_genomes")


    def load_genome_mapping(self, project):
        if project["ref_genome"] == "local":
            ref_genome_dir = os.path.join(project["dir"], "local_ref_genome", project["local_ref_genome"])
        else:
            ref_genome_dir = os.path.join(self.ref_genomes_root_dir, project["ref_genome"])

        with open(self.gene_info_in_ref_genome(ref_genome_dir)) as jsonData:
            loaded = json.load(jsonData)
            # clean absoulte path: TODO: should be solved by project creator
            res = {}
            if project["ref_genome"] == "local":
                for index, key in enumerate(loaded):
                    path = loaded[key]
                    res[key] = os.path.join(ref_genome_dir,
                                             os.path.basename(os.path.dirname(path)),
                                             os.path.basename(path))
            else:
                res = loaded

            jsonData.close()
        return res


    @staticmethod
    def merge_exon_ranges(exons):
        # maybe not the most efficient -- TODO: check for better method
        itree = IntervalTree()
        for _,exon in exons.iteritems():
            itree[exon["start"]:exon["end"]] = exon["id"]
        itree.split_overlaps()
        merged_ranges = []
        cur_min = -1
        cur_max = -1
        cur_names = []
        cur_exon_count = 0
        for x in sorted(list(itree)):
            if cur_max == -1:  # for init
                cur_min = x.begin
                cur_max = x.end
                cur_names.append(x.data)
            else:
                if x.begin <= cur_max + 1:  # if connected
                    cur_max = x.end
                    cur_names.append(x.data)
                else:  # if not connected
                    exid = "ex" + format(cur_exon_count, '05d')
                    merged_ranges.append({"start": cur_min, "end": cur_max, "names": list(set(cur_names)), "id": exid})

                    cur_exon_count += 1
                    cur_min = x.begin
                    cur_max = x.end
                    cur_names = [x.data]

        # don't forget the last one
        if cur_names.__len__() > 0:
            exid = "ex" + format(cur_exon_count, '05d')
            merged_ranges.append({"start": cur_min, "end": cur_max, "names": list(set(cur_names)), "id": exid})
        return merged_ranges


    def get_genes(self, project, datagroup):
        genes = self.load_genome_mapping(project)
        return genes;


    def generate_meta_info(self, geneName, project, datagroup):
        # map gene name to a file
        gene_file_mapping = self.load_genome_mapping(project)
        gene_file = gene_file_mapping[geneName]
        # parse gene info
        tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
        mRNAs, strand, chromID = parseGene(gene_file, geneName)
        # create exons
        # exons = []
        exonMap = {}
        for x in gene_obj.parts:
            ex_info = {"start": x.start, "end": x.end, "id": x.label}
            # exons.append(ex_info)
            exonMap[x.label] = ex_info

        # merge exons
        merged_ranges = self.merge_exon_ranges(exonMap)
        # create a mapping between merged exon version and new name
        exonExonMap = {}
        for newExon in merged_ranges:
            for oldExon in newExon["names"]:
                exonExonMap[oldExon] = newExon

        # add (reference) isoforms
        isoforms = {}
        self.isoformIDMap ={}
        for x in gene_obj.isoforms:
            related_exons = x.desc #map(lambda exonName: exonExonMap[exonName], x.desc)
            isoforms[x.label] = {"start": x.genomic_start, "end": x.genomic_end, "exons": related_exons, "id": x.label}
            self.isoformIDMap["_".join(related_exons)] = x.label

        self.exonExonMap = exonExonMap
        self.exonMap = exonMap

        return chromID, strand, tx_end, tx_start, exonMap, isoforms

    def read_data(self, geneName, all_jxns_ends, all_jxns_starts, all_sapmple_infos, datagroup, isoform_measured, jxns, sample_reads, project):
        # TODO: modularize -- now only BAM
        for sample, sample_info in datagroup['samples'].iteritems():
            all_sapmple_infos[sample] = {"id": sample, "type": self.data_type, "origin": sample_info}

            # MISO information
            summary_file_name = self.summary_file_in_project(project['dir'], sample)
            with open(summary_file_name, 'rb') as miso_file:
                miso_info = json.load(miso_file);
                mean_count = 0
                for isodesc in miso_info["isoforms"]:
                    isodesc = isodesc.strip(" '")
                    mean = miso_info['miso_posterior_mean'][mean_count]
                    mean_count += 1

                    # ex_orig = []
                    # ex_new = []
                    # ex_infos = []
                    # for isopart in isodesc.strip("'").split("_"):
                        # ex_orig.append(isopart)
                        # ex_infos.append(self.exonMap[isopart])
                        # ex_new.append(self.exonExonMap[isopart]["id"])

                    if isodesc in self.isoformIDMap:
                        isodesc = self.isoformIDMap[isodesc]
                    isoform_measured.append(
                        {"id": isodesc,
                         # "orig": ex_orig,
                         # "new": ex_new,
                         # "info": ex_infos,
                         "weight": mean,
                         "sample": sample
                        })

            # junction information
            # pysam won't take a unicode string, probably should fix
            sample_jxns_file_name = os.path.join(project['dir'], sample, geneName,
                                                 "jxns.json")  # TODO: create a head function for that
            with open(sample_jxns_file_name) as jxns_file:
                jxns_info = json.load(jxns_file)
                for range, value in jxns_info.iteritems():
                    start, end = range.split(":")
                    jxns.append({"start": start, "end": end, "weight": value, "sample": sample})
                    all_jxns_starts.append(start)
                    all_jxns_ends.append(end)

            sample_reads_file_name = os.path.join(project['dir'], sample, geneName,
                                                 "wiggles_10k.json")  # TODO: create a head function for that
            with open(sample_reads_file_name) as reads_file:
                all_samples = json.load(reads_file)
                sample_reads.append({"sample": sample,
                                     "min": min(all_samples),
                                     "max": max(all_samples),
                                     "weights": all_samples})

    # def fillGeneIDs