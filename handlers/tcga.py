import os
from intervaltree import IntervalTree
from flask import json
from misopy.parse_gene import parseGene
from bodymap import BodyMapHandler

__author__ = 'hen'
class TCGAHandler:

    def __init__(self, data_root):
        self.name = "TCGAHandler"
        # self.summary_file_in_project = lambda project, sample: os.path.join(project, "summary", "%s.miso_summary.json" % sample)

        self.data_type = "TCGA"
        self.data_root = data_root
        self.gene_info_in_ref_genome = lambda x: os.path.join(x, "genes_to_filenames.json")
        self.ref_genomes_root_dir = os.path.join(data_root, "ref_genomes")
        self.exon_info_file_in_project = lambda proj: os.path.join(proj["dir"],"exon_info.json")
        self.jxns_in_project_and_sample = lambda proj, sample: os.path.join(proj["dir"],sample,"jxns.json")
        self.exons_in_project_and_sample = lambda proj, sample: os.path.join(proj["dir"],sample,"exons.json")
        self.isoforms_in_project_and_sample = lambda proj, sample: os.path.join(proj["dir"],sample,"isoforms.json")
        self.isoforms_info_file_in_project = lambda proj: os.path.join(proj["dir"], "isoform_info.json")

    def load_genome_mapping(self, project):
        if project["ref_genome"] == "local":
            ref_genome_dir = os.path.join(project["dir"], "local_ref_genome", project["local_ref_genome"])
        else:
            ref_genome_dir = os.path.join(self.ref_genomes_root_dir, project["ref_genome"])

        with open(self.gene_info_in_ref_genome(ref_genome_dir)) as jsonData:
            loaded = json.load(jsonData)
            # clean absoulte path: TODO: should be solved by project creator
            res = {}
            for index, key in enumerate(loaded):
                    path = loaded[key]
                    res[key] = os.path.join(ref_genome_dir,
                                             os.path.basename(os.path.dirname(path)),
                                             os.path.basename(path))

            jsonData.close()
        return res
    def get_genes(self, project, datagroup):
        return datagroup["geneIDs"]


    def unifyExonID(self, oldID):
        chromID, range, strand = oldID.split(":")
        chromID = chromID.replace("chr", "")
        range = range.replace("-", "_")
        return chromID+"_" + range + "_" + strand

    def generate_meta_info(self, geneName, project, datagroup):
        # map gene name to a file
        gene_file_mapping = self.load_genome_mapping(project)
        gene_file = gene_file_mapping[geneName]
        # parse gene info
        tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
        mRNAs, strand, chromID = parseGene(gene_file, geneName)

        exons = {}
        with open(self.exon_info_file_in_project(project)) as exon_inf_file:
            for _, ex in json.load(exon_inf_file).iteritems():
                if ex["gene_id"] == geneName:
                    newID = self.unifyExonID(ex["id"])
                    exons[newID] = {"start": ex["start"], "end": ex["end"], "id": newID, "name": ex["name"]}


        merged_ranges = BodyMapHandler.merge_exon_ranges(exons)




        isoforms = {}
        with open(self.isoforms_info_file_in_project(project)) as iso_info_file:
            for id, isoform in json.load(iso_info_file).iteritems():
                if isoform['gene_id'] == geneName:
                    isoform["exons"] = map(self.unifyExonID, isoform["exons"])
                    isoforms[id] = isoform

        return chromID, strand, tx_end, tx_start, exons, isoforms, merged_ranges

    def read_data(self, geneName, add_reads,  all_jxns_ends, all_jxns_starts, all_sapmple_infos, datagroup,
                  isoform_measured, jxns, sample_reads, project):

        for sample in datagroup['samples']:

            with open(self.jxns_in_project_and_sample(project,sample)) as jxns_file:
                allJunctions = json.load(jxns_file)

                for position, value in allJunctions.iteritems():
                    _, start, end, _ = position.split("_")
                    weight = {
                        "start": int(start),
                        "end": int(end),
                        "weight": float(value),
                        "sample": sample
                    }
                    all_jxns_starts.append(weight["start"])
                    all_jxns_ends.append(weight["end"])
                    jxns.append(weight)

            if add_reads:
                with open(self.exons_in_project_and_sample(project, sample)) as exons_file:
                    allExons = json.load(exons_file)
                    sample_exon_reads = {"sample": sample, "weights": []}

                    # tree to hold exon values
                    exon_tree = IntervalTree()
                    # tree to hold region values
                    region_tree = IntervalTree()

                    for exon_ID, weight in allExons.iteritems():
                        exon_start, exon_end = exon_ID.split("_")[1:3]
                        exon_tree[exon_start:exon_end] = float(weight)
                        region_tree[exon_start:exon_end] = 0
                    region_tree.split_overlaps()

                    for region in region_tree:
                        region_weight = 0
                        for exon in exon_tree[region.begin:region.end]:
                            region_weight += exon.data

                        sample_exon_reads["weights"].append({
                                "exon": exon_ID,
                                "start": region.begin,
                                "end": region.end,
                                "sample": sample,
                                "weight": region_weight
                            })

                    sample_exon_reads["weights"] = sorted(sample_exon_reads["weights"], key=lambda x: x["start"])
                    sample_reads.append(sample_exon_reads)


            # -- DATA ---
            # "uc003tqi.2": {
            #   "scaled_estimate": "0",
            #   "raw_count": "0.00"
            # }

            with open(self.isoforms_in_project_and_sample(project, sample)) as iso_file:
                all_iso = json.load(iso_file)

                for id, infos in all_iso.iteritems():
                    isoform_measured.append({
                        "id": id,
                        "weight": infos['scaled_estimate'],
                        "sample":sample
                    })





            all_sapmple_infos[sample]= {"id": sample, "type": "TCGA", "origin": {}}
