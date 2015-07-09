from flask.ext.cors import CORS
# from flask.ext.restful.utils.cors import crossdomain
import pickle

__author__ = 'Hendrik Strobelt'

import pysam
from flask.ext.restplus import Resource
from caleydo_server.apiutil import create_api

import configparser
from intervaltree import IntervalTree
# from misopy import index_gff, parse_gene
from misopy.sashimi_plot.plot_utils.plot_gene import readsToWiggle_pysam
from misopy.parse_gene import parseGene
import os
# import shelve
import json
import pprint

from handlers.bodymap import BodyMapHandler
from handlers.tcga import TCGAHandler

# create the api application
app, api = create_api(__name__, version='1.0', title='Caleydo AltSplice API', description='splicing data handler')


class Helpers:
    @staticmethod
    def get_subdirs(path):
        res = []
        for x in os.listdir(path):
            subdir = os.path.join(path, x)
            if os.path.isdir(subdir):
                res.append(subdir)
        return res

    @staticmethod
    def project_info(path):
        res = {}
        with open(config_file_in_project(path)) as jsonData:
            res = json.load(jsonData)
            jsonData.close()
        res["dir"] = path
        return os.path.basename(path), res  # TODO: hack until name in config file

    @staticmethod
    def ref_genome_info(path):
        res = {}
        with open(config_file_in_ref_genome(path)) as jsonData:
            res = json.load(jsonData)
            jsonData.close()
        res["dir"] = path
        return os.path.basename(path), res  # TODO: hack until name in config file



    @staticmethod
    def get_data_handler(data_type, data_root):
        if data_type == "BodyMap":
            return BodyMapHandler(data_root)
        elif data_type == "TCGA":
            return TCGAHandler(data_root)





exec_root_dir = os.path.dirname(os.path.realpath(__file__))


# define root directory for all AltSplice importations
import caleydo_server.config
data_root = caleydo_server.config.view('vials_server').data_root
if not os.path.exists(data_root):
    exit(1)

# define paths to subdirectories & files
projects_root_dir = os.path.join(data_root, "projects")
ref_genomes_root_dir = os.path.join(data_root, "ref_genomes")
samples_dir_in_project = lambda x: os.path.join(x, "samples")
config_file_in_project = lambda x: os.path.join(x, "config.json")

bam_file_in_project_and_sample = lambda project, sample_info:  os.path.join(project["host_base"], sample_info["rel_path"])

gene_info_in_ref_genome = lambda x: os.path.join(x, "genes_to_filenames.json")
config_file_in_ref_genome = lambda x: os.path.join(x, "meta.json")


# create a dictionary of projects
all_projects = dict(map(Helpers.project_info, Helpers.get_subdirs(projects_root_dir)))
currentProject = all_projects[all_projects.keys()[0]]

# dict for all reference genomes
all_ref_genomes = dict(map(Helpers.ref_genome_info, Helpers.get_subdirs(ref_genomes_root_dir)))
current_ref_genome = all_ref_genomes[all_ref_genomes.keys()[0]]

parser = api.parser()
# parser.add_argument('chromID', type=str, help='chromosome ID')
parser.add_argument('geneName', type=str, help='gene name')
parser.add_argument('pos', type=int, help='gene position')
parser.add_argument('baseWidth', type=float, help='number of bases')
parser.add_argument('projectID', type=str, help='projectID')
parser.add_argument('noReads', type=bool, help='projectID')

# project_dir = config['server']['project_dir']
#
# with open(os.path.join(project_dir, "config.json")) as config_json:
# project_config = json.load(config_json)
#
# samples = project_config["samples"].keys()
# samples_dir = os.path.join(project_dir, "samples")
# pickle_dir = os.path.join(project_dir, "indexed_gff")
# genes_filename = os.path.join(pickle_dir, "genes_to_filenames.json")
# genes_info_file = os.path.join(pickle_dir, "genes_info.json")
# with open(genes_info_file) as genes_info_in:
#     genes_info = json.load(genes_info_in)
#
# # create the api application
# app, api = createAPI(__name__, version='1.0', title='Caleydo Web BAM API', description='BAM operations')
#
# cors = CORS(app)
#
# parser = api.parser()
# # parser.add_argument('chromID', type=str, help='chromosome ID')
# parser.add_argument('geneName', type=str, help='gene name')
# parser.add_argument('pos', type=int, help='gene position')
# parser.add_argument('baseWidth', type=float, help='number of bases')
#
#
# @api.route("/pileup")
# class BAMInfo(Resource):
#     def get(self):
#         args = parser.parse_args()
#         # chromID = args['chromID'] or chromID
#         geneName = args['geneName']  # or geneName
#
#         try:
#             with open(genes_info_file) as genes_info_in:
#                 genes_info = json.load(genes_info_in)
#             gene_info = genes_info[geneName]
#         except Exception as e:
#             print "Gene {0} not found.".format(geneName)
#             raise e
#
#         chromID = gene_info["chromID"]
#         pos = args['pos'] or gene_info["tx_start"]
#         base_width = args['baseWidth'] or gene_info["tx_end"] - gene_info["tx_start"] + 1
#
#         # find exons present in view range
#         y = [tuple(list(exon) + [i]) for i, exon in enumerate(gene_info["exons"])]
#         exon_tree = IntervalTree.from_tuples(y)
#         full_overlaps = IntervalTree(exon_tree[pos:pos + base_width + 1])
#         curExons = []
#         exon_group_ids = {}
#         visited = set()
#         while full_overlaps:
#             # get next exon in view
#             cur_exon = full_overlaps.pop()
#             exon_group = [cur_exon]
#             # find exons that overlap this one
#             cur_overlaps = full_overlaps[cur_exon.begin:cur_exon.end]
#             exon_start, exon_end = cur_exon.begin, cur_exon.end
#             while cur_overlaps:
#                 # mark overlapping exons as visited and remove from cur_overlaps
#                 overlap_exon = cur_overlaps.pop()
#                 exon_group.append(overlap_exon)
#                 visited.add(overlap_exon)
#                 # update group range
#                 exon_start = min(exon_start, overlap_exon.begin)
#                 exon_end = max(exon_end, overlap_exon.end)
#                 # widen group to include overlapping exons not seen before
#                 new_overlaps = [exon for exon in full_overlaps[overlap_exon.begin:overlap_exon.end]
#                                 if exon not in visited]
#                 cur_overlaps = cur_overlaps.union(new_overlaps)
#             # remove full group of exons
#             full_overlaps.remove_overlap(exon_start, exon_end)
#             curExons.append([exon_start, exon_end])
#
#             # map combined exon to single id
#             for x in exon_group:
#                 exon_group_ids[x.data] = len(curExons) - 1
#
#         # this could fail if two exons in the same group but don't themselves overlap
#         curRNAs = []
#         for RNA in gene_info["mRNAs"]:
#             if set(exon_group_ids.keys()).intersection(RNA):
#                 curRNAs.append([exon_group_ids[i] for i in RNA])
#
#         data = {'geneInfo': {
#             'curExons': curExons,
#             'curRNAs': curRNAs,
#             'geneSpan': [gene_info['tx_start'], gene_info['tx_end']]
#         },
#                 'samples': {}}
#
#         max_wiggle = 0
#         for sample, sample_info in project_config["samples"].iteritems():
#             sample_dir = os.path.join(project_dir, "samples", sample)
#             sample_positions = []
#
#             # pysam won't take a unicode string, probably should fix
#             chromID = chromID.encode('ascii', 'ignore')
#             bam_file_vagrant_path = os.path.join(project_config["vagrant_base"], sample_info["rel_path"])
#             samfile = pysam.Samfile(bam_file_vagrant_path, "rb")
#             subset_reads = samfile.fetch(reference=chromID, start=pos, end=pos + base_width + 1)
#             wiggle, jxn_reads = readsToWiggle_pysam(subset_reads, pos, pos + base_width + 1)
#             wiggle = (wiggle / sample_info["coverage"]).tolist()
#             max_wiggle = max(max(wiggle), max_wiggle)
#
#             sample_jxns = []
#             for jxn_range_str, jxn_count in jxn_reads.iteritems():
#                 jxn_range = map(int, jxn_range_str.split(':'))
#                 if jxn_range[0] >= pos and jxn_range[1] <= pos + base_width:
#                     sample_jxns.append((jxn_range, jxn_count))
#
#             for pileupcolumn in samfile.pileup(chromID, pos, pos + base_width):
#                 seqs = []
#                 # for pileupread in pileupcolumn.pileups:
#                 # seqs.append({'name':pileupread.alignment.qname,'val':pileupread.alignment.seq[pileupread.qpos]})
#                 if ((pileupcolumn.pos >= pos) & (pileupcolumn.pos - base_width < pos)):
#                     sample_positions.append(
#                         {'pos': pileupcolumn.pos, 'no': pileupcolumn.n, 'seq': seqs, 'wiggle': wiggle.pop(0)})
#
#             # amended from miso's processing code
#             sample_psis = []
#             try:
#                 with open(genes_filename) as genes_in:
#                     genes_to_filenames = json.load(genes_in)
#                     gene_filename = os.path.splitext(genes_to_filenames[geneName])[0] + ".miso"
#                     miso_path = os.path.join(sample_dir, "miso-data", gene_filename)
#                     with open(miso_path) as miso_in:
#                         for line in miso_in.xreadlines():
#                             if not line.startswith("#") and not line.startswith("sampled"):
#                                 psi, logodds = line.strip().split("\t")
#                                 sample_psis.append(map(float, psi.split(",")))
#                         sample_psis = [sum(psis) / len(psis) for psis in zip(*sample_psis)]
#             except:
#                 print "Psis could not be found for sample {0}.".format(sample)
#
#             samfile.close()
#             data["samples"][sample] = {'positions': sample_positions, 'jxns': sample_jxns, 'psis': sample_psis}
#
#         # normalize wiggles
#         for sample in samples:
#             for d in data["samples"][sample]["positions"]:
#                 d["wiggle"] = d["wiggle"] / float(max_wiggle)
#
#         return data
#
#         # def post(self):
#         # api.abort(403)
#
#
# @api.route("/header")
# class BAMHeaderInfo(Resource):
#     def get(self):
#         headers = {}
#         for sample, sample_info in project_config["samples"].iteritems():
#             bam_file_vagrant_path = os.path.join(project_config["vagrant_base"], sample_info["rel_path"])
#             samfile = pysam.Samfile(bam_file_vagrant_path, "rb")
#             headers[bam_file] = samfile.header
#             samfile.close()
#         return headers
#
#         # def post(self):
#         # api.abort(403)
#
#


@api.route("/gene")
class GeneInfo(Resource):
    def get(self):
        args = parser.parse_args()

        # parse gene name (necessary)
        geneName = args["geneName"]
        if not geneName:
            return {}

        # parse project name (should)
        project = currentProject
        if args["projectID"]:
            project = all_projects[args["projectID"]]


        add_reads = True
        if args["noReads"]:
            add_reads = False

        isoform_measured = []
        jxns = []
        all_sapmple_infos = {}
        all_jxns_starts = []
        all_jxns_ends = []
        all_exons = {}
        all_isoforms = {}
        sample_reads = []


        chromID, strand, tx_end, tx_start = (0, "+", 0, 100)

        for datagroup in project["data"]:
            data_type = datagroup['data_type']
            handler = Helpers.get_data_handler(data_type, data_root)

            # --------------------------
            # reference data
            # --------------------------
            chromID,  strand, tx_end, tx_start, exons, isoforms, merged_ranges = handler.generate_meta_info(geneName, project, datagroup)
            all_exons.update(exons)
            all_isoforms.update(isoforms)

            # --------------------------
            # project data
            # --------------------------
            handler.read_data(geneName, add_reads, all_jxns_ends, all_jxns_starts, all_sapmple_infos,  datagroup,
                              isoform_measured, jxns, sample_reads, project, tx_start, tx_end, chromID)

        # settify and sort
        all_jxns_starts = sorted(list(set(all_jxns_starts)))
        all_jxns_ends = sorted(list(set(all_jxns_ends)))

        # ex_orig = []
        # ex_new = []
        # ex_infos = []

        # for isopart in isodesc.strip("'").split("_"):
        #     ex_orig.append(isopart)
        #     ex_infos.append(exonMap[isopart])
        #     ex_new.append(exonExonMap[isopart]["id"])

        # isoform_measured.append({"orig":ex_orig, "new":ex_new, "info": ex_infos, "mean": mean, "sample":sample})

        chromID = chromID.replace("chr","")

        theGene = {'chromID': chromID,
                   'start': tx_start,
                   'end': tx_end,
                   'strand': strand,
                   'exons': all_exons,
                   'isoforms': all_isoforms,
                   "name": geneName,
                   "merged_ranges": merged_ranges
                   }

        theData = {"jxns":
                       {"all_starts": all_jxns_starts,
                        "all_ends": all_jxns_ends,
                        "weights": jxns},
                   "isoforms":isoform_measured,
                   "reads": sample_reads,
                   "data_type": datagroup['data_type'],# TODO: allow projects with different data types
                    "isoform_unit": datagroup['isoform_unit']
                   }

        return {'gene': theGene,
                'measures': theData,
                'samples': all_sapmple_infos
                }


#  Test: ENSG00000168769
# test2: ENSG00000185345
# tP53: ENSG00000141510 // 1649196	TP53




@api.route("/genes")
class GeneOverview(Resource):
    def get(self):
        args = parser.parse_args()

        project = currentProject
        # ref_genome_dir = current_ref_genome["dir"]
        if args["projectID"]:
            project = all_projects[args["projectID"]]

        geneIDs = []
        for datagroup in project["data"]:
            data_type = datagroup['data_type']
            handler = Helpers.get_data_handler(data_type, data_root)
            genes = handler.get_genes(project, datagroup)
            geneIDs.extend(genes)

        # geneIDs = Helpers.load_genome_mapping(project).keys()

        return geneIDs


@api.route("/projects")
class ProjectInfos(Resource):
    def get(self):
        # allProjects = dict(map(Helpers.project_info, Helpers.get_subdirs()))
        # anyway run by each API call
        projectinfo = {project: info for project, info in all_projects.iteritems()}
        return projectinfo


if __name__ == '__main__':
    app.run()


def create(*args, **kwargs):
    return app

