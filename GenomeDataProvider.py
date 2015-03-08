from flask.ext.cors import CORS
# from flask.ext.restful.utils.cors import crossdomain

__author__ = 'Hendrik Strobelt'

import pysam
from flask.ext.restplus import Resource
from caleydo.apiutil import createAPI

import configparser
from intervaltree import IntervalTree
from misopy import index_gff, parse_gene
from misopy.sashimi_plot.plot_utils.plot_gene import readsToWiggle_pysam
import os
# import shelve
import json

config = configparser.ConfigParser()
dn = os.path.dirname(os.path.realpath(__file__))
config.read(dn + '/config.ini')

project_dir = config['SERVER']['project_dir']
with open(os.path.join(project_dir, "config.json")) as config_json:
    project_config = json.load(config_json)

samples = project_config["samples"].keys()
samples_dir = os.path.join(project_dir, "samples")
pickle_dir = os.path.join(project_dir, "indexed_gff")
genes_filename = os.path.join(pickle_dir, "genes_to_filenames.json")
genes_info_file = os.path.join(pickle_dir, "genes_info.json")
with open(genes_info_file) as genes_info_in:
    genes_info = json.load(genes_info_in)

# create the api application
app, api = createAPI(__name__, version='1.0', title='Caleydo Web BAM API', description='BAM operations')

cors = CORS(app)

parser = api.parser()
# parser.add_argument('chromID', type=str, help='chromosome ID')
parser.add_argument('geneName', type=str, help='gene name')
parser.add_argument('pos', type=int, help='gene position')
parser.add_argument('baseWidth', type=float, help='number of bases')

@api.route("/pileup")
class BAMInfo(Resource):
    def get(self):
        args = parser.parse_args()
        # chromID = args['chromID'] or chromID
        geneName = args['geneName']# or geneName

        try:
            with open(genes_info_file) as genes_info_in:
                genes_info = json.load(genes_info_in)
            gene_info = genes_info[geneName]
        except Exception as e:
            print "Gene {0} not found.".format(geneName)
            raise e


        chromID = gene_info["chromID"]
        pos = args['pos'] or gene_info["tx_start"]
        base_width = args['baseWidth'] or gene_info["tx_end"]-gene_info["tx_start"]+1

        # find exons present in view range
        y = [tuple(list(exon) + [i]) for i,exon in enumerate(gene_info["exons"])]
        exon_tree = IntervalTree.from_tuples(y)
        full_overlaps = IntervalTree(exon_tree[pos:pos+base_width+1])
        curExons = []
        exon_group_ids = {}
        visited = set()
        while full_overlaps:
            # get next exon in view
            cur_exon = full_overlaps.pop()
            exon_group = [cur_exon]
            # find exons that overlap this one
            cur_overlaps = full_overlaps[cur_exon.begin:cur_exon.end]
            exon_start, exon_end = cur_exon.begin, cur_exon.end
            while cur_overlaps:
                # mark overlapping exons as visited and remove from cur_overlaps
                overlap_exon = cur_overlaps.pop()
                exon_group.append(overlap_exon)
                visited.add(overlap_exon)
                # update group range
                exon_start = min(exon_start, overlap_exon.begin)
                exon_end = max(exon_end, overlap_exon.end)
                # widen group to include overlapping exons not seen before
                new_overlaps = [exon for exon in full_overlaps[overlap_exon.begin:overlap_exon.end]
                                     if  exon not in visited]
                cur_overlaps = cur_overlaps.union(new_overlaps)
            # remove full group of exons
            full_overlaps.remove_overlap(exon_start, exon_end)            
            curExons.append([exon_start, exon_end])

            # map combined exon to single id
            for x in exon_group:
                exon_group_ids[x.data] = len(curExons) - 1

        # this could fail if two exons in the same group but don't themselves overlap
        curRNAs = []
        for RNA in gene_info["mRNAs"]:
            if set(exon_group_ids.keys()).intersection(RNA):
                curRNAs.append([exon_group_ids[i] for i in RNA])

        data = {'geneInfo': {
                                'curExons': curExons,
                                'curRNAs': curRNAs,
                                'geneSpan': [gene_info['tx_start'], gene_info['tx_end']]
                            },
                'samples': {}}

        max_wiggle = 0
        for sample, sample_info in project_config["samples"].iteritems():
            sample_dir = os.path.join(project_dir, "samples", sample)
            sample_positions = []

            # pysam won't take a unicode string, probably should fix
            chromID = chromID.encode('ascii','ignore')
            bam_file_vagrant_path = os.path.join(project_config["vagrant_base"], sample_info["rel_path"])
            samfile = pysam.Samfile(bam_file_vagrant_path, "rb")
            subset_reads = samfile.fetch(reference=chromID, start=pos, end=pos+base_width+1)
            wiggle, jxn_reads = readsToWiggle_pysam(subset_reads, pos, pos+base_width+1)
            wiggle = (wiggle / sample_info["coverage"]).tolist()
            max_wiggle = max(max(wiggle), max_wiggle)

            sample_jxns = []
            for jxn_range_str, jxn_count in jxn_reads.iteritems():
                jxn_range = map(int, jxn_range_str.split(':'))
                if jxn_range[0] >= pos and jxn_range[1] <= pos+base_width:
                    sample_jxns.append((jxn_range, jxn_count))

            for pileupcolumn in samfile.pileup(chromID, pos, pos+base_width):
                seqs =[]
                # for pileupread in pileupcolumn.pileups:
                #     seqs.append({'name':pileupread.alignment.qname,'val':pileupread.alignment.seq[pileupread.qpos]})
                if ((pileupcolumn.pos>=pos) & (pileupcolumn.pos-base_width<pos)):
                    sample_positions.append({'pos': pileupcolumn.pos, 'no': pileupcolumn.n, 'seq': seqs, 'wiggle': wiggle.pop(0)})

            # amended from miso's processing code
            sample_psis = []
            try:
                with open(genes_filename) as genes_in:
                    genes_to_filenames = json.load(genes_in)
                    gene_filename = os.path.splitext(genes_to_filenames[geneName])[0] + ".miso"
                    miso_path = os.path.join(sample_dir, "miso-data", gene_filename)
                    with open(miso_path) as miso_in:
                        for line in miso_in.xreadlines():
                            if not line.startswith("#") and not line.startswith("sampled"):
                                psi, logodds = line.strip().split("\t")
                                sample_psis.append(map(float, psi.split(",")))
                        sample_psis = [sum(psis)/len(psis) for psis in zip(*sample_psis)]
            except:
                print "Psis could not be found for sample {0}.".format(sample)

            samfile.close()
            data["samples"][sample] = {'positions': sample_positions, 'jxns': sample_jxns, 'psis': sample_psis}

        # normalize wiggles
        for sample in samples:
            for d in data["samples"][sample]["positions"]:
                d["wiggle"] = d["wiggle"] / float(max_wiggle)

        return data

    # def post(self):
    #     api.abort(403)

@api.route("/header")
class BAMHeaderInfo(Resource):
    def get(self):
        headers = {}
        for sample, sample_info in project_config["samples"].iteritems():
            bam_file_vagrant_path = os.path.join(project_config["vagrant_base"], sample_info["rel_path"])
            samfile = pysam.Samfile(bam_file_vagrant_path, "rb")
            headers[bam_file] = samfile.header
            samfile.close()
        return headers

    # def post(self):
    #     api.abort(403)

@api.route("/genes")
class BAMGenesInfo(Resource):
    def get(self):
        geneinfo = {gene: info for gene, info in genes_info.iteritems()}
        # geneinfo = {gene: info for gene, info in genes_info.iteritems() if info["chromID"] == "11"}
        return geneinfo


if __name__ == '__main__':
    app.run()

def create():
  return app

