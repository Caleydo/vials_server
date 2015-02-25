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
samples = json.loads(config['SERVER']['samples'])
bam_dir = config['SERVER']['bam_dir']
sample_files = {s: os.path.join(bam_dir, s + ".sorted.bam") for s in samples}
miso_dir = config['SERVER']['miso_dir']
gff_file = config['SERVER']['gff_file']
pickle_dir = os.path.join(os.path.dirname(gff_file), "indexed_gff")
genes_filename = os.path.join(pickle_dir, "genes_to_filenames.json")


# create the api application
app, api = createAPI(__name__, version='1.0', title='Caleydo Web BAM API', description='BAM operations')

cors = CORS(app)

parser = api.parser()
# parser.add_argument('chromID', type=str, help='chromosome ID')
parser.add_argument('geneName', type=str, help='gene name')
parser.add_argument('pos', type=int, help='gene position')
parser.add_argument('baseWidth', type=float, help='number of bases')

def gene_to_dict(gene, filename):
    tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
       mRNAs, strand, chromID = parse_gene.parseGene(filename, gene);
    exons = sorted(set([tuple(exon) for RNA in mRNAs for exon in RNA]), key=lambda x: x[0])
    mRNAs_span = mRNAs
    mRNAs = [sorted([exons.index(tuple(exon)) for exon in RNA]) for RNA in mRNAs]

    # psis = {}
    # for sample in samples:
    #     sample_psis = []
    #     for line in open(os.path.join(miso_dir, sample, chromID, gene + ".miso")):
    #         if not line.startswith("#") and not line.startswith("sampled"):
    #             psi, logodds = line.strip().split("\t")
    #             sample_psis.append(float(psi.split(",")[0]))
    #     psis[sample] = [sum(psis)/len(psis) for psis in zip(sample_psis)]

    return {'chromID': chromID, 'tx_start': tx_start, 'tx_end': tx_end,
            'exons': exons, 'mRNAs': mRNAs, 'strand': strand}

@api.route("/pileup")
class BAMInfo(Resource):
    def get(self):
        args = parser.parse_args()
        # chromID = args['chromID'] or chromID
        geneName = args['geneName']# or geneName

        try:
            with open(genes_filename) as genes_fp:
                genes = json.load(genes_fp)
            gene_info = gene_to_dict(geneName, genes[geneName])
        except Exception as e:
            print "Gene {0} not found.".format(geneName)
            raise e


        chromID = gene_info["chromID"]
        pos = args['pos'] or gene_info["tx_start"]
        base_width = args['baseWidth'] or gene_info["tx_end"]-gene_info["tx_start"]+1

        # find exons present in view range
        y = [tuple(list(exon) + [i]) for i,exon in enumerate(gene_info["exons"])]
        # x = [tuple(exon + [i]) for i, exon in enumerate(gene_info["exons"])]
        exon_tree = IntervalTree.from_tuples(y)
        curExonIDs = [exon.data for exon in exon_tree[pos:pos+base_width+1]]
        curRNAs = [RNA for RNA in gene_info["mRNAs"] if set(curExonIDs).intersection(RNA)]
        curExons = [gene_info["exons"][i] for i in curExonIDs]

        data = {'geneInfo': {
                                'curExons': curExons,
                                'curRNAs': curRNAs,
                                'geneSpan': [gene_info['tx_start'], gene_info['tx_end']]
                            },
                'samples': {}}
        for sample, bam_file in sample_files.iteritems():
            sample_positions = []

            # pysam won't take a unicode string, probably should fix
            chromID = chromID.encode('ascii','ignore')
            samfile = pysam.Samfile(bam_file, "rb")
            subset_reads = samfile.fetch(reference=chromID, start=pos, end=pos+base_width+1)
            wiggle, jxn_reads = readsToWiggle_pysam(subset_reads, pos, pos+base_width+1)
            wiggle = wiggle.tolist()

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

            samfile.close()
            data["samples"][sample] = {'positions': sample_positions, 'jxns': sample_jxns}

        return data

    # def post(self):
    #     api.abort(403)

@api.route("/header")
class BAMHeaderInfo(Resource):
    def get(self):
        headers = {}
        for bam_file in sample_files.values():
            samfile = pysam.Samfile(bam_file, "rb")
            headers[bam_file] = samfile.header
            samfile.close()
        return headers

    # def post(self):
    #     api.abort(403)

@api.route("/genes")
class BAMGenesInfo(Resource):
    def get(self):
        # pickle_dir = os.path.dirname(gff_file)
        # genes_filename = os.path.join(pickle_dir, "genes_to_filenames.json")
        #
        # genes = {}
        # with open(genes_filename) as jsonF:
        #     genes = json.load(jsonF)
        #
        # def gene_to_dict(gene, filename):
        #     tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
        #        mRNAs, strand, chromID = parse_gene.parseGene(filename, gene);
        #
        #     exons = sorted(set([tuple(exon) for RNA in mRNAs for exon in RNA]), key=lambda x: x[0])
        #     mRNAs = [sorted([exons.index(tuple(exon)) for exon in RNA]) for RNA in mRNAs]
        #
        #     psis = {}
        #     for sample in samples:
        #         sample_psis = []
        #         for line in open(os.path.join(miso_dir, sample, chromID, gene + ".miso")):
        #             if not line.startswith("#") and not line.startswith("sampled"):
        #                 psi, logodds = line.strip().split("\t")
        #                 sample_psis.append(float(psi.split(",")[0]))
        #         psi_avg = float(sum(sample_psis))/len(sample_psis)
        #         psis[sample] = [psi_avg, 1-psi_avg]
        #
        #     return {'chromID': chromID, 'tx_start': tx_start, 'tx_end': tx_end,
        #             'exons': exons, 'mRNAs': mRNAs, 'psis': psis, 'strand': strand}

        # return {gene: gene_to_dict(gene, filename) for gene, filename in genes.iteritems()}
        print "genes file", genes_filename

        gene_info = {}
        with open(genes_filename) as genes_fp:
            gene_info = json.load(genes_fp)

        res = {gene: gene_to_dict(gene, filename) for gene, filename in gene_info.iteritems()}

        return res
        # return gene_info.keys()


if __name__ == '__main__':
    # app.debug1 = True
    app.run()

def create():
  return app

