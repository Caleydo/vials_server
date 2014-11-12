from flask.ext.cors import CORS
# from flask.ext.restful.utils.cors import crossdomain

__author__ = 'Hendrik Strobelt'

import pysam
from flask import Flask, request
from flask.ext.restplus import Api, Resource
import configparser
from misopy import index_gff, parse_gene
from misopy.sashimi_plot.plot_utils.plot_gene import readsToWiggle_pysam
import os
import shelve
import json

config = configparser.ConfigParser()
config.read('config.ini')
samples = json.loads(config['SERVER']['samples'])
bam_dir = config['SERVER']['bam_dir']
sample_files = {s: os.path.join(bam_dir, s + ".sorted.bam") for s in samples}
miso_dir = config['SERVER']['miso_dir']
gff_file = config['SERVER']['gff_file']

app = Flask(__name__)
api = Api(app)
cors = CORS(app)

parser = api.parser()
parser.add_argument('chromID', type=str, help='chromosome ID')
parser.add_argument('pos', type=int, help='gene position')
parser.add_argument('baseWidth', type=float, help='number of bases')
ns = api.namespace('bam', description='BAM operations')

@ns.route("/pileup")
class BAMInfo(Resource):
    def get(self):
        chromID = '1'
        pos=0
        base_width = 128

        args = parser.parse_args()
        chromID = args['chromID'] or chromID
        pos = args['pos'] or pos
        base_width = args['baseWidth'] or base_width
        print args

        sample_data = {}
        for sample, bam_file in sample_files.iteritems():
            sample_positions = []
            samfile = pysam.Samfile(bam_file, "rb")
            print samples
            subset_reads = samfile.fetch(reference=chromID, start=pos,end=pos+base_width)
            wiggle, jxn_reads = readsToWiggle_pysam(subset_reads, pos, pos+base_width)
            wiggle = wiggle.tolist()

            sample_jxns = []
            for jxn_range_str, jxn_count in jxn_reads.iteritems():
                jxn_range = map(int, jxn_range_str.split(':'))
                sample_jxns.append((jxn_range, jxn_count))

            for pileupcolumn in samfile.pileup(chromID, pos, pos+base_width):
                seqs =[]
                for pileupread in pileupcolumn.pileups:
                    seqs.append({'name':pileupread.alignment.qname,'val':pileupread.alignment.seq[pileupread.qpos]})
                if ((pileupcolumn.pos>=pos) & (pileupcolumn.pos-base_width<pos)):
                    sample_positions.append({'pos': pileupcolumn.pos, 'no': pileupcolumn.n, 'seq': seqs, 'wiggle': wiggle.pop(0)})

            samfile.close()
            sample_data[sample] = {'positions': sample_positions, 'jxns': sample_jxns}

        return sample_data

    # def post(self):
    #     api.abort(403)

@ns.route("/header")
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

@ns.route("/genes")
class BAMGenesInfo(Resource):
    def get(self):
        pickle_dir = os.path.dirname(gff_file)
        genes_filename = os.path.join(pickle_dir, "genes_to_filenames.shelve")
        
        genes = dict(shelve.open(genes_filename))

        def gene_to_dict(gene, filename):
            tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
               mRNAs, strand, chromID = parse_gene.parseGene(filename, gene);

            exons = sorted(set([tuple(exon) for RNA in mRNAs for exon in RNA]), key=lambda x: x[0])
            mRNAs = [sorted([exons.index(tuple(exon)) for exon in RNA]) for RNA in mRNAs]

            psis = {}
            for sample in samples:
                sample_psis = []
                for line in open(os.path.join(miso_dir, sample, chromID, gene + ".miso")):
                    if not line.startswith("#") and not line.startswith("sampled"):
                        psi, logodds = line.strip().split("\t")
                        sample_psis.append(float(psi.split(",")[0]))
                psis[sample] = float(sum(sample_psis))/len(sample_psis)

            return {'chromID': chromID, 'tx_start': tx_start, 'tx_end': tx_end,
                    'exons': exons, 'mRNAs': mRNAs, 'psis': psis, 'strand': strand}

        return {gene: gene_to_dict(gene, filename) for gene, filename in genes.iteritems()}

if __name__ == '__main__':
    # app.debug1 = True
    app.run()
