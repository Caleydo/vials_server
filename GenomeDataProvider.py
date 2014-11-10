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
bamFiles = json.loads(config['SERVER']['bamFiles'])
gffFileLocation = config['SERVER']['gffFile']

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

        positions, jxns = [], []
        for bamFile in bamFiles:
            sample_positions, sample_jxns = [], []
            samfile = pysam.Samfile( bamFile, "rb")
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
            positions.append(sample_positions)
            jxns.append(sample_jxns)

        return {'samples': bamFiles, 'positions': positions, 'jxns': jxns}

    # def post(self):
    #     api.abort(403)

@ns.route("/wiggle")
class BAMWiggle(Resource):
    def get(self):
        chromID = '1'
        pos=0
        base_width = 128

        args = parser.parse_args()
        chromID = args['chromID'] or chromID
        pos = args['pos'] or pos
        base_width = args['baseWidth'] or base_width

        samfile = pysam.Samfile(bamFileLocation, 'rb')
        subset_reads = samfile.fetch(reference=chromID, start=pos,end=pos+base_width)
        print subset_reads
        wiggle, _ = readsToWiggle_pysam(subset_reads, pos, pos+base_width)

        return zip([pos + idx for idx in range(len(wiggle))], wiggle.tolist())

@ns.route("/header")
class BAMHeaderInfo(Resource):
    def get(self):
        headers = {}
        for bamFile in bamFiles:
            samfile = pysam.Samfile( bamFile, "rb" )
            headers[bamFile] = samfile.header
            samfile.close()
        return headers

    # def post(self):
    #     api.abort(403)

@ns.route("/genes")
class BAMGenesInfo(Resource):
    def get(self):
        pickle_dir = os.path.dirname(gffFileLocation)
        genes_filename = os.path.join(pickle_dir,
                              "genes_to_filenames.shelve")
        
        genes = dict(shelve.open(genes_filename))

        def gene_to_dict(gene, filename):
            tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
               mRNAs, strand, chromID = parse_gene.parseGene(filename, gene);

            exons = sorted(set([tuple(exon) for RNA in mRNAs for exon in RNA]), key=lambda x: x[0])
            mRNAs = [sorted([exons.index(tuple(exon)) for exon in RNA]) for RNA in mRNAs]

            return {'chromID': chromID, 'tx_start': tx_start, 'tx_end': tx_end,
                    'exons': exons, 'mRNAs': mRNAs, 'strand': strand}

        return {gene: gene_to_dict(gene, filename) for gene, filename in genes.iteritems()}

if __name__ == '__main__':
    app.debug = True
    app.run()
