from flask.ext.cors import CORS
# from flask.ext.restful.utils.cors import crossdomain

__author__ = 'Hendrik Strobelt'

import pysam
from flask import Flask, request
from flask.ext.restplus import Api, Resource
import configparser

config = configparser.ConfigParser()
config.read('config.ini')
bamFileLocation = config['SERVER']['bamFile']


app = Flask(__name__)
api = Api(app)
cors = CORS(app)

parser = api.parser()
parser.add_argument('pos', type=int, help='gene position')
ns = api.namespace('bam', description='BAM operations')



@ns.route("/pileup")
class BAMInfo(Resource):
    def get(self):
        pos=8033880

        args = parser.parse_args()
        pos = args['pos'] or pos
        print args
        positions =[]
        samfile = pysam.Samfile( bamFileLocation, "rb")

        for pileupcolumn in samfile.pileup('11', pos, pos+100):
            seqs =[]
            for pileupread in pileupcolumn.pileups:
                seqs.append({'name':pileupread.alignment.qname,'val':pileupread.alignment.seq[pileupread.qpos]})
            if ((pileupcolumn.pos>=pos) & (pileupcolumn.pos-100<=pos)):
                positions.append({'pos': pileupcolumn.pos, 'no': pileupcolumn.n, 'seq': seqs})
        samfile.close()
        return positions

    # def post(self):
    #     api.abort(403)

@ns.route("/header")
class BAMHeaderInfo(Resource):
    def get(self):
        samfile = pysam.Samfile( bamFileLocation, "rb" )
        h= samfile.header
        samfile.close()
        return h

    # def post(self):
    #     api.abort(403)

if __name__ == '__main__':
    app.run()
