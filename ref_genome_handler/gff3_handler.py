import sqlite3 as sqlite
import os

__author__ = 'Hendrik Strobelt'


class GFFHandler:
    def __init__(self, rg_info):
        self.rg_info = rg_info

    def get_gene_info(self, gene_id):
        con = sqlite.connect(self.rg_info['db_file'])
        con.row_factory = sqlite.Row

        cur = con.cursor()

        cur.execute('SELECT * from features WHERE ID= ?', (gene_id,))

        gene_info = cur.fetchone()

        cur.execute("select parent as isoformID,group_concat(ID, '_') as exonNames,"
                    "group_concat(startpos, '_') as exonStarts,"
                    "group_concat(endpos, '_') as exonEnds "
                    "from features where feature = 'exon' "
                    "AND parent in (select ID from features where parent = ? AND feature ='mRNA') "
                    "GROUP BY parent", (gene_id,))

        isoform_info = cur.fetchall()

        con.close()
        #
        # rrr = {}

        allIso = {}

        for isoi in isoform_info:
            rrr = {}
            for k in isoi.keys():
                rrr[k] = isoi[k]
            allIso[isoi['isoformID']] = rrr







        gene = {
            'name': gene_info['ID'],
            'chromID': gene_info['chromID'],
            'end': gene_info['endpos'],
            'start': gene_info['startpos'],
            'strand': gene_info['strand'],
            'isoforms': allIso
        }

        return gene
