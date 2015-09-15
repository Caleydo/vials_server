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

        res = cur.fetchall()

        con.close()

        return {'x':res}
