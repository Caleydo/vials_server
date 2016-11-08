import sqlite3 as sqlite
import os
from intervaltree import IntervalTree

__author__ = 'Hendrik Strobelt'


class GFFHandler:
    def __init__(self, rg_info):
        self.rg_info = rg_info

    @staticmethod
    def merge_exon_ranges(exons):
        # maybe not the most efficient -- TODO: check for better method
        itree = IntervalTree()
        for exon in exons:
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

    def get_gene_info(self, gene_id):
        con = sqlite.connect(self.rg_info['db_file'])
        con.row_factory = sqlite.Row

        cur = con.cursor()

        cur.execute('SELECT * FROM features WHERE ID= ?', (gene_id,))

        gene_info = cur.fetchone()

        cur.execute("SELECT parent AS isoformID,group_concat(ID, '_') AS exonNames,"
                    "group_concat(startpos, '_') AS exonStarts,"
                    "group_concat(endpos, '_') AS exonEnds "
                    "FROM features WHERE feature = 'exon' "
                    "AND parent IN (SELECT ID FROM features WHERE parent = ? AND feature ='mRNA') "
                    "GROUP BY parent", (gene_id,))

        isoform_info = cur.fetchall()

        con.close()
        #
        # rrr = {}

        allIso = {}
        allExons = []

        for isoi in isoform_info:
            rrr = {}
            for k in isoi.keys():
                rrr[k] = isoi[k]
            allIso[isoi['isoformID']] = rrr

            # fill all_exons:
            exons = isoi['exonNames'].split('_')
            starts = isoi['exonStarts'].split('_')
            ends = isoi['exonEnds'].split('_')

            if len(exons) == len(starts) == len(ends):
                for index, exon in enumerate(exons):
                    allExons.append({
                        "id": exon,
                        "start": int(starts[index]),
                        "end": int(ends[index])
                    })

        if len(allExons) > 0:
            merged_ranges = GFFHandler.merge_exon_ranges(allExons)
        else:
            merged_ranges = []

        gene = {
            'name': gene_info['ID'],
            'chromID': gene_info['chromID'],
            'end': gene_info['endpos'],
            'start': gene_info['startpos'],
            'strand': gene_info['strand'],
            'isoforms': allIso,
            'merged_ranges': merged_ranges
        }

        return gene
