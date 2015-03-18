import os
from flask import json

__author__ = 'hen'
class BodyMapHandler:

    def __init__(self):
        self.name = "BodyMapHandler"
        self.summary_file_in_project = lambda project, sample: os.path.join(project, "summary", "%s.miso_summary.json" % sample)
        self.data_type = "BodyMap"

    def read_data(self, all_jxns_ends, all_jxns_starts, all_sapmple_infos, datagroup, exonExonMap, exonMap,
                  isoform_measured, jxns, project):
        # TODO: modularize -- now only BAM
        for sample, sample_info in datagroup['samples'].iteritems():
            all_sapmple_infos[sample] = {"id": sample, "type": self.data_type, "origin": sample_info}

            # MISO information
            summary_file_name = self.summary_file_in_project(project['dir'], sample)
            with open(summary_file_name, 'rb') as miso_file:
                miso_info = json.load(miso_file);
                mean_count = 0
                for isodesc in miso_info["isoforms"]:
                    mean = miso_info['miso_posterior_mean'][mean_count]
                    mean_count += 1

                    ex_orig = []
                    ex_new = []
                    ex_infos = []
                    for isopart in isodesc.strip("'").split("_"):
                        ex_orig.append(isopart)
                        ex_infos.append(exonMap[isopart])
                        ex_new.append(exonExonMap[isopart]["id"])

                    isoform_measured.append(
                        {"orig": ex_orig, "new": ex_new, "info": ex_infos, "mean": mean, "sample": sample})

            # junction information
            # pysam won't take a unicode string, probably should fix
            sample_jxns_file_name = os.path.join(project['dir'], sample,
                                                 "jxns.json")  # TODO: create a head function for that
            with open(sample_jxns_file_name) as jxns_file:
                jxns_info = json.load(jxns_file)
                jxns_file.close()
                for range, value in jxns_info.iteritems():
                    start, end = range.split(":")
                    jxns.append({"start": start, "end": end, "weight": value, "sample": sample})
                    all_jxns_starts.append(start)
                    all_jxns_ends.append(end)