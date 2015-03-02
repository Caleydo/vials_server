import os
import sys
import time
import glob
import json
import shelve

from misopy.index_gff import index_gff
from misopy.settings import Settings
from misopy.miso import get_main_logger, compute_all_genes_psi

def build_project(output_dir, config_file):
    Settings.load("miso_settings.txt")

    with open(config_file) as config_in:
        config = json.load(config_in)

    # naively prevent overwriting existing projects
    if os.path.isdir(output_dir):
        raise Exception("Directory already exists. Please choose a new location.")
    else:
        os.makedirs(output_dir)

    with open(os.path.join(output_dir, "config.json"), "w") as config_out:
        json.dump(config, config_out)

    gff_dir = os.path.join(output_dir, "indexed_gff")
    index_gff(config["gff_file"],
              gff_dir,
              compress_id=False,
              use_json=True,
              rel_paths=True)

    for sample, bam_file in config["samples_to_bam"].iteritems():
        # directory for each sample contains symlinked bamfile and miso output
        sample_dir = os.path.join(output_dir, "samples", sample)
        os.makedirs(sample_dir)
        os.symlink(bam_file, os.path.join(sample_dir, os.path.basename(bam_file)))

        # run miso on bamfile
        main_logger = get_main_logger(os.path.join(sample_dir, "logs"))
        compute_all_genes_psi(gff_dir, bam_file, config["read_len"],
                              os.path.join(sample_dir, "miso-data"), main_logger)

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-o", "--outdir", dest="output_dir", nargs=1, default=None,
                      help="Output directory for new project.")
    parser.add_option("-c", "--config", dest="config_file", nargs=1, default=None,
                      help="Configuration file for new project.")
    (options, args) = parser.parse_args()

    if not (options.output_dir and options.config_file):
        print "Must specify output directory with --outdir [-o] and config file with --config [-c]."
    else:
        build_project(options.output_dir, options.config_file)


if __name__ == '__main__':
    main()
