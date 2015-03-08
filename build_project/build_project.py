import os
import sys
import time
import glob
import json
import shelve

from misopy.index_gff import index_gff
from misopy.parse_gene import parseGene
from misopy.settings import Settings
from misopy.miso import get_main_logger, compute_all_genes_psi

def gene_to_dict(gene, filename):
    tx_start, tx_end, exon_starts, exon_ends, gene_obj, \
       mRNAs, strand, chromID = parseGene(filename, gene);
    exons = sorted(set([tuple(exon) for RNA in mRNAs for exon in RNA]), key=lambda x: x[0])
    mRNAs_span = mRNAs
    mRNAs = [sorted([exons.index(tuple(exon)) for exon in RNA]) for RNA in mRNAs]

    return {'chromID': chromID, 'tx_start': tx_start, 'tx_end': tx_end,
            'exons': exons, 'mRNAs': mRNAs, 'strand': strand}

def build_project(output_dir, config_file, indexed_gff):
    gff_dir = os.path.join(output_dir, "indexed_gff")
    genes_filename = os.path.join(gff_dir, "genes_to_filenames.json")
    genes_info_file = os.path.join(gff_dir, "genes_info.json")

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

    index_gff(os.path.join(config["host_base"], config["gff_file"]),
              gff_dir,
              compress_id=False,
              use_json=True,
              rel_paths=True)

    for sample, sample_info in config["samples"].iteritems():
        # directory for each sample contains symlinked bamfile and miso output
        sample_dir = os.path.join(output_dir, "samples", sample)
        os.makedirs(sample_dir)

        bam_file_host_path = os.path.join(config["host_base"], sample_info["rel_path"])

        # run miso on bamfile
        miso_dir = os.path.join(sample_dir, "miso-data")
        main_logger = get_main_logger(os.path.join(sample_dir, "logs"))
        compute_all_genes_psi(gff_dir, bam_file_host_path, sample_info["read_len"],
                              miso_dir, main_logger)

    with open(genes_filename) as genes_fp:
        genes = json.load(genes_fp)

    genes_info = {gene: gene_to_dict(gene, os.path.join(gff_dir, filename))
                            for gene, filename in genes.iteritems()}

    with open(genes_info_file, "w") as genes_info_out:
        json.dump(genes_info, genes_info_out)

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
