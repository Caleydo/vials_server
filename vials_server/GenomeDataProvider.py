from phovea_server.apiutil import create_api, Resource
import os
import json
from .data_handlers.bodymap import BodyMapHandler
from .data_handlers.tcga import TCGAHandler
# define root directory for all AltSplice importations
import phovea_server.config

__author__ = 'Hendrik Strobelt'

# create the api application
app, api = create_api(__name__, version='1.0', title='Caleydo AltSplice API', description='splicing data handler')


class Helpers:
  @staticmethod
  def get_subdirs(path):
    res = []
    for x in os.listdir(path):
      subdir = os.path.join(path, x)
      if os.path.isdir(subdir):
        res.append(subdir)
    return res

  @staticmethod
  def project_info(path):
    res = {}
    with open(config_file_in_project(path)) as jsonData:
      res = json.load(jsonData)
      jsonData.close()
    res["dir"] = path
    return os.path.basename(path), res  # TODO: hack until name in config file

  @staticmethod
  def ref_genome_info(path):
    res = {}
    with open(config_file_in_ref_genome(path)) as jsonData:
      res = json.load(jsonData)
      jsonData.close()
    res["dir"] = path
    return os.path.basename(path), res  # TODO: hack until name in config file

  @staticmethod
  def get_data_handler(data_type, data_root):
    if data_type == "BodyMap":
      return BodyMapHandler(data_root)
    elif data_type == "TCGA":
      return TCGAHandler(data_root)


exec_root_dir = os.path.dirname(os.path.realpath(__file__))

data_root = phovea_server.config.view('vials_server').data_root
if not os.path.exists(data_root):
  exit(1)

# define paths to subdirectories & files
projects_root_dir = os.path.join(data_root, "projects")
ref_genomes_root_dir = os.path.join(data_root, "ref_genomes")


def samples_dir_in_project(x):
  return os.path.join(x, "samples")


def config_file_in_project(x):
  return os.path.join(x, "config.json")


def bam_file_in_project_and_sample(project, sample_info):
  return os.path.join(project["host_base"], sample_info["rel_path"])


def gene_info_in_ref_genome(x):
  return os.path.join(x, "genes_to_filenames.json")


def config_file_in_ref_genome(x):
  return os.path.join(x, "meta.json")


# create a dictionary of projects
all_projects = dict(map(Helpers.project_info, Helpers.get_subdirs(projects_root_dir)))
currentProject = all_projects[all_projects.keys()[0]]

# dict for all reference genomes
all_ref_genomes = dict(map(Helpers.ref_genome_info, Helpers.get_subdirs(ref_genomes_root_dir)))
current_ref_genome = all_ref_genomes[all_ref_genomes.keys()[0]]

parser = api.parser()
# parser.add_argument('chromID', type=str, help='chromosome ID')
parser.add_argument('geneName', type=str, help='gene name')
parser.add_argument('pos', type=int, help='gene position')
parser.add_argument('baseWidth', type=float, help='number of bases')
parser.add_argument('projectID', type=str, help='projectID')
parser.add_argument('noReads', type=bool, help='projectID')


@api.route("/gene")
class GeneInfo(Resource):
  def get(self):
    args = parser.parse_args()

    # parse gene name (necessary)
    gene_name = args["geneName"]
    if not gene_name:
      return {}

    # parse project name (should)
    project = currentProject
    if args["projectID"]:
      project = all_projects[args["projectID"]]

    add_reads = True
    if args["noReads"]:
      add_reads = False

    isoform_measured = []
    jxns = []
    all_sapmple_infos = {}
    all_jxns_starts = []
    all_jxns_ends = []
    all_exons = {}
    all_isoforms = {}
    sample_reads = []

    chrom_id, strand, tx_end, tx_start = (0, "+", 0, 100)

    for datagroup in project["data"]:
      data_type = datagroup['data_type']
      handler = Helpers.get_data_handler(data_type, data_root)

      # --------------------------
      # reference data
      # --------------------------
      chrom_id, strand, tx_end, tx_start, exons, isoforms, merged_ranges = handler.generate_meta_info(gene_name,
                                                                                                      project,
                                                                                                      datagroup)
      all_exons.update(exons)
      all_isoforms.update(isoforms)

      # --------------------------
      # project data
      # --------------------------
      handler.read_data(gene_name, add_reads, all_jxns_ends, all_jxns_starts, all_sapmple_infos, datagroup,
                        isoform_measured, jxns, sample_reads, project, tx_start, tx_end, chrom_id)

    # settify and sort
    all_jxns_starts = sorted(list(set(all_jxns_starts)))
    all_jxns_ends = sorted(list(set(all_jxns_ends)))

    # ex_orig = []
    # ex_new = []
    # ex_infos = []

    # for isopart in isodesc.strip("'").split("_"):
    #     ex_orig.append(isopart)
    #     ex_infos.append(exonMap[isopart])
    #     ex_new.append(exonExonMap[isopart]["id"])

    # isoform_measured.append({"orig":ex_orig, "new":ex_new, "info": ex_infos, "mean": mean, "sample":sample})

    chrom_id = chrom_id.replace("chr", "")

    the_gene = {'chromID': chrom_id,
                'start': tx_start,
                'end': tx_end,
                'strand': strand,
                'exons': all_exons,
                'isoforms': all_isoforms,
                "name": gene_name,
                "merged_ranges": merged_ranges
                }

    the_data = {"jxns": {"all_starts": all_jxns_starts,
                         "all_ends": all_jxns_ends,
                         "weights": jxns},
                "isoforms": isoform_measured,
                "reads": sample_reads,
                "data_type": datagroup['data_type'],  # TODO: allow projects with different data types
                "isoform_unit": datagroup['isoform_unit']}

    return {'gene': the_gene,
            'measures': the_data,
            'samples': all_sapmple_infos
            }


# Test: ENSG00000168769
# test2: ENSG00000185345
# tP53: ENSG00000141510 // 1649196	TP53

@api.route("/genes")
class GeneOverview(Resource):
  def get(self):
    args = parser.parse_args()

    project = currentProject
    # ref_genome_dir = current_ref_genome["dir"]
    if args["projectID"]:
      project = all_projects[args["projectID"]]

    gene_ids = []
    for datagroup in project["data"]:
      data_type = datagroup['data_type']
      handler = Helpers.get_data_handler(data_type, data_root)
      genes = handler.get_genes(project, datagroup)
      gene_ids.extend(genes)

    # gene_ids = Helpers.load_genome_mapping(project).keys()

    return gene_ids


@api.route("/projects")
class ProjectInfos(Resource):
  def get(self):
    # allProjects = dict(map(Helpers.project_info, Helpers.get_subdirs()))
    # anyway run by each API call
    projectinfo = {project: info for project, info in all_projects.iteritems()}
    return projectinfo


if __name__ == '__main__':
  app.run()


def create(*args, **kwargs):
  return app
