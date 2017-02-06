from __future__ import print_function
import os
import json
from phovea_server.apiutil import create_api, Resource
import phovea_server.config
from .data_handlers.miso_handler import MisoHandler
from .ref_genome_handler.gff3_handler import GFFHandler

__author__ = 'Hendrik Strobelt'

vials_config_file_name = 'vials_project.json'
samples_meta_file_name = 'samples.json'
vials_project_dir_ending = ".vials_project"

# create the api application
app, api = create_api(__name__, version='1.5', title='Caleydo Vials API', description='splicing data handler')
parser = api.parser()
parser.add_argument('geneID', type=str, help='gene ID')
parser.add_argument('projectID', type=str, help='projectID')
parser.add_argument('noReads', type=bool, help='do not show reads')
parser.add_argument('selectFilter', type=str, help='selectFilter')
parser.add_argument('exactMatch', type=bool, help='does match exactly')

projects_dir = phovea_server.config.view('vials_server').projects_dir
if not os.path.exists(projects_dir):
  exit(1)

# global variables for simple caching
project_info = {}
ref_genome_info = {}

data_handler_mapping = {'miso': MisoHandler}

ref_genome_handler_mapping = {'gff3': GFFHandler}


class Helpers:
  @staticmethod
  def get_data_handler(project):
    p_type = project['info']['project_type']
    if p_type in data_handler_mapping:
      return data_handler_mapping[p_type](project)
    return None

  @staticmethod
  def get_ref_genome_handler(project):
    if ref_genome_info == {}:
      # collect all reference genomes
      ref_genome_dir = phovea_server.config.view('vials_server').ref_genomes_dir
      for root, dirs, files in os.walk(ref_genome_dir):
        if 'meta.json' in files:
          with open(os.path.join(root, 'meta.json')) as meta_file:
            m_info = json.load(meta_file)
            ref_genome_info[m_info['name']] = {'info': m_info,
                                               'db_file': os.path.join(root, m_info['sqlite_file']),
                                               'origin': m_info['origin_type']
                                               }

    if project['info']['ref_genome'] in ref_genome_info:
      rg = ref_genome_info[project['info']['ref_genome']]
      if rg['origin'] in ref_genome_handler_mapping:
        return ref_genome_handler_mapping[rg['origin']](rg)

        # return ref_genome_info[project['info']['ref_genome']]

    return None

  @staticmethod
  def get_all_projects(refresh=False):
    global project_info
    if refresh or project_info == {}:
      project_info = {}
      for root, dirs, files in os.walk(projects_dir):
        print(root)
        if root.endswith(vials_project_dir_ending) and samples_meta_file_name in files:
          # info = {}
          with open(os.path.join(root, vials_config_file_name)) as p_file:
            info = json.load(p_file)

          project_id = os.path.basename(root).replace(vials_project_dir_ending, '')
          try:
            name = info.name
          except AttributeError:
            name = project_id

          project_info[project_id] = {'name': name, 'project_id': project_id, 'info': info, 'dir': root}

    return project_info


@api.route("/projects")
class ProjectInfos(Resource):
  def get(self):
    return Helpers.get_all_projects(True)


@api.route("/genes")
class GeneOverview(Resource):
  def get(self):
    args = parser.parse_args()
    project_id = args.projectID
    all_projects = Helpers.get_all_projects()

    if project_id and project_id in all_projects:
      project = all_projects[project_id]
      handler = Helpers.get_data_handler(project)
      return handler.get_genes_in_project()

    return None


@api.route("/geneselect")
class GeneSelectView(Resource):
  def get(self):
    args = parser.parse_args()
    project_id = args.projectID
    all_projects = Helpers.get_all_projects()
    exact_match = False
    if args.exactMatch:
      exact_match = True

    if project_id and project_id in all_projects:
      project = all_projects[project_id]
      handler = Helpers.get_data_handler(project)
      return handler.get_genes_in_project_filtered(args.selectFilter, exact_match)

    return None


@api.route("/geneinfo")
class GeneInfo(Resource):
  def get(self):
    args = parser.parse_args()
    project_id = args.projectID
    all_projects = Helpers.get_all_projects()

    if project_id and project_id in all_projects:
      project = all_projects[project_id]
      data_handler = Helpers.get_data_handler(project)

      ref_genome_handler = Helpers.get_ref_genome_handler(project)

      gene_id = args.geneID
      gene_meta = ref_genome_handler.get_gene_info(gene_id)

      samples, iso_measures, jxns_measures, wiggles, wiggle_sample_size = data_handler.get_samples_and_measures(gene_id,
                                                                                                                gene_meta)

      return {'gene': gene_meta,
              'samples': samples,
              'measures': {'data_type': 'miso',
                           'isoform_unit': 'scaled estimate (TPM)',
                           'wiggle_sample_size': wiggle_sample_size,
                           'isoforms': iso_measures,
                           'jxns': jxns_measures,
                           'wiggles': wiggles
                           }
              }

    return None


if __name__ == '__main__':
  app.run()


def create(*args, **kwargs):
  return app
