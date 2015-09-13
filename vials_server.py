import os
import json

from flask.ext.restplus import Resource

from caleydo_server.apiutil import create_api
import caleydo_server.config
from handlers.miso_handler import MisoHandler

__author__ = 'Hendrik Strobelt'


vials_config_file_name = 'vials_project.json'
samples_meta_file_name = 'samples.json'
vials_project_dir_ending = ".vials_project"


# create the api application
app, api = create_api(__name__, version='1.5', title='Caleydo Vials API', description='splicing data handler')
parser = api.parser()
parser.add_argument('geneID', type=str, help='gene ID')
parser.add_argument('noReads', type=bool, help='do not show reads')
parser.add_argument('projectID', type=str, help='projectID')
parser.add_argument('selectFilter', type=str, help='selectFilter')

projects_dir = caleydo_server.config.view('vials_server').projects_dir
if not os.path.exists(projects_dir):
    exit(1)

# global variable for simple caching
project_info = {}


class Helpers:
    @staticmethod
    def get_data_handler(project):
        if project['info']['project_type']== "miso":
            return MisoHandler(project)
        return None

    @staticmethod
    def get_all_projects(refresh=False):
        global project_info
        if refresh or project_info == {}:
            project_info = {}
            for root, dirs, files in os.walk(projects_dir):
                print root
                if root.endswith(vials_project_dir_ending) and samples_meta_file_name in files:
                    # info = {}
                    with open(os.path.join(root, vials_config_file_name)) as p_file:
                        info = json.load(p_file)

                    project_id = os.path.basename(root).replace(vials_project_dir_ending, '')
                    try:
                        name = info.name
                    except AttributeError:
                        name = project_id

                    project_info[project_id] = {
                        'name': name,
                        'project_id': project_id,
                        'info': info,
                        'dir': root
                    }

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

        if project_id and project_id in all_projects:
            project = all_projects[project_id]
            handler = Helpers.get_data_handler(project)
            return handler.get_genes_in_project_filtered(args.selectFilter)

        return None


if __name__ == '__main__':
    app.run()


def create(*args, **kwargs):
    return app