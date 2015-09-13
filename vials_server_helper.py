import json
import os

__author__ = 'Hendrik Strobelt'

vials_config_file_name = 'vials_project.json'


def update_project_file(project):
    with open(os.path.join(project['dir'], vials_config_file_name), 'wb') as p_file:
        json.dump(project['info'], p_file, indent=4)
