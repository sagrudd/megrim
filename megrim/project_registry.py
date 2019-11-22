import datetime
import json
import os
import pytz
from .accessory_functions import *


class project_registry:
    def __init__(self, source):
        target_dir(source)
        self.path = os.path.join(source, "projects.json")
        self.read()

    def read(self):
        if not os.path.exists(self.path):
            print("Registry has not been instantiated ...")
            self.create()
        else:
            with open(self.path) as f:
                self.registry = json.load(f)

    def write(self):
        with open(self.path, 'w') as json_file:
            json.dump(self.registry, json_file)

    def create(self):
        self.registry = {
            "name": "nanopoRe project registry",
            "projects": []
        }
        self.write()

    def add_project(self, name, path, project_type):
        newProject = {
            "timestamp": "" + datetime.datetime.now(tz=pytz.utc).isoformat(),
            "name": name,
            "path": path,
            "project_type": project_type
        }
        self.registry.get('projects').append(newProject)
        self.write()

    def list_projects(self, project_type):
        items = []
        for item in self.registry.get('projects'):
            if item.get("project_type") == project_type:
                items.append(item.get("name"))
        return items

    def get_project_path(self, project_name):
        for item in self.registry.get('projects'):
            if item.get("name") == project_name:
                return item.get("path")
        return None
