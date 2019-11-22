import datetime
import json
import os
import pytz


class ProjectManager:
    def __init__(self, source, project_type, project_name):
        self.path = os.path.join(source, "nanopoRe_project.json")
        self.type = project_type
        self.name = project_name
        self.read()

    def read(self):
        if not os.path.exists(self.path):
            print("project setup not been instantiated ...")
            self.create()
        else:
            with open(self.path) as f:
                self.project = json.load(f)

    def create(self):
        self.project = {
            "meta-type": "nanopoRe project",
            "timestamp": "" + datetime.datetime.now(tz=pytz.utc).isoformat(),
            "project-type": self.type,
            "project-path": os.path.dirname(self.path),
            "project-name": self.name
        }
        self.write()

    def write(self):
        with open(self.path, 'w') as json_file:
            json.dump(self.project, json_file)

    def has_defined_value(self, key):
        return key in self.project

    def get_defined_value(self, key):
        return self.project[key]

    def set_defined_value(self, key, value):
        self.project[key] = value
        self.write()

    def get_path(self):
        return self.project["project-path"]

    def get_name(self):
        return self.project["project-name"]
