import datetime
import json
import os
import pytz


class reference_data_registry:
    def __init__(self, source, filename, key):
        self.path = os.path.join(source, filename)
        self.keys = key
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
            "name": "nanopoRe data registry",
            self.keys: []
        }
        self.write()

    def add_item(self, name, url, path):
        newGenome = {
            "timestamp": "" + datetime.datetime.now(tz=pytz.utc).isoformat(),
            "name": name,
            "url": url,
            "path": path
        }
        self.registry.get(self.keys).append(newGenome)
        self.write()

    def list_references(self):
        items = []
        for item in self.registry.get(self.keys):
            items.append(os.path.basename(item.get("path")) + "::" + item.get("name"))
        return items

    def extract_existing(self, option):
        filename, name = option.split("::")
        for item in self.registry.get(self.keys):
            if item.get("name") == name:
                return item.get("path")
        return None

    def getIdFromPath(self, path):
        print(path)
        for item in self.registry.get(self.keys):
            if item.get("path") == path:
                return os.path.basename(item.get("path")) + "::" + item.get("name")
