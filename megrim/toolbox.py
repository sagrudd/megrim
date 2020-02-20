#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:22:55 2020

@author: srudd
"""
import inspect
import os
import pkgutil
from megrim.environment import MegrimPlugin
import argparse
from importlib import reload
import logging
from megrim.environment import Flounder
import tempfile
        

class MegrimToolBox:
    
    def __init__(self, plugin_package, target="plugins"):
        self.plugin_package = plugin_package
        self.target = target
        self.plugins = []
        self.reload_plugins()

        
    def reload_plugins(self):
        logging.debug(f"looking for megrim plugins in package {self.plugin_package}")
        if self.target is None:
            self.walk_package(self.plugin_package)
        else:
            self.walk_package(f"{self.plugin_package}.{self.target}")
        
    def walk_package(self, package):
        # print(f'Walk-package ... {package}')
        imported_package = __import__(package, fromlist=[""])
        for _, pluginname, ispkg in pkgutil.iter_modules(
                imported_package.__path__, imported_package.__name__ + "."):
            if not ispkg:
                # print(f"checking {pluginname}")
                plugin_module = __import__(pluginname, fromlist=[""])
                # print(plugin_module)
                clsmembers = inspect.getmembers(plugin_module, inspect.isclass)
                for (_, c) in clsmembers:
                    # print(f"{_}\t\t{c}")
                    if issubclass(c, MegrimPlugin) and (c is not MegrimPlugin):
                        logging.debug(f'\timporting plugin class: {c.__name__}')
                        self.plugins.append(c())

    def list(self):
        code_words = []
        for plugin in self.plugins:
            code_words.append(plugin.tool)
        return ", ".join(code_words)

    def execute(self, args):
        try:
            fubar = True
            for plugin in self.plugins:
                if plugin.tool == args.method:
                    fubar = False
                    plugin.execute(args)
            if fubar:
                raise ValueError(f"the requested method [{args.method}] is not known")
        except ValueError as e:
            print("Exception!", e)

    def arg_params(self, subparsers, parent_parser):
        logging.info("merging annotations ...")
        for plugin in self.plugins:
            plugin.arg_params(subparsers, parent_parser)



def main():
    reload(logging)
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:%(message)s',
        level=logging.INFO, datefmt='%I:%M:%S')

    megrim_plugins = MegrimToolBox("megrim")
    parser = argparse.ArgumentParser()
    # parser.add_argument("method", help=f"Define the megrim method to run]")

    subparsers = parser.add_subparsers(title='subcommand help')
    subparsers.required = True
    subparsers.dest = 'method'

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--cache', metavar="/tmp", action='store', help='Path to location for storing cached and temporary files.', dest="cache", default=tempfile.gettempdir())

    megrim_plugins.arg_params(subparsers, parent_parser)

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")

    args = parser.parse_args()

    # setup a Flounder for this workflow ...
    flounder = Flounder()
    flounder.cache_path = args.cache
    print(flounder.cache_path)

    megrim_plugins.execute(args)
