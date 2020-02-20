#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:22:55 2020

@author: srudd
"""
import inspect
import os
import pkgutil

class MegrimPlugin:
    
    def __init__(self):
        self.tool = "toolname"
        
        
        

class MegrimToolBox(object):
    
    def __init__(self, plugin_package):
        self.plugin_package = plugin_package
        self.reload_plugins()
        
    def reload_plugins(self):
        self.plugins = []
        self.seen_paths = []
        print(f"looking for plugins in package {self.plugin_package}")
        self.walk_package(self.plugin_package)
        
    def walk_package(self, package):
        imported_package = __import__(package, fromlist=[])
        
        for _, pluginname, ispkg in pkgutil.iter_modules(
                imported_package.__path__, imported_package.__name__ + "."):
            if not ispkg:
                plugin_module = __import__(pluginname, fromlist=[])
                clsmembers = inspect.getmembers(plugin_module, inspect.isclass)
                for (_, c) in clsmembers:
                    if issubclass(c, MegrimPlugin) and (c is not MegrimPlugin):
                        print(f'\tFound plugin class: {c.__module__}.{c.__name__}')
                        self.plugins.append(c())
                        
        all_current_paths = []
        if isinstance(imported_package.__path__, str):
            all_current_paths.append(imported_package.__path__)
        else:
            all_current_paths.extend([x for x in imported_package.__path__])
            
        for pkg_path in all_current_paths:
            if pkg_path not in self.seen_paths:
                self.seen_paths.append(pkg_path)
                
                child_pkgs = [p for p in os.listdir(pkg_path) if os.path.isdir(os.path.join(pkg_path, p))]
        
                for child_pkg in child_pkgs:
                    self.walk_package(package+"."+child_pkg)

def main():
    print("Hello world")