#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 08:56:25 2019

@author: srudd

This is a functional copy of the environment controls provided in the 
nanopoRe package - for occasions such as in the development of EPI2ME-labs,
a user may be working in a non-persistent workspace with primary data,
scripts and results in a peripheral location
"""

import logging
import re
from pkg_resources import resource_string, get_distribution, resource_filename
from IPython.display import Image, display, Markdown 
from bokeh.plotting import show
from bokeh.io import export_png
from bokeh.io.export import get_screenshot_as_png
import pathlib
import sys
import os
import tempfile


class Flounder:
    """
    FloundeR was the name of the R package that managed tutorial content,
    workflows and logic - deprecated with move to Python - retained here
    as a memorial to a great package name.
    
    This now unintuitively named package is responsible for the orchestration
    of workspace paths, directory structures and is aiming towards
    persistence of analysis states where possible
    """

    def __init__(self, plot_width=640, plot_height=480, plot_type="native", 
                 plot_tools="save,reset", plot_dpi=96):
        self.location = None
        self.results_dir = None
        self.plot_width = plot_width
        self.plot_height = plot_height
        self.plot_type = plot_type
        self.plot_tools = plot_tools
        self.plot_dpi = plot_dpi

    def set_path(self, path):
        self.location = path

    def get_path(self):
        return self.location
    
    def set_results_dir(self, results_dir):
        self.results_dir = results_dir
        
    def get_results_dir(self):
        return self.results_dir
    
    def set_plot_width(self, plot_width):
        self.plot_width = plot_width
        
    def set_plot_height(self, plot_height):
        self.plot_height = plot_height
        
    def set_plot_type(self, plot_type):
        self.plot_type = plot_type
        
    def set_plot_tools(self, plot_tools):
        self.plot_tools = plot_tools
        
    def set_plot_dpi(self, plot_dpi):
        self.plot_dpi = plot_dpi
        
    def get_plot_height(self):
        return self.plot_height
    
    def get_plot_width(self):
        return self.plot_width
    
    def get_plot_type(self):
        return self.plot_type
    
    def get_plot_tools(self):
        return self.plot_tools
    
    def get_plot_dpi(self):
        return self.plot_dpi
    
    def get_destination(self):
        destination_dir = self.get_path()
        if destination_dir is None:
            logging.error("Flounder requires a path for writing files to")
            sys.exit(0)
        if self.get_results_dir() is not None:
            destination_dir = os.path.join(destination_dir, self.get_results_dir())
        logging.info("destination --> %s" % destination_dir)
        if not os.path.exists(destination_dir):
            logging.warning("creating DIR==%s" % destination_dir)
            os.makedirs(destination_dir)
        return destination_dir
        
    def handle_output(self, p, plot_type, prefix="bokeh_"):
        if plot_type == 'native':
            return p
        elif plot_type == 'jupyter':
            show(p)
            return None
        elif plot_type == 'png':
            logging.debug("exporting figure to png")
            fd, temp_path = tempfile.mkstemp(dir=self.get_destination(), suffix=".png", prefix=prefix)
            os.close(fd)
            logging.info("tmppath@%s" % temp_path)
            export_png(p, filename=temp_path)
            display(Image(temp_path))
            return temp_path
        elif plot_type == 'screenshot':
            logging.debug("exporting figure to screenshot")
            image = get_screenshot_as_png(p)
            return image
        else:
            logging.error("unknown plottype")
            sys.exit(0)
            return "unknown plottype"
    
    def handle_kwargs(self, fields, **kwargs):
        stuff = []
        for f in fields:
            if f in kwargs.keys():
                stuff.append(kwargs.get(f))
            else:
                if f == "plot_width":
                    stuff.append(self.get_plot_width())
                elif f == "plot_height":
                    stuff.append(self.get_plot_height())
                elif f == "plot_type":
                    stuff.append(self.get_plot_type())
                elif f == "plot_tools":
                    stuff.append(self.get_plot_tools())
                elif f == "plot_dpi":
                    stuff.append(self.get_plot_dpi())
                else:
                    logging.warning("[%s] is an unknown and undefined variable" % (f))
                    stuff.append(None)
        tres = tuple(stuff)
        logging.debug(tres)
        return tres
    

def get_megrim_version():
    try:
        distribution = get_distribution("megrim")
        print("%s(%s)" % (distribution.key, distribution.version))
    except:
        logging.warning("megrim package may not be installed")
        txt_chunk = open((pathlib.Path(resource_filename('megrim', 'data')).parent.parent / "setup.py").as_posix()).read()
        txt_lines = txt_chunk.splitlines()
        print("%s(%s)" % (
            re.findall('(?<=name=\')[^\']+', list(filter(lambda x: 'name=' in x, txt_lines))[0])[0], 
            re.findall('(?<=version=\')[^\']+', list(filter(lambda x: 'version=' in x, txt_lines))[0])[0]))


def tutorial_branding(tutorial=None, legend=None, 
                      telemetry=None):
    
    def get_framework():
        try:
            ipy_str = str(type(get_ipython()))
            if 'zmqshell' in ipy_str:
                return "jupyter"
            if 'terminal' in ipy_str:
                return 'ipython'
        except:
            return "terminal"
        
    framework = get_framework()
    if framework == "jupyter":
        display(Image(resource_string(__name__, 'data/ONT_logo.png')))
        if legend != None:
            display(Markdown("# "+legend))


            

        