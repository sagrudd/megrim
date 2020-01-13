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

from pkg_resources import resource_string

class Flounder:
    """
    FloundeR was the name of the R package that managed tutorial content,
    workflows and logic - deprecated with move to Python - retained here
    as a memorial to a great package name.
    
    This now unintuitively named package is responsible for the orchestration
    of workspace paths, directory structures and is aiming towards
    persistence of analysis states where possible
    """

    def __init__(self):
        self.location = ''

    def set_path(self, path):
        self.location = path

    def get_path(self):
        return self.location
    


def get_branding_logo():
    return resource_string(__name__, 'data/ONT_logo.png')