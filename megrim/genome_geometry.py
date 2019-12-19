#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:44:30 2019

@author: srudd
"""

import logging
from math import log10
import pandas as pd


class GenomeGeometry:
    
    def __init__(self, list_data):
        self.list_data = list_data.sort_values(ascending=False).reset_index(drop=True)
        self.accumulated = self.list_data.cumsum()
        
    def get_longest_read(self):
        return self.list_data.max()
        
    def get_mean_length(self):
        return self.list_data.mean()
    
    def get_lengths(self):
        return self.list_data
    
    def get_n_value(self, n=50):
        n_sum = int(self.list_data.sum())
        n_targ = n_sum * (n/100)
        logging.debug("n: %d" % n_sum)
        logging.debug("n_val: %d" % n_targ)
        index_pos = self.accumulated.loc[(self.accumulated >= n_targ)].index[0]
        logging.debug("index_pos: %d" % index_pos)
        logging.debug(self.list_data)
        logging.debug("N%d: %d" % (n, self.list_data[index_pos]))
        return self.list_data[index_pos]
    
    def calculate_mean_quality(self, series):
        series = series.dropna()
        return -10 * log10((10 ** (series / -10)).mean())

        