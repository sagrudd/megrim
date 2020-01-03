#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:44:30 2019

@author: srudd
"""

import logging
from math import log10
import pandas as pd
import pyranges as pr


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




class BedHandler:
    
    def __init__(self, bedfile):
        self.bedfile = bedfile
        self.proximity = 0
        self.ref = None
        self.load_bed_file()
        
    def load_bed_file(self):
        self.targets = pd.read_csv(self.bedfile, header=None, 
                                   names="Chromosome Start End Name".split(), 
                                   sep='\t')
        
    def get_bed_targets(self):
        print(self.targets)
        
    def set_target_proximity(self, proximity):
        self.proximity = proximity
        
    def get_target_proximity(self):
        return self.proximity
    
    def set_reference(self, ref):
        self.ref = ref
        
    def get_reference(self):
        return self.ref
    
    def get_target_ranges(self):
        return pr.PyRanges(self.targets)
    
    def get_target_proximal_ranges(self):
        downstream = self.get_target_ranges()
        upstream = self.get_target_ranges()
       
        upstream.End = upstream.Start - 1
        upstream.Start += -(self.get_target_proximity() + 1)
        # there may be edge exceptions where the Start coordinate < 0?
        
        downstream.Start = downstream.End + 1
        downstream.End += self.get_target_proximity() + 1
        # there may be edge exceptions where End coordinate drops off chromo.
        
        merged = pr.concat([upstream, downstream])
        return pr.gf.genome_bounds(merged, self.ref.get_reference_ranges(), clip=True)
    
    def get_untargeted_ranges(self):
        whole_genome_ranges = self.ref.get_reference_ranges()
        
        target_proximal_ranges = self.get_target_ranges()
        target_proximal_ranges.Start += -(self.get_target_proximity() + 1)
        target_proximal_ranges.End += self.get_target_proximity() + 1
        # there is a possibility for edge cases where boundaries of
        # target_proximal_ranges extend beyond the limits of the chromosome
        untargeted = whole_genome_ranges.subtract(target_proximal_ranges)
        return pr.gf.genome_bounds(untargeted, self.ref.get_reference_ranges(), clip=True)

    

    
    
class BamHandler:
    
    def __init__(self, bam):
        self.bam = bam

    def get_bam_ranges(self):
        return pr.read_bam(self.bam)
    
    def get_bam_coverage(self):
        return self.get_bam_ranges().to_rle(strand=False)

    
    
        