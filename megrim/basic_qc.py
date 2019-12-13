#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:11:57 2019

@author: srudd
"""

import logging
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from time import time
from infographic_plots import InfographicPlot, InfographicNode

class Timer():
    """
    https://stackoverflow.com/questions/3620943/measuring-elapsed-time-with-the-time-module/46544199
    """
    
    def __init__(self, message):
        self.message = message
        
    def __enter__(self):
        self.start = time()
        return None  # could return anything, to be used like this: with Timer("Message") as value:

    def __exit__(self, type, value, traceback):
        e = int(time() - self.start)
        m = '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60)
        print(self.message.format(m))
        


class SequenceSummaryHandler:
    
    def __init__(self, target_file):
        self.target_file = target_file
        self._import_data()
        
    def _import_data(self):
        """
        method parses the provided Sequencing_summary.txt file and loads 
        the reduced contents into memory - a DASK dataframe is used to
        allow for the scaling to large PromethION runs (or equivalent)

        Returns
        -------
        None.

        """
        self.seq_sum = dd.read_csv(
            self.target_file, 
            delimiter='\t',
            blocksize=64000000,
        )
        # start excluding dask columns that are not of core interest
        keep = ['channel', 'start_time', 'duration', 'num_events', 
                'sequence_length_template', 'mean_qscore_template', 
                'passes_filtering',
        ]
        for col in self.seq_sum.columns:
            if not col in keep:
                print("dropping %s" % col)
                self.seq_sum = self.seq_sum.drop(col, axis=1)
        pbar = ProgressBar()
        pbar.register()
        with Timer("Elapsed time to extract sequence data {}"):
                self.seq_sum = self.seq_sum.persist()
        pbar.unregister()

    def executive_summary(self):
        read_count = len(self.seq_sum)
        total_bases = self.seq_sum['sequence_length_template'].sum().compute()
        # passed_bases = self.seq_sum[self.seq_sum['passes_filtering']]['sequence_length_template'].sum().compute()
        
        
        flowcell_node = InfographicNode(legend="flowcell", 
                                value="FAK85195", 
                                graphic='fingerprint')

        readcount_node = InfographicNode(legend="Reads produced", 
                                value="{:,}".format(read_count),
                                graphic='filter')

        gb_seq_val = InfographicNode(legend="Gigabases called", 
                                value="{:.2f}".format(total_bases/1e9),
                                graphic='flag-checkered')
        infographic_data = [flowcell_node, readcount_node, gb_seq_val]
        ip = InfographicPlot(infographic_data, rows=1, columns=3)  
        ip.plot_infographic() 


nul = """
ssh = SequenceSummaryHandler("/Users/srudd/Desktop/Dolores_PromethION_1M.txt")
ssh.executive_summary()
"""


class SequencingSummaryGetChannelMap:
    """
    This class is responsible for the handling of flowcell channel maps -
    the prototype of this class was written in R (see @sagrudd/nanopoRe);
    code has been transposed and simplified
    """
    
    def __init__(self, seq_sum):
        self.seq_sum = seq_sum
        
    def get_platform(self):
        """
        method scores the defined channels in the provided sequencing_summary
        to look for the largest defined channel - based on the number observed
        reports whether it is most likely to be Flongle / MinION or PromethION

        Returns
        -------
        String representation of flowcell type (MinION/Flongle/PromethION)

        """

        platform = "MinION"
        max_channel = self.seq_sum['channel'].max()

        if max_channel < 130:
            # this is likely to be a Flongle ...
            platform = "Flongle"

        if max_channel > 1000:
            # this is likely to be a PromethION
            platform = "PromethION"
        logging.debug("flowcell_type identified as %s" % platform)  
        return platform

        
    
    