#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:11:57 2019

@author: srudd
"""


from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from time import time

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
        
        i = 1
        


ssh = SequenceSummaryHandler("/Users/srudd/Desktop/Dolores_PromethION_1M.txt")




 
    
def executive_summary():
    with Timer("Time to extract passed state information {}"):
        seq_sum = Client.persist(seq_sum)
    readCount <- formatC(nrow(sequencedata), big.mark=",")
    totalBases = sum(sequencedata$sequence_length_template,na.rm=T)/10^9
    passedBases = sum(passedSeqs$sequence_length_template,na.rm=T)/10^9
    gigabases <- round(totalBases,2)