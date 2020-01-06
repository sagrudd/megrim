#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:52:28 2020

@author: srudd
"""


from pysam import FastaFile
import pyranges as pr
import numpy as np
import pandas as pd
from environment import Flounder
from tqdm import tqdm
#from tqdm.auto import tqdm

class ReferenceGenome:
    
    def __init__(self, referenceFasta):
        self.fasta = FastaFile(filename=referenceFasta)
        self.allowed = []
        self.skipped = []
        try:
            globals()['flounder']
        except KeyError:
            print("flounder environment not yet available ..")
            global flounder
            flounder = Flounder()
        
    def info(self):
        print(self.fasta.is_open())
        
        print(self.fasta.filename)
        
        references = self.fasta.references
        for reference in references:
            size = self.fasta.get_reference_length(reference)
            print("ref %s == %s" % (reference, size))
        
    def get_chromosomes(self):
        chromosomes = self.fasta.references
        approved_chromosomes = []
        for chromosome in chromosomes:
            if not chromosome in self.skipped:
                if len(self.allowed) == 0:
                    approved_chromosomes.append(chromosome)
                elif chromosome in self.allowed:
                    approved_chromosomes.append(chromosome)
        return approved_chromosomes
        
    def get_chromosome_lengths(self, chromosomes):
        lengths = []
        for chromosome in chromosomes:
            lengths.append(self.fasta.get_reference_length(chromosome))
        return lengths
    
    def get_reference_ranges(self):
        chromosomes = self.get_chromosomes()
        starts = np.repeat(0, len(chromosomes))
        ends = self.get_chromosome_lengths(chromosomes)
        gr = pr.PyRanges(chromosomes=chromosomes, starts=starts, ends=ends)
        return gr
    
    def skip_chromosome(self, chr):
        if not chr in self.skipped:
            self.skipped.append(chr)
            
    def skip_chromosomes(self, chrs):
        for chr in chrs:
            self.skip_chromosome(chr)
            
    def get_skipped(self):
        print(" ".join(self.skipped))
        
    def get_tiled_genome(self, tile_size=100):
        return pr.gf.tile_genome(self.get_reference_ranges(), tile_size, tile_last=False)
    
    def get_tiled_coverage(self, bam_ranges, tile_size=100):
        return self.get_tiled_genome(tile_size=tile_size).coverage(bam_ranges)
    
    def get_tiled_mean_coverage(self, bam_ranges, bam_rle, tile_size=100):
        tiled_coverage = self.get_tiled_coverage(bam_ranges, tile_size)
        print(tiled_coverage)
        meanList = []
        for row in tqdm(tiled_coverage.df.itertuples(), total=tiled_coverage.df.shape[0]):
            n = row.NumberOverlaps
            if n > 0:
                #print(row)
                chromosome = row.Chromosome
                start = row.Start
                end = row.End
                
                rle = bam_rle[chromosome][start:end]
                mean = np.repeat(rle.values,rle.runs).mean()

                meanList.append(mean)
            else:
                meanList.append(0)
        tiled_coverage.MeanCov = meanList
        return tiled_coverage
            
    def get_tiled_mean_coverage2(self, bam_ranges, bam_rle, tile_size=100):
        # this doesn't have any performance benefit over a df.itertuples
        tiled_coverage = self.get_tiled_coverage(bam_ranges, tile_size)            
        def chunk(row):
            if row.NumberOverlaps == 0:
                return 0                
            rle = bam_rle[row.Chromosome][row.Start : row.End]
            return np.repeat(rle.values,rle.runs).mean()
        tqdm.pandas()
        tiled_coverage.MeanCoverage = tiled_coverage.df.progress_apply(chunk, axis=1)
        return tiled_coverage
    
    def get_tiled_mean_coverage3(self, bam_ranges, bam_rle, tile_size=100):
        tiled_coverage = self.get_tiled_coverage(bam_ranges, tile_size)
        # precompute the per-base depths of coverage ...
        base_data = {}
        chromosomes = self.get_chromosomes()
        for chromosome in tqdm(chromosomes):
            base_data[chromosome] = np.repeat(
                bam_rle[chromosome].values, bam_rle[chromosome].runs)
        def chunk(row):
            if row.NumberOverlaps == 0:
                return 0                
            return base_data[row.Chromosome][row.Start : row.End].mean()
        tqdm.pandas()
        tiled_coverage.MeanCoverage = tiled_coverage.df.progress_apply(chunk, axis=1)
        return tiled_coverage
    
    def get_tiled_mean_coverage4(self, bam_ranges, bam_rle, tile_size=100):
        tiled_coverage = self.get_tiled_coverage(bam_ranges, tile_size)
        tiled_coverage.MeanCoverage = 0
        
        df_data = tiled_coverage.df
        chromosomes = self.get_chromosomes()
        for chromosome in tqdm(chromosomes):
            base_data = np.repeat(bam_rle[chromosome].values, 
                                  bam_rle[chromosome].runs)
            def chunk(row):
                return base_data[row.Start : row.End].mean()
            
            scores = df_data.loc[((df_data['Chromosome']==chromosome) & (df_data['NumberOverlaps']>0)),:].apply(chunk, axis=1)
            df_data.loc[((df_data['Chromosome']==chromosome) & (df_data['NumberOverlaps']>0)),('MeanCoverage')]=scores
        return pr.PyRanges(df_data)
            
