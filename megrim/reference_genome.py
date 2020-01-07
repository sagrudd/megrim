#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:52:28 2020

@author: srudd
"""


from pysam import FastaFile
import pyranges as pr
import numpy as np
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
    
    def get_tiled_mean_coverage(self, bam, ranges=None, tile_size=100):
        if ranges is None:
            ranges = self.get_tiled_coverage(bam.get_bam_ranges(), tile_size)
        else:
            if 'NumberOverlaps' not in ranges.columns:
                ranges = ranges.coverage(bam.get_bam_ranges())
        # create a clear dataset for coverage ...
        ranges.MeanCoverage = 0
        # create data objects required
        bam_rle = bam.get_bam_coverage()
        df_data = ranges.df
        # pull chromosome id from the provided ranges
        chromosomes = df_data['Chromosome'].unique()
        for chromosome in tqdm(chromosomes):
            base_data = np.repeat(bam_rle[chromosome].values, 
                                  bam_rle[chromosome].runs)
            def chunk(row):
                return base_data[row.Start : row.End].mean()
            scores = df_data.loc[((df_data['Chromosome']==chromosome) & (df_data['NumberOverlaps']>0)),:].apply(chunk, axis=1)
            df_data.loc[((df_data['Chromosome']==chromosome) & (df_data['NumberOverlaps']>0)),('MeanCoverage')]=scores
        return pr.PyRanges(df_data)
    
    
    
    
    
    def augment_annotation(self, bam, ranges):
        # rstart - the number of reads that start within the given interval
        
        # basesstart - the number of bases contained within rstart
        
        # meanreadlen - mean read length for any reads within this interval
        
        # startreadlen - mean read length for reads that start within interval
        
        # strandp
        
        # strandn
        
        # mapq - mapq for reads starting in segment
        
        # map0 - mapq for reads overlapping the segment
        
        # readq - per read q score for reads starting in segment 
        
        # read0 - per read q score for reads overlapping segment 
        
        # nm - this is the #NM mismatch count; reads starting in segment
        
        # cigar_del
        
        # cigar_ins
        
        # cigar_mapped
        
        ##### and some local sequence context annotations
        
        # gccount
        
        # ncount
        
        return ranges
    
        
            
