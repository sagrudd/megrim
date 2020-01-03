#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:52:28 2020

@author: srudd
"""


from pysam import FastaFile
import pyranges as pr
import numpy as np

class ReferenceGenome:
    
    def __init__(self, referenceFasta):
        self.fasta = FastaFile(filename=referenceFasta)
        self.allowed = []
        self.skipped = []
        
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
        
