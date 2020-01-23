#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:52:28 2020

@author: srudd
"""

from pysam import FastaFile
import pyranges as pr
import numpy as np
from megrim.environment import Flounder
from tqdm.notebook import tqdm
from math import log10
import pandas as pd
import logging
import warnings
import swifter
import weightedcalcs as wc

# from tqdm.auto import tqdm

flounder = None


class ReferenceGenome(Flounder):

    def __init__(self, reference_fasta):
        Flounder.__init__(self)
        self.fasta = FastaFile(filename=reference_fasta)
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
    
    def get_genome_size(self):
        size = 0
        for chromosome in self.get_chromosomes():
            size += self.fasta.get_reference_length(chromosome)
        return size
        

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
        return pr.gf.tile_genome(
            self.get_reference_ranges(), tile_size, tile_last=False)

    def get_tiled_coverage(self, bam_ranges, tile_size=100):
        return self.get_tiled_genome(tile_size=tile_size).coverage(bam_ranges)

    def get_tiled_mean_coverage(self,bam,ranges=None, tile_size=100):
        if ranges is None:
            ranges = self.get_tiled_genome(tile_size)
        rle = bam.get_bam_coverage()

        df = rle[ranges].df

        df['VR'] = df['Run'] * df['Value']
        df = df.groupby(["Chromosome", "Start"]).agg(
            {"Chromosome": "first", "Start": "first", "End": "first", 
             "Run":np.sum, "VR":np.sum})
        df['MeanCoverage'] = df['VR'] / df['Run']
        return pr.PyRanges(
            df.reset_index(drop=True).drop(["Run", "VR"], axis=1))

    def stranded_dive(self, bam, arange, target, target_proximity=5000, window_size=10):
        targets = arange.df
        targets = targets.loc[targets.Name == target,]
        
        chromos = str(targets.Chromosome.tolist()[0])
        starts = int(targets.Start.tolist()[0]) - target_proximity
        ends = int(targets.End.tolist()[0]) + target_proximity
        
        myrange = pr.PyRanges(pd.concat([pr.PyRanges(chromosomes=[chromos], 
                               starts=[starts], 
                               ends=[ends],
                               strands=["+"]).df,
                             pr.PyRanges(chromosomes=[chromos], 
                               starts=[starts], 
                               ends=[ends],
                               strands=["-"]).df]))
        xrange = pr.gf.tile_genome(
                myrange, window_size, tile_last=False)
    
        bam_data = bam.get_bam_ranges().df
        #bam_data = bam_data.loc[bam_data.Strand=="+",]
        bam_data = pr.PyRanges(bam_data)

        rle = bam_data.to_rle(strand=True)
        df = rle[xrange].df
        
        df['VR'] = df['Run'] * df['Value']
        df = df.groupby(["Chromosome", "Start", "Strand"]).agg(
            {"Chromosome": "first", "Start": "first", "End": "first", 
             "Strand": "first", "Run":np.sum, "VR":np.sum})
        df['MeanCoverage'] = df['VR'] / df['Run']
        # pivot the data ...
        df = df.reset_index(drop=True)
        df = pd.pivot_table(
            df, values = 'MeanCoverage', 
            index=['Chromosome','Start', 'End'], 
            columns = 'Strand').reset_index()
        df['MeanCoverage'] = df["+"] + df["-"]
        return pr.PyRanges(df)
        
def weighted_percentile(universe, q=0.5):
    df = pd.DataFrame({'v': universe.MeanCoverage.tolist(), 'w': (universe.End - universe.Start).tolist()})
    calc = wc.Calculator('w')  # w designates weight
    return calc.quantile(df, 'v', q)


def augment_annotation(bam, ranges):
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    
        def extract_annot(row):
            # bam_data['reference_start'] >= 155179779
            # start_data = bam_data[bam_data['reference_start'] >= 155179779]
            bam_data = bam.get_sam_annotation(row.Chromosome, row.Start, row.End)
            start_data = bam_data[bam_data['reference_start'] >= row.Start]
            # rstart - the number of reads that start within the given interval
            rstart = len(start_data)
            # basesstart - the number of bases contained within rstart
            bases_start = start_data.reference_length.sum()
            # meanreadlen - mean read length for any reads within this interval
            mean_read_len = bam_data.reference_length.mean()
            # startreadlen - mean read length for reads that start within interval
            start_read_len = start_data.reference_length.mean()
            # strandp
            strand_p = (bam_data.strand == '+').sum()
            # strandn
            strand_n = (bam_data.strand == '-').sum()
            # mapq - mapq for reads starting in segment
            mapq = (-10 * log10((10 ** (start_data.mapping_quality / -10)).mean()))
            # map0 - mapq for reads overlapping the segment
            map0 = (-10 * log10((10 ** (bam_data.mapping_quality / -10)).mean()))
            # readq - per read q score for reads starting in segment
            readq = (-10 * log10((10 ** (start_data.mapped_read_q / -10)).mean()))
            # read0 - per read q score for reads overlapping segment
            read0 = (-10 * log10((10 ** (bam_data.mapped_read_q / -10)).mean()))
            # nm - this is the #NM mismatch count; reads starting in segment
            nm = start_data.nm.sum()
            # cigar_del
            cigar_d = start_data.cigar_d.sum()
            # cigar_ins
            cigar_i = start_data.cigar_i.sum()
            # cigar_mapped
            cigar_m = start_data.cigar_m.sum()
            ##### and some local sequence context annotations
    
            # gccount
    
            # ncount
    
            return rstart, bases_start, mean_read_len, start_read_len, \
                   strand_p, strand_n, mapq, map0, readq, read0, nm, cigar_m, \
                   cigar_i, cigar_d
    
        df_data = ranges.df
        
        df_data[['rstart', 'bases_start', 'mean_read_len', 'start_read_len',
                 'strand_p', 'strand_n', 'mapq', 'map0', 'readq', 'read0',
                 'nm', 'cigar_m', 'cigar_i', 'cigar_d']
                ] = df_data.swifter.progress_bar(enable=True, desc=None).apply(
                    extract_annot, axis=1, result_type='expand')
        return pr.PyRanges(df_data)
