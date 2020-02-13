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
import functools
import pysam
import sys


class BedHandler:

    def __init__(self, bedfile):
        self.bedfile = bedfile
        self.targets = None
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
        self.samfile = pysam.AlignmentFile(bam, "rb")

    @functools.lru_cache()
    def get_bam_ranges(self, filter_flag=3844):
        logging.info("Extracting BAM ranges")
        # note that read_bam by default has a specific set of SAM flags
        # this code needs to be updated to allow for selection of +/-, primary
        # secondary etc ... - this method is also independent of pysam//samtools
        return pr.read_bam(self.bam, filter_flag=3844)

    def get_bam_coverage(self):
        logging.debug("Extracting BAM coverage")
        return self.get_bam_ranges().to_rle(strand=False)

    def get_sam_reads(self, chromosome, start, end):
        # logging.info("getSamReads ({}) {}:{}".format(chromosome, start, end))
        result = None
        try:
            result = self.samfile.fetch(chromosome, int(start), int(end))
        except Exception:
            result = None
        return result

    def get_sam_core(self, chromosome, start, end):
        annot = self.get_sam_reads(chromosome, start, end)
        return annot

    def get_sam_annotation(self, chromosome, start, end):
        annot = self.get_sam_reads(chromosome, start, end)
        if annot is None:
            return None
        start = []
        reference_length = []
        mapping_quality = []
        strand = []
        mapped_read_q = []
        cigar_m = []
        cigar_i = []
        cigar_d = []
        nm = []
        read_counter = 0
        for read in annot:
            if not any([read.is_secondary, read.is_supplementary,
                        read.is_qcfail, read.is_duplicate]):
                read_counter += 1
                start.append(read.reference_start)
                reference_length.append(read.reference_length)
                mapping_quality.append(read.mapping_quality)
                if read.is_reverse:
                    strand.append("-")
                else:
                    strand.append("+")
                mapped_read_q.append(-10 * log10((10 ** (pd.Series(
                    read.query_alignment_qualities) / -10)).mean()))
                cigar_stats = read.get_cigar_stats()
                cigar_m.append(cigar_stats[0][0])
                cigar_i.append(cigar_stats[0][1])
                cigar_d.append(cigar_stats[0][2])
                nm.append(cigar_stats[0][10])
            annot = {'reference_start': start,
                     'reference_length': reference_length,
                     'mapping_quality': mapping_quality,
                     'strand': strand,
                     'mapped_read_q': mapped_read_q,
                     'cigar_m': cigar_m,
                     'cigar_i': cigar_i,
                     'cigar_d': cigar_d,
                     'nm': nm}
        if read_counter > 0: 
            return pd.DataFrame.from_dict(annot)
        return None
