#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:44:30 2019

@author: srudd
"""

import logging
from math import log10
from megrim.environment import Flounder
import pandas as pd
import pyranges as pr
import functools
import pysam
from tqdm import tqdm
import numpy as np
import os
import sys
import argparse
from bokeh.plotting import figure
from bokeh.models import Span, NumeralTickFormatter, LinearAxis

flounder = None


def include_flounder(args):
    # setup a Flounder for this workflow ...
    global flounder
    flounder = Flounder()
    if isinstance(args, argparse.Namespace):
        flounder.argparse(args)
    if isinstance(args, dict):
        flounder.dictparse(args)


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


class BamHandler(Flounder):

    def __init__(self, bam):
        Flounder.__init__(self)
        self.bam = bam
        self.samfile = pysam.AlignmentFile(bam, "rb")

    @functools.lru_cache()
    def get_bam_ranges(self, filter_flag=3844):
        logging.info("Extracting BAM ranges")
        # note that read_bam by default has a specific set of SAM flags
        # this code needs to be updated to allow for selection of +/-, primary
        # secondary etc ... - this method is also independent of pysam//samtools
        return pr.read_bam(self.bam, filter_flag=3844)

    def bam_stats(self, force=False):
        data = None
        if (not force) & ("flounder" in globals()):
            data = flounder.read_cache(
                self.bam, pd.DataFrame())
        if data is None:
            read_count = 0
            for l in pysam.idxstats(self.bam).splitlines():
                ll = l.split("\t")
                if ll[0] != "*":
                    read_count += int(ll[2])
            logging.debug(f"BamFile readcount == {read_count}")
            assignments = []
            # parse the BAM entries for info on primary, secondary, supplementary
            #    - pick out qualitative information for summarising mapping performance
            for read in tqdm(self.samfile.fetch(), total=read_count):
                flag = "1"
                if read.is_qcfail:
                    flag = "F"
                elif read.is_duplicate:
                    flag = "D"
                elif read.is_secondary:
                    flag = "2"
                elif read.is_supplementary:
                    flag = "S"
                elif read.is_unmapped:
                    flag = "U"
                readq = -10 * log10((10 ** (pd.Series(
                        read.query_alignment_qualities) / -10)).mean())
                datapoint = [flag, read.query_length, read.reference_length, readq, read.mapping_quality, read.reference_name, read.reference_start, read.reference_end]
                assignments.append(datapoint)
            # mung the datapoints into a DataFrame

            data = pd.DataFrame(
                assignments,
                columns=["flag", "query_length", "reference_length", "read_qual", "map_qual",
                         "Chromosome", "Start", "End"])

            if "flounder" in globals():
                flounder.write_cache(
                    self.bam, data)
        return data

    def bam_index_tiled_ranges(self, tile_size=1000):
        ranges = []
        for l in pysam.idxstats(self.bam).splitlines():
            ll = l.split("\t")
            if ll[0] != "*":
                ranges.append([ll[0], 0, ll[1]])
        tiled_ranges = pr.gf.tile_genome(
            pr.PyRanges(pd.DataFrame(ranges, columns=["Chromosome", "Start", "End"])),
            tile_size, tile_last=False).df
        return pr.PyRanges(tiled_ranges)

    def mapping_summary(self, tile_size=1000, long=False):
        mapping_data = self.bam_stats()
        tiled_ranges = self.bam_index_tiled_ranges(tile_size=tile_size).df
        tiled_ranges['cov'] = 0
        tiled_ranges = pr.PyRanges(tiled_ranges)

        def stratify_bam_coverage(dataframe, key):
            subdata = dataframe[dataframe.flag == key].copy()
            map_ranges = pr.PyRanges(
                subdata.loc[:, ["Chromosome", "Start", "End"]])
            bg_rle = tiled_ranges.to_rle("cov") # background of zero
            rle = map_ranges.to_rle() + bg_rle
            df = rle[tiled_ranges].df
            df['VR'] = df['Run'] * df['Value']
            df = df.groupby(["Chromosome", "Start"]).agg(
                {"Chromosome": "first", "Start": "first", "End": "first",
                 "Run": np.sum, "VR": np.sum})
            df['MeanCoverage'] = df['VR'] / df['Run']
            df = pr.PyRanges(
                df.reset_index(drop=True).drop(["Run", "VR"], axis=1))

            readq = -10 * log10((10 ** (pd.Series(subdata.read_qual) / -10)).mean())
            mapq = -10 * log10((10 ** (pd.Series(subdata.map_qual) / -10)).mean())
            return pd.Series({"reads": "{:,}".format(len(subdata.index)),
                             "bases": "{:,}".format(subdata.query_length.sum()),
                             "mapped_bases": "{:,}".format(subdata.reference_length.sum()),
                             "mean_length": "{:.2f}".format(subdata.query_length.mean()),
                             "mean_quality": "{:.2f}".format(readq),
                             "map_quality": "{:.2f}".format(mapq),
                             "mean_coverage": "{:.2f}".format(df.df.MeanCoverage.mean())})

        mapping_summary = pd.DataFrame.from_dict({"primary": stratify_bam_coverage(mapping_data, "1"),
                                      "secondary": stratify_bam_coverage(mapping_data, "2"),
                                      "supplementary": stratify_bam_coverage(mapping_data, "S")}, orient="columns")

        mapping_summary['secondary'][['bases', 'mean_quality', 'map_quality', 'mean_length']] = ["", "", "", ""]

        if long:
            mapping_summary.reset_index(inplace=True)
            mapping_summary = pd.melt(
                mapping_summary, id_vars=["index"],
                value_vars=["primary", "secondary", "supplementary"], var_name="map_type")
            arrays = [mapping_summary.map_type, mapping_summary['index']]
            mapping_summary = pd.DataFrame(mapping_summary['value'].tolist(), index=arrays)
            mapping_summary.set_axis([os.path.basename(self.bam)], axis=1, inplace=True)
        else:
            col = pd.MultiIndex.from_arrays(
                [[os.path.basename(self.bam), os.path.basename(self.bam), os.path.basename(self.bam)],
                 ["primary", "secondary", "supplementary"]])
            mapping_summary = pd.DataFrame(mapping_summary.values.tolist(), columns=col, index=mapping_summary.index)
        return mapping_summary


    def plot_coverage_distribution(self, tile_size=1000, bins=30, deepest_bin=None, **kwargs):

        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        mapping_data = self.bam_stats()
        tiled_ranges = self.bam_index_tiled_ranges(tile_size=tile_size).df
        tiled_ranges['cov'] = 0
        tiled_ranges = pr.PyRanges(tiled_ranges)

        subdata = mapping_data.copy()

        map_ranges = pr.PyRanges(
            subdata.loc[:, ["Chromosome", "Start", "End"]])
        bg_rle = tiled_ranges.to_rle("cov")  # background of zero
        rle = map_ranges.to_rle() + bg_rle
        df = rle[tiled_ranges].df
        df['VR'] = df['Run'] * df['Value']
        df = df.groupby(["Chromosome", "Start"]).agg(
            {"Chromosome": "first", "Start": "first", "End": "first",
             "Run": np.sum, "VR": np.sum})
        df['MeanCoverage'] = df['VR'] / df['Run']
        coverage = pr.PyRanges(
            df.reset_index(drop=True).drop(["Run", "VR"], axis=1))

        mc = coverage.MeanCoverage
        if deepest_bin is None:
            deepest_bin = mc.max() + 1
        boundaries = np.linspace(
            0, deepest_bin, num=bins, endpoint=True, retstep=False)
        assignments = np.digitize(mc, boundaries)

        xxx = pd.DataFrame({"coverage": mc,
                            "assignment": assignments,
                            "tally": tile_size}).groupby(
            ["assignment"]).agg(
            {"assignment": ["first"],
             "tally": [np.sum]})
        xxx.columns = xxx.columns.droplevel()
        yyy = pd.DataFrame({"start": boundaries[:-1],
                            "end": boundaries[1:]},
                           index=np.arange(1, len(boundaries)))
        coverage_dist = yyy.merge(
            xxx, how="outer", left_index=True, right_index=True)
        coverage_dist.columns = ["start", "end", "batch", "count"]
        coverage_dist["batch"] = coverage_dist.index
        coverage_dist = coverage_dist.fillna(0)
        coverage_dist["colour"] = "#1F78B4"

        p = figure(
            title="Histogram showing distribution of coverage",
            background_fill_color="lightgrey", plot_width=plot_width,
            plot_height=plot_height, tools=plot_tools)
        p.quad(
            source=coverage_dist, top="count", bottom=0, left='start',
            right='end', fill_color='colour', line_color="white", alpha=0.7)

        p.xaxis.axis_label = 'Depth-of-coverage (X-fold)'
        p.yaxis.axis_label = 'Bases of genome (n)'

        return self.handle_output(p, plot_type)



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
