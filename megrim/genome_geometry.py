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
    """
    Set the working Flounder environment for methods in this module.

    Flounder provides a collection of tools and environment handles for the
    caching of results and the orchestration of bokeh figure preparation.
    Use this module level function to augment the environment for calls
    outside of a Flounder subclass.

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # setup a Flounder for this workflow ...
    global flounder
    flounder = Flounder()
    if isinstance(args, argparse.Namespace):
        flounder.argparse(args)
    if isinstance(args, dict):
        flounder.dictparse(args)


class BedHandler:
    """
    BedHandler provides a class for quick parsing of BED coordinates.

    The functionality within the class is currently driven by tutorial
    requirements and everything is likely to change.
    """

    def __init__(self, bedfile):
        self.bedfile = bedfile
        self.targets = None
        self.proximity = 0
        self.ref = None
        self.load_bed_file()

    def load_bed_file(self):
        """
        Load the class associated bed file.

        This loads the class associated bed file into memory.

        Returns
        -------
        None.

        """
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
        return pr.gf.genome_bounds(
            merged, self.ref.get_reference_ranges(), clip=True)

    def get_untargeted_ranges(self):
        whole_genome_ranges = self.ref.get_reference_ranges()

        target_proximal_ranges = self.get_target_ranges()
        target_proximal_ranges.Start += -(self.get_target_proximity() + 1)
        target_proximal_ranges.End += self.get_target_proximity() + 1
        # there is a possibility for edge cases where boundaries of
        # target_proximal_ranges extend beyond the limits of the chromosome
        untargeted = whole_genome_ranges.subtract(target_proximal_ranges)
        return pr.gf.genome_bounds(
            untargeted, self.ref.get_reference_ranges(), clip=True)


class BamHandler(Flounder):

    def __init__(self, bam, args=None):
        Flounder.__init__(self)
        self.bam = bam
        self.samfile = pysam.AlignmentFile(bam, "rb")
        if args is not None:
            if isinstance(args, argparse.Namespace):
                self.argparse(args)
            if isinstance(args, dict):
                self.dictparse(args)

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
        """
        Plot a histogram showing frequency of different depths-of-coverage.

        Parameters
        ----------
        tile_size : TYPE, optional
            DESCRIPTION. The default is 1000.
        bins : TYPE, optional
            DESCRIPTION. The default is 30.
        deepest_bin : TYPE, optional
            DESCRIPTION. The default is None.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

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

        p.xaxis.axis_label = 'Depth-of-coverage (X)'
        p.yaxis.axis_label = 'Frequency (bases)'
        p.yaxis.formatter = NumeralTickFormatter(format="0,0")
        p.grid.grid_line_color = "white"

        return self.handle_output(p, plot_type)

    def chunk_generator(self, tile_size=5000000, force=False):
        """
        Prepare a chunk generator of BAM content across whole genome.

        For whole-genome analysis, it makes sense to process BAM files in
        chunks of content - this can enable both memory tractability for
        tutorials and can enable robust parallelisation for other scenarios.
        This method will consume the BAM file and emit chunks of salient
        content for parsing by calling functions.

        Parameters
        ----------
        tile_size : TYPE, optional
            DESCRIPTION. The default is 5000000.
        force : TYPE, optional
            DESCRIPTION. The default is False.

        Yields
        ------
        bam_chunk: pandas.DataFrame
            DataFrame where each row corresponds to a mapped sequence entry.

        """
        map = self.bam_index_tiled_ranges(tile_size=tile_size)
        for i in map.df.index:
            chromo = map.df.iloc[i].Chromosome
            start = map.df.iloc[i].Start
            end = map.df.iloc[i].End
            logging.info(
                f"extracting reads from chromosome {chromo} chunk {i+1}/{len(map.df.index)} [{start}:{end}]")
            bam_chunk = self.extract_bam_chunk(chromo, start, end, force)
            yield bam_chunk

    def extract_bam_chunk(self, chromosome, start, end, force=False):
        """
        Extract minimal bam associated annotation for given coordinates.

        This method will parse bamfile, bam, and return a pandas DataFrame
        containing key observations to facilitate a base-modification workflow.

        Parameters
        ----------
        bam: bam
            An object of class bam.
        chromosome: Str
            The chromosome object of interest.
        start: int
            start position in chromosome.
        end: int
            The chromosomal end position
        force: boolean, optional
            Whether to force the analysis, or whether to allow for a cached
            result to be returned. The default is False.

        Returns
        -------
        data: pd.DataFrame
            Pandas DataFrame containing the parsed entries.

        """
        # recover the cached data if possible
        if not force:
            data = self.read_cache(
                self.bam, pd.DataFrame(), chromosome, start, end)

        if data is None:
            data = []
            # extract reads from a bam file
            reads = self.get_sam_core(chromosome, start, end)
            for read in reads:
                # select primary mappings only
                if not any([read.is_secondary, read.is_supplementary,
                            read.is_qcfail, read.is_duplicate]):
                    row = [read.query_name, read.query_length,
                           read.reference_name, read.reference_start,
                           read.reference_end, "-" if read.is_reverse else "+",
                           read.cigartuples]
                    data.append(row)
            # convert data into a DataFrame
            data = pd.DataFrame(
                data,
                columns=["query_name", "query_length", "reference_name",
                         "reference_start", "reference_end", "strand",
                         "cigar"])
            data['block_start'] = start
            data['block_end'] = end
            # prettify the data
            data.set_index("query_name", drop=False, inplace=True)
            self.write_cache(self.bam, data, chromosome, start, end)
        return data

    def get_bam_coverage(self):
        """
        Return BAM coverage information in PyRanges format.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        logging.debug("Extracting BAM coverage")
        return self.get_bam_ranges().to_rle(strand=False)

    def get_sam_reads(self, chromosome, start, end):
        """
        Get SAM sequence reads from defined genomic interval.

        Parameters
        ----------
        chromosome : TYPE
            DESCRIPTION.
        start : TYPE
            DESCRIPTION.
        end : TYPE
            DESCRIPTION.

        Returns
        -------
        result : TYPE
            DESCRIPTION.

        """
        # logging.info("getSamReads ({}) {}:{}".format(chromosome, start, end))
        result = None
        try:
            result = self.samfile.fetch(chromosome, int(start), int(end))
        except Exception:
            result = None
        return result

    def get_sam_core(self, chromosome, start, end):
        """
        Get SAM core.

        Parameters
        ----------
        chromosome : TYPE
            DESCRIPTION.
        start : TYPE
            DESCRIPTION.
        end : TYPE
            DESCRIPTION.

        Returns
        -------
        annot : TYPE
            DESCRIPTION.

        """
        annot = self.get_sam_reads(chromosome, start, end)
        return annot

    def get_sam_annotation(self, chromosome, start, end):
        """
        Get SAM annotation.

        Parameters
        ----------
        chromosome : TYPE
            DESCRIPTION.
        start : TYPE
            DESCRIPTION.
        end : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
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
