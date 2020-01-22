#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 18:09:10 2020

@author: srudd
"""


from megrim.environment import Flounder
from megrim.genome_geometry import BamHandler, BedHandler
from megrim.reference_genome import ReferenceGenome, augment_annotation, weighted_percentile
from megrim.infographic_plots import InfographicPlot, InfographicNode
import numpy as np
import pandas as pd
import logging
import functools
import pyranges as pr
from bokeh.plotting import figure
from bokeh.models import Span, NumeralTickFormatter


class TargetEnrichment(Flounder):

    def __init__(self, reference_genome=None, bed_file=None, bam_file=None,
                 target_proximity=5000, background_threshold=20, 
                 tile_size=100):
        Flounder.__init__(self)
        
        logging.info("Preparing reference genome")
        self.ref = ReferenceGenome(reference_genome)
        
        logging.info("Preparing BED coordinates")
        self.bed = BedHandler(bed_file)
        self.bed.set_reference(self.ref)
        self.bed.set_target_proximity(target_proximity)
        
        logging.info("setting BAM file")
        self.bam = BamHandler(bam_file)
        self.background_threshold = background_threshold
        self.target_proximity = target_proximity
        self.tile_size = 100
        
        
        
    @functools.lru_cache()
    def get_on_target(self):
        """
        on-target is defined in the host bed file

        Parameters
        ----------
        tile_size : TYPE, optional
            DESCRIPTION. The default is 100.

        Returns
        -------
        on_target_universe : TYPE
            DESCRIPTION.

        """
        logging.info("loading on target coordinates")
        on_target_universe = self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_target_ranges())
        on_target_universe.Name = self.bed.get_target_ranges().Name
        return on_target_universe
    
    @functools.lru_cache()
    def get_on_target_annotation(self):
        logging.info("loading on target annotation")
        return augment_annotation(
            self.bam, self.get_on_target())
        
    @functools.lru_cache()
    def get_target_proximal(self):
        logging.info("loading target proximal coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_target_proximal_ranges())
    
    @functools.lru_cache()
    def get_target_proximal_annotation(self):
        logging.info("loading target proximal annotation")
        return augment_annotation(
            self.bam, self.get_target_proximal())

    @functools.lru_cache()
    def get_off_target(self, tile_size=None):    
        if tile_size is None:
            tile_size = self.tile_size
        logging.info(
            "loading off target coordinates (tile_size={})".format(tile_size))
        filtered_coverage = self.ref.get_tiled_mean_coverage(
            self.bam, tile_size=tile_size).subtract(
                self.get_on_target()).subtract(self.get_target_proximal())
        return self.ref.get_tiled_mean_coverage(
            self.bam, 
            ranges=filtered_coverage[
                filtered_coverage.MeanCoverage >= \
                    self.background_threshold * self.get_background_coverage()].slack(10).merge().slack(-10))
    
    @functools.lru_cache()
    def get_background_coverage(self, tile_size=None):
        if tile_size is None:
            tile_size = self.tile_size
        logging.info(
            "calculating background coverage (tile_size={})".format(tile_size))
        filtered_coverage = self.ref.get_tiled_mean_coverage(
            self.bam, tile_size=tile_size).subtract(
                self.get_on_target()).subtract(self.get_target_proximal())
        # this is an approximation since not all windows will be the same size
        mean_coverage = filtered_coverage.MeanCoverage.mean()
        return mean_coverage
            
    @functools.lru_cache()    
    def get_off_target_annotation(self, tile_size=None):    
        logging.info("loading off target annotation")
        return augment_annotation(
            self.bam, self.get_off_target(tile_size=tile_size))
        
    @functools.lru_cache()
    def get_background(self):
        logging.info("loading background coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_untargeted_ranges().subtract(self.get_off_target()))
    
    @functools.lru_cache()
    def get_background_annotation(self):
        logging.info("loading background annotation")
        return augment_annotation(
            self.bam, self.get_background())
        
    def get_fine_coverage_aggregate(self, universe, target_proximity=None):
        if target_proximity is None:
            target_proximity = self.target_proximity
        aggregated_cov = self.ref.deep_dive(
            self.bam, universe, target_proximity=target_proximity)
        return aggregated_cov
    
    def get_mapped_bases(self): 
        return self.get_off_target_annotation().df.bases_start.sum() + \
            self.get_background_annotation().df.bases_start.sum() + \
            self.get_target_proximal_annotation().df.bases_start.sum() + \
            self.get_on_target_annotation().df.bases_start.sum()
            
    def get_on_target_perc(self):
        return self.get_on_target_annotation().df.rstart.sum() / \
            (self.get_on_target_annotation().df.rstart.sum() + \
             self.get_off_target_annotation().df.rstart.sum() + \
                 self.get_background_annotation().df.rstart.sum() + \
                     self.get_target_proximal_annotation().df.rstart.sum()) * 100
                
    def get_depletion_factor(self):
        # depletion_label = "FAILED"
        return weighted_percentile(self.get_on_target_annotation()) / \
            weighted_percentile(self.get_background_annotation())
            
    def get_mean_target_coverage(self):
        return self.get_on_target_annotation().df.MeanCoverage.repeat(
            self.get_on_target_annotation().End - \
                self.get_on_target_annotation().Start).mean()
            
    def get_cas9_exec_infographic(self, **kwargs):
        
        (plot_width, plot_dpi) = self.handle_kwargs(
            ["plot_width", "plot_dpi"], **kwargs)
                
        throughput = InfographicNode(legend="Throughput",
                                        value="{:.2f} Gb".format(self.get_mapped_bases()/1e9),
                                        graphic='calculator')
        
        on_target_node = InfographicNode(legend="Reads on Target",
                                         value="{:.2f}%".format(self.get_on_target_perc()),
                                         graphic='cut')
        
        coverage_node = InfographicNode(legend="Mean Target Coverage",
                                     value="{:.2f}X".format(self.get_mean_target_coverage()),
                                     graphic='map-marked')
        
        depletion_node = InfographicNode(legend="Non-target Depletion",
                                     value="{:.1f}X".format(self.get_depletion_factor()),
                                     graphic='fill-drip')
        
        infographic_data = [throughput, on_target_node, coverage_node, depletion_node]
        ip = InfographicPlot(infographic_data, rows=1, columns=4)
        ip.plot_infographic(plot_width, plot_dpi)

        
    def get_mapping_stats(self, annotated_ranges, scale="Megabases"):
        
        scaleVal = 1
        if scale == "Gigabases":
            scaleVal = 1e9
        elif scale == "Megabases":
            scaleVal = 1e6
        elif scale == "Kilobases":
            scaleVal = 1e3
        
        total_reads = annotated_ranges.df.rstart.sum()
        fastq_nts = annotated_ranges.df.bases_start.sum() / scaleVal
        mapped_clipped = annotated_ranges.df.cigar_m.sum() / scaleVal
        genome_size = self.ref.get_genome_size()
        mapped_space = np.sum(annotated_ranges.End - annotated_ranges.Start)
        ref_space = mapped_space / genome_size * 100
        mean_cov = np.sum((annotated_ranges.End - annotated_ranges.Start) * \
            annotated_ranges.MeanCoverage) / \
            np.sum(annotated_ranges.End - annotated_ranges.Start)
        # return a pd.Series of observations
        return pd.Series({"Total reads": "{:.0f}".format(total_reads),
                          scale: "{:.2f}".format(fastq_nts),
                          "Mapped {}".format(scale): "{:.2f}".format(mapped_clipped),
                          "Genome fraction (%)": "{:.3f}".format(ref_space),
                          "Mean coverage (X)": "{:.2f}".format(mean_cov)})
    
    def get_mapping_summary_stats(self):
        on_t = self.get_mapping_stats(self.get_on_target_annotation())
        off_t = self.get_mapping_stats(self.get_off_target_annotation())
        t_p = self.get_mapping_stats(self.get_target_proximal_annotation())
        bg = self.get_mapping_stats(self.get_background_annotation())
        
        df = pd.concat([on_t, t_p, off_t, bg], axis=1, keys=["on_target", "target_proximal", "off_target", "background"])
        df.iloc[0] = df.iloc[0].astype(int)
        return df
            
    
    def get_target_performance(self, target, scale="Megabases"):
        scaleVal = 1
        if scale == "Gigabases":
            scaleVal = 1e9
        elif scale == "Megabases":
            scaleVal = 1e6
        elif scale == "Kilobases":
            scaleVal = 1e3
        df = self.get_on_target_annotation().df
        df = df.loc[df.Name==target,]
        
        return target, \
            "{:,}".format(int(int(df.End) - int(df.Start))), \
            "{:.2f}".format(float(df.MeanCoverage)), \
            "{:.0f}".format(float(df.rstart)), \
            "{:.2f}".format(float(df.bases_start) / scaleVal), \
            "{:.1f}".format(float(df.start_read_len)), \
            "{:.2f}".format(float(df.readq)), \
            "{:.2f}".format(float(df.mapq)), \
            "{:.1f}".format(float(df.strand_p / (df.strand_n + df.strand_p) * 100))    
        
    def get_target_performance_summary(self, scale="Megabases"):
        targets = self.get_on_target().Name
        return targets.apply(self.get_target_performance).apply(
            pd.Series, 
            index=["Target", "Target size", "Mean coverage", "Read count", 
                   scale, "MeanReadLen", "MeanReadQ", "MeanMapQ", "(+)Strand"])
        
        
        
    def get_target_plot(self, target, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        targets = self.get_on_target().df
        targets = targets.loc[targets.Name == target,]
        
        aggregated_cov = self.ref.deep_dive(
            bam=self.bam,
            arange=self.get_on_target(),
            target=target)

        print(aggregated_cov)
        
        plot = figure(
            title="Plot showing depth of coverage across ({}) target region".format(targets.Name.tolist()[0]), 
            x_axis_label="Position on Chr({})".format(str(targets.Chromosome.tolist()[0])),
            y_axis_label='Depth of coverage (X)', 
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)
        
        plot.line(aggregated_cov.Start, aggregated_cov.MeanCoverage, 
                  line_width=2, line_color='#1F78B4', 
                  legend_label='Depth of Coverage')
        
        start_line = Span(location=targets.Start.tolist()[0], dimension='height', line_color='red', line_width=2, line_alpha=0.7)
        end_line = Span(location=targets.End.tolist()[0], dimension='height', line_color='red', line_width=2, line_alpha=0.7)
        bg_line = Span(location=self.get_background_coverage(), dimension='width', line_color='orange', line_width=2, line_alpha=0.7)
        
        plot.renderers.extend([start_line, end_line, bg_line])
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)
        
    