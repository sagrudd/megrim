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
import logging
import functools


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
                self.get_on_target())
        return self.ref.get_tiled_mean_coverage(
            self.bam, 
            ranges=filtered_coverage[
                filtered_coverage.MeanCoverage >= \
                    self.background_threshold].slack(10).merge().slack(-10))
     
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

        

    