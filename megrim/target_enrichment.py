#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 18:09:10 2020

@author: srudd
"""


from megrim.environment import Flounder
from megrim.genome_geometry import BamHandler, BedHandler
from megrim.reference_genome import ReferenceGenome, augment_annotation
import logging
import functools


class TargetEnrichment(Flounder):

    def __init__(self, reference_genome=None, bed_file=None, bam_file=None,
                 target_proximity=5000, background_threshold=20):
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
        
        
        
    @functools.lru_cache()
    def get_on_target(self, tile_size=100):
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
    def get_on_target_annotation(self, tile_size=100):
        logging.info("loading on target annotation")
        return augment_annotation(
            self.bam, self.get_on_target(tile_size=tile_size))
        
    @functools.lru_cache()
    def get_target_proximal(self, tile_size=100):
        logging.info("loading target proximal coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_target_proximal_ranges())
    
    @functools.lru_cache()
    def get_target_proximal_annotation(self, tile_size=100):
        logging.info("loading target proximal annotation")
        return augment_annotation(
            self.bam, self.get_target_proximal(tile_size=tile_size))

    @functools.lru_cache()
    def get_off_target(self, tile_size=100):     
        logging.info("loading off target coordinates")
        filtered_coverage = self.ref.get_tiled_mean_coverage(
            self.bam, tile_size=tile_size).subtract(
                self.get_on_target(tile_size=tile_size))
        return self.ref.get_tiled_mean_coverage(
            self.bam, 
            ranges=filtered_coverage[
                filtered_coverage.MeanCoverage >= \
                    self.background_threshold].slack(10).merge().slack(-10))
     
    @functools.lru_cache()    
    def get_off_target_annotation(self, tile_size=100):    
        logging.info("loading off target annotation")
        return augment_annotation(
            self.bam, self.get_off_target(tile_size=tile_size))
        
    @functools.lru_cache()
    def get_background(self, tile_size=100):
        logging.info("loading background coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_untargeted_ranges().subtract(self.get_off_target()))
    
    @functools.lru_cache()
    def get_background_annotation(self, tile_size=100):
        logging.info("loading background annotation")
        return augment_annotation(
            self.bam, self.get_background(tile_size=tile_size))
        
    def get_fine_coverage_aggregate(self, universe, target_proximity=None):
        if target_proximity is None:
            target_proximity = self.target_proximity
        aggregated_cov = self.ref.deep_dive(
            self.bam, universe, target_proximity=target_proximity)
        return aggregated_cov

    