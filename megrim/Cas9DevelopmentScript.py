#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:40:04 2020

@author: srudd
"""

import sys
from IPython.display import Image
from megrim.environment import tutorial_branding, Flounder
from megrim.genome_geometry import BamHandler, BedHandler
from megrim.reference_genome import ReferenceGenome, augment_annotation
from megrim.reproducible_research import SessionInfo
from importlib import reload
import logging
reload(logging)
logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',  
    level=logging.INFO, datefmt='%I:%M:%S')


tutorial_branding("Cas9", "Oxford Nanopore Tutorial: Cas9")

target_proximity = 5000

# define a reference genome - this requires a fasta file
logging.info("Preparing reference genome")
ref = ReferenceGenome("/Users/srudd/Desktop/Human_genome.fasta")
ref.skip_chromosome("MT")

# define a bed file of target regions of interest
logging.info("Preparing BED coordinates")
bed = BedHandler("/Users/srudd/Desktop/enrichment_targets.bed")
bed.set_reference(ref)
bed.set_target_proximity(target_proximity)

# use the bed information to create pyranges coordinates for on-target,
# target-proximal and not-target regions of the genome 
on_target = bed.get_target_ranges()
target_proximal = bed.get_target_proximal_ranges()
untargeted = bed.get_untargeted_ranges()

# define the BAM file of interest
logging.info("setting BAM file")
bam = BamHandler("/Users/srudd/Desktop/cas9_FAK76554.bam")

# bam.get_sam_annotation('1', 155179779, 155195266)

# create a tiled_genome_representation of coverage
tiled_coverage_means = ref.get_tiled_mean_coverage(bam, tile_size=100)
print(tiled_coverage_means)

# prepare coverage update for the on_target ranges
on_target_universe = ref.get_tiled_mean_coverage(bam, ranges=on_target)
on_target_universe.Name = on_target.Name
target_proximal_universe = ref.get_tiled_mean_coverage(
    bam, ranges=target_proximal)
print(target_proximal_universe)

# look for the off-target regions of the genome
off_target_scale = 20
# we have local depths of coverage within the tiled_coverage_means ...
# this also includes the on_target genomic regions ... first step is thus
# to filter out the regions of the genome that are on_target ...
filtered_coverage = tiled_coverage_means.subtract(on_target_universe)
background_threshold = filtered_coverage.MeanCoverage.mean() * off_target_scale
off_target_universe = ref.get_tiled_mean_coverage(bam, ranges=filtered_coverage[
    filtered_coverage.MeanCoverage >= background_threshold].slack(10).merge().slack(-10))
background_universe = ref.get_tiled_mean_coverage(bam, ranges=untargeted.subtract(off_target_universe))

on_target_universe = augment_annotation(bam, on_target_universe)
off_target_universe = augment_annotation(bam, off_target_universe)
background_universe = augment_annotation(bam, background_universe)
target_proximal_universe = augment_annotation(bam, target_proximal_universe)
aggregated_cov = ref.deep_dive(bam, on_target_universe, target_proximity=target_proximity)

print(SessionInfo())