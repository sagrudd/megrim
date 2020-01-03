#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:40:04 2020

@author: srudd
"""


from environment import Flounder
from genome_geometry import BamHandler, BedHandler
from reference_genome import ReferenceGenome

# create an instance of flounder for path handling
flounder = Flounder()
flounder.set_path("/tmp/floundeR")

# define a reference genome - this requires a fasta file
ref = ReferenceGenome("/Users/srudd/Desktop/Human_genome.fasta")
ref.info()
ref.skip_chromosome("MT")

# define a bed file of target regions of interest
bed = BedHandler("/Users/srudd/Desktop/enrichment_targets.bed")
bed.set_reference(ref)
bed.set_target_proximity(5000)
print(bed.get_bed_targets())

# use the bed information to create pyranges coordinates for on-target,
# target-proximal and not-target regions of the genome 
on_target = bed.get_target_ranges()
print(on_target)
target_proximal = bed.get_target_proximal_ranges()
print(target_proximal)
untargeted = bed.get_untargeted_ranges()
print(untargeted)

# define the BAM file of interest
bam = BamHandler("/Users/srudd/Desktop/cas9_FAK76554.bam")

