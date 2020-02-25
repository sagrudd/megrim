from megrim.environment import MegrimPlugin
from megrim.reference_genome import ReferenceGenome
from megrim.genome_geometry import BamHandler, include_flounder
import pandas as pd
import argparse
import logging
import warnings


class BamStats(MegrimPlugin):
    def __init__(self):
        super().__init__()
        self.tool = "BamStats"

    def execute(self, args):
        include_flounder(args)

        warnings.simplefilter(action='ignore', category=FutureWarning)

        bam = BamHandler(args.bam)
        bam.mapping_summary()


    def arg_params(self, subparsers, parent_parser):
        argparser = subparsers.add_parser(self.tool, help="bam stats help", parents=[parent_parser])
        argparser.add_argument('-b', '--bam', metavar="/path/to/BAM", action='store', help='Path to the BAM-format mapping file',
                               required=True, dest="bam")
