from megrim.environment import MegrimPlugin
from megrim.base_modifications import fast5s_to_basemods, map_methylation_signal, include_flounder, reduce_mapped_methylation_signal, augment_reduced_methylation_signal
from megrim.reference_genome import ReferenceGenome
from megrim.genome_geometry import BamHandler
import pandas as pd
import argparse
import logging
import warnings


class BaseModifications(MegrimPlugin):
    def __init__(self):
        super().__init__()
        self.tool = "BaseModifications"

    def execute(self, args):
        include_flounder(args)

        warnings.simplefilter(action='ignore', category=FutureWarning)

        modifications = fast5s_to_basemods(args.fast5, modification=args.modifcation,
                                           threshold=args.probability, context=args.context)
        modifications.set_index("read_id", drop=False, inplace=True)
        # print out the modifications - quick reality check
        print(modifications)

        # define a reference genome and bam file
        reference = ReferenceGenome(args.fasta)
        bam = BamHandler(args.bam)

        pd.set_option("display.max_columns", None)
        pd.set_option("display.expand_frame_repr", False)
        pd.set_option("max_colwidth", -1)

        # associated mapped bases with the available modifications
        methylation_signal = map_methylation_signal(reference, bam, modifications)
        print(methylation_signal)

        reduced_reads = reduce_mapped_methylation_signal(methylation_signal)
        print(reduced_reads)

        augment_reduced_methylation_signal(reduced_reads, bam)
        # we should save this "result-file" as a deliverable

    def arg_params(self, subparsers, parent_parser):
        argparser = subparsers.add_parser(self.tool, help="base modifications help", parents=[parent_parser])
        argparser.add_argument('-5', '--fast5', metavar="/path/to/FAST5", action='store', help='Path to the FAST5-format sequences', required=True, dest="fast5")
        argparser.add_argument('-b', '--bam', metavar="/path/to/BAM", action='store', help='Path to the BAM-format mapping file',
                               required=True, dest="bam")
        argparser.add_argument('-f', '--fasta', metavar="/path/to/fasta", action='store', help='Path to the fasta format reference sequence',
                               required=True, dest="fasta")
        argparser.add_argument('-p', '--probability', metavar="[0..1]", action='store',
                               help='Base-modification probability. This is a floating point number between 0 and 1. A stringent selection will be closer to 1. [The default is 0.90]',
                               dest="probability", default=0.90, type=float)
        argparser.add_argument('-c', '--context', metavar="CG", action='store',
                               help='Base-modification context. Only CpG has been implemented at present. [The default is CG]',
                               dest="context", default="CG")
        argparser.add_argument('-m', '--modification', metavar="[5mC|6mA]", action='store',
                               help='The base modification to score for - this may be either 5mC or 6mA in this version of the software. [The default is 5mC]',
                               dest="modifcation", default="5mC")