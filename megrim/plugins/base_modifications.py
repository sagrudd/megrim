from megrim.environment import MegrimPlugin
from megrim.base_modifications import BaseModifications
from megrim.reference_genome import ReferenceGenome
from megrim.genome_geometry import BamHandler
import pandas as pd
import logging
import warnings
import os
import multiprocessing


class BaseModificationsPlugin(MegrimPlugin):
    def __init__(self):
        super().__init__()
        self.tool = "BaseModifications"

    def execute(self, args):
        warnings.simplefilter(action='ignore', category=FutureWarning)
        os.environ["NUMEXPR_MAX_THREADS"] = str(multiprocessing.cpu_count())
        fast5 = args.fast5
        bam = BamHandler(args.bam, args)
        reference = ReferenceGenome(args.fasta)
        base_mods = BaseModifications(fast5, bam, reference, args)
        if args.index:
            logging.degug(f"saving base-mod coordinates to CSV file [{args.output}]")
            base_mods.fast5s_to_basemods().to_csv(args.output, sep="\t")
        else:
            logging.debug(f"saving data as CSV file [{args.output}]")
            base_mods.reduce_mapped_methylation_signal().to_csv(args.output, sep="\t", index=False, chunksize=1e6)
            # use the chunksize here = from default (None) to 1e6 reduces run time by ~ 15X
        logging.debug(f"fin ...")


    def arg_params(self, subparsers, parent_parser):
        argparser = subparsers.add_parser(self.tool, help="base modifications help", parents=[parent_parser])
        argparser.add_argument('-5', '--fast5', metavar="/path/to/FAST5", action='store',
                               help='Path to the FAST5-format sequences', required=True, dest="fast5")
        argparser.add_argument('-b', '--bam', metavar="/path/to/BAM", action='store',
                               help='Path to the BAM-format mapping file', required=True, dest="bam")
        argparser.add_argument('-f', '--fasta', metavar="/path/to/fasta", action='store',
                               help='Path to the fasta format reference sequence', required=True, dest="fasta")
        argparser.add_argument('-p', '--probability', metavar="[0..1]", action='store',
                               help='Base-modification probability. This is a floating point number between 0 and 1. '
                                    'A stringent selection will be closer to 1. [The default is 0.90]',
                               dest="probability", default=0.90, type=float)
        argparser.add_argument('-c', '--context', metavar="CG", action='store',
                               help='Base-modification context. Only CpG has been implemented at present. '
                                    '[The default is CG]', dest="context", default="CG")
        argparser.add_argument('-m', '--modification', metavar="[5mC|6mA]", action='store',
                               help='The base modification to score for - this may be either 5mC or 6mA in this '
                                    'version of the software. [The default is 5mC]', dest="modifcation", default="5mC")
        argparser.add_argument('-x', '--index', action='store_true', dest="index",
                               help="only index the FAST5; do not process the bam files.", default=False)
        argparser.add_argument('-o', '--output', metavar="results-file", action='store', dest="output", required=True,
                               help="file path to a file location where the results will be stored. The results will "
                                    "be stored in a TSV format.")