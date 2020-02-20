from megrim.environment import MegrimPlugin
import argparse

class BaseModifications(MegrimPlugin):
    def __init__(self):
        super().__init__()
        self.tool = "BaseModifications"

    def execute(self, args):
        print(42)

    def arg_params(self, subparsers, parent_parser):
        argparser = subparsers.add_parser(self.tool, help="base modifications help", parents=[parent_parser])
        argparser.add_argument('-5', '--fast5', metavar="/path/to/FAST5", action='store', help='Path to the FAST5-format sequences', required=True, dest="fast5")
        argparser.add_argument('-b', '--bam', metavar="/path/to/BAM", action='store', help='Path to the BAM-format mapping file',
                               required=True, dest="bam")
        argparser.add_argument('-f', '--fasta', metavar="/path/to/fasta", action='store', help='Path to the fasta format reference sequence',
                               required=True, dest="fasta")
        argparser.add_argument('-p', '--probability', metavar="[0..1]", action='store',
                               help='Base-modification probability. This is a floating point number between 0 and 1. A stringent selection will be closer to 1.',
                               dest="probability", default=0.90, type=float)
        argparser.add_argument('-c', '--context', metavar="CpG", action='store',
                               help='Base-modification context. Only CpG has been implemented at present',
                               dest="context", default="CpG")
        argparser.add_argument('-m', '--modification', metavar="[5mC|6mA]", action='store',
                               help='The base modification to score for - this may be either 5mC or 6mA in this version of the software',
                               dest="modifcation", default="5mC")