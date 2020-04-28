
from megrim.environment import MegrimPlugin
from megrim.conda_reader import CondaGit
import warnings
import multiprocessing
import os


class Conda2Rpm(MegrimPlugin):
    """Class for mapping conda specs into an RPM framework"""

    def __init__(self):
        super().__init__()
        self.tool = "BamStats"

    def execute(self, args):
        warnings.simplefilter(action='ignore', category=FutureWarning)
        os.environ["NUMEXPR_MAX_THREADS"] = str(multiprocessing.cpu_count())

        conda = CondaGit(args.bam)

        conda.lookup()

    def arg_params(self, subparsers, parent_parser):

        argparser = subparsers.add_parser(
            self.tool, help="conda2rpm help", parents=[parent_parser])

        argparser.add_argument(
            '-s', '--source', metavar="bioconda git repo", action='store',
            help='Path to the git cloned content from bioconda', required=True,
            dest="bioconda")
        argparser.add_argument(
            '-t', '--target', metavar="target method", action='store',
            dest="target", required=True, help='name of the package to process')
