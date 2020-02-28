"""
A plugin module to provide a command line interface to BAM statistics.

This module provides a single plugin class to be used with the megrim
toolbox application. This is intended to present core mapping statistics
from the BAM class of the genome_geometry plugin. Please see the toolbox
module for information on the plugin requirements.
"""

from megrim.environment import MegrimPlugin
from megrim.genome_geometry import BamHandler, include_flounder
import warnings
import multiprocessing
import os


class BamStats(MegrimPlugin):
    """BamStats class - to provide basic BAM mapping statistics at cmdline."""

    def __init__(self):
        super().__init__()
        self.tool = "BamStats"

    def execute(self, args):
        """
        Execute the plugin functionality.

        This method is a plugin requirement from the toolbox module - this
        orchestrates the logic that is contained within the function; in this
        case this means pulling some basic statistics from the BAM file and
        returning these in a pandas format, printed to console.

        Parameters
        ----------
        args: argparse derived object
            The requirement is an argparse object containing the minimally
            required parameters and optional parameters for a run-through of
            the workflow.

        Returns
        -------
        Nothing at all - stuff may be presented to screen.

        """
        include_flounder(args)

        warnings.simplefilter(action='ignore', category=FutureWarning)
        os.environ["NUMEXPR_MAX_THREADS"] = str(multiprocessing.cpu_count())

        bam = BamHandler(args.bam)
        map_stats = bam.mapping_summary()
        map_stats.to_csv(args.output, sep="\t")

    def arg_params(self, subparsers, parent_parser):
        """
        Append an argparse subparser object with salient parameters.

        This method is another requirement from the toolbox plugin framework.
        The method appends workflow specific parameters to the generic
        toolbox parameters. This is used to customise and parameterise
        workflows as required.

        Parameters
        ----------
        subparsers: argparse
            The argparse object to which the parameters will be added.
        parent_parser: argparse
            The parental argparse object (which at the time of implementation)
            seemed required for functionality.

        Returns
        -------
        None. - this stuff is appended to object.
        """
        argparser = subparsers.add_parser(
            self.tool, help="bam stats help", parents=[parent_parser])

        argparser.add_argument(
            '-b', '--bam', metavar="/path/to/BAM", action='store',
            help='Path to the BAM-format mapping file', required=True,
            dest="bam")
        argparser.add_argument(
            '-o', '--output', metavar="results-file", action='store',
            dest="output", required=True, help='file path to a file location '
            'where the results will be stored. The results will be stored in '
            'a TSV format.')
