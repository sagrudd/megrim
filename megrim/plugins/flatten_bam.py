"""
A plugin module to provide a command line interface to FASTQ seq streams.

This module provides a single plugin class to be used with the megrim
toolbox application. This is intended to enable streaming of FASTQ format
sequence streams from a variety of input formats - designed to assist and
simplify logic in Snakemake workflows
"""

from megrim.environment import MegrimPlugin
import gzip
import bz2
from mimetypes import guess_type
from functools import partial
import os
import sys
import warnings
from megrim.genome_geometry import BamHandler











class FlattenBam(MegrimPlugin):
    """ToFastq class - to provide fastq streaming functionality at cmdline."""

    def __init__(self):
        super().__init__()
        self.tool = "FlattenBam"
        self.args = None





    def execute(self, args):
        """
        Execute the plugin functionality.

        This method is a plugin requirement from the toolbox module - this
        orchestrates the logic that is contained within the function; in this
        case this means pulling fastq sequences from the provided file(s) and
        streaming sequence to stdout.

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
        warnings.simplefilter(action='ignore', category=FutureWarning)
        self.args = args
        
        print(f"Parsing BAM [{args.bam}] at depth @{args.depth}")
        bam = BamHandler(args.bam)
        
        for chunk in bam.chunk_generator():
            print(chunk)



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
            self.tool, help="FlattenBam help", parents=[parent_parser])
        argparser.add_argument(
            '-b', '--bam', metavar="/path/to/bam_file", action='store',
            help='Path to the source BAM file', required=True, dest="bam")
        argparser.add_argument(
            '-d', '--depth', action='store', help='Depth to filter to',
            required=True, dest="depth", default=False, type=int)
        



if __name__ == '__main__':
    print("Hellow world")