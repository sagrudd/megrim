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


FASTQ = "FASTQ"
FQ = "FQ"
FASTA = "FASTA"
FA = "FA"
BAM = "BAM"
SAM = "SAM"
FAST5 = "FAST5"
DIRECTORY = "DIRECTORY"
candidates = {FASTQ: FASTQ,  FQ: FASTQ, FASTA: FASTA, FA: FASTA, BAM: BAM,
              SAM: SAM, FAST5: FAST5}


def has_compress_suffix(file):
    """
    Check is the file gives any indication that it might be compressed.

    Simple method ... does a file have e.g. a gzip suffix?

    Parameters
    ----------
    file: str
        Filename of the file that might be compressed.

    Returns
    -------
    bool
        Is the file potentially compressed using a (bioinformatics) standard
        compression technique?

    """
    if file.endswith(".GZIP"):
        return True
    elif file.endswith(".GZ"):
        return True
    elif file.endswith(".BZ2"):
        return True
    elif file.endswith(".BZIP2"):
        return True
    return False


def strip_compress_suffix(file):
    """
    Remove the compression-format suffix from a file.

    Parameters
    ----------
    file: str
        Filename of the file that might be compressed.

    Returns
    -------
    str
        The filename stripped of canonical compression suffixes.

    """
    if file.endswith(".GZIP"):
        return file[:-5]
    if file.endswith(".GZ"):
        return file[:-3]
    elif file.endswith(".BZ2"):
        return file[:-4]
    elif file.endswith(".BZIP2"):
        return file[:-6]
    return file


def process_fastq(file):
    """
    Stream output from a fastq format input file.

    This method is intended to stream fastq format sequence data from a
    variety of input sources; this is the handler to process sequence from
    a fastq source file (that may be compressed).

    Parameters
    ----------
    file: str
        The input fastq file.

    Returns
    -------
    None - but a load of sequence should be streamed to the stdin.

    """
    encoding = guess_type(file)[1]
    # print(f"file encoding checked == {encoding}")
    _open = open
    if encoding == "gzip":
        _open = partial(gzip.open, mode="rt")
    if encoding == "bzip2":
        _open = partial(bz2.open, mode="rt")
    with _open(file) as f:
        # for record in SeqIO.parse(f, 'fastq'):
        #    sys.stdout.write(record.format("fastq"))
        for line in f:
            sys.stdout.write(line)


def process_file(file, file_type):
    """
    Sequence handler to process a validated sequence containing input file.

    Parameters
    ----------
    file: str
        Path to the input file.
    file_type: str
        An internally defined file-type for processing. This value is used
        to determine the handler that will be used.

    Returns
    -------
     None - but a load of sequence should be streamed to the stdin.

    """
    if file_type == FASTQ:
        process_fastq(file)


class ToFastq(MegrimPlugin):
    """ToFastq class - to provide fastq streaming functionality at cmdline."""

    def __init__(self):
        super().__init__()
        self.tool = "ToFastq"
        self.args = None

    def guess_format(self, file=None, iter=False):
        """
        Guess the sequence format for the provided sequence file.

        This method tries to guess the file sequence format on the basis of
        typical bioinformatics file formats.

        Parameters
        ----------
        file: str, optional
            The path to the sequence file to process. The default is None, but
            the None value should only be used for directory parsing purposes.
        iter: boolean, optional
            Whether if is a directory setting the files within directory
            should be iterated over to batch processing of multiple files.
            The default is False.

        Raises
        ------
        ValueError
            This will be thrown if the file does not correspond to a valid
            bioinformatics sequence format.

        Returns
        -------
        str
            The sequence format of the test file..

        """
        if os.path.isdir(self.args.src) & self.args.dir & (file is None):
            return DIRECTORY
        if file is None:
            file = os.path.basename(self.args.src).upper()
        if has_compress_suffix(file):
            file = strip_compress_suffix(file)
        else:
            file = file.upper()
        for key in candidates.keys():
            if file.endswith(key):
                return candidates[key]
        if iter:
            return None
        raise ValueError("file does not contain an acceptable sequence format")

    def parse_directory(self):
        """
        Parse available filesystem directory for sequence files.

        Returns
        -------
        None.

        """
        for file in os.listdir(self.args.src):
            file_type = self.guess_format(file, True)
            if file_type is not None:
                if (self.args.pickup is None) or (
                        file_type == self.args.pickup.upper()):
                    process_file(os.path.join(self.args.src, file), file_type)

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

        file_format = self.guess_format()
        if file_format == DIRECTORY:
            self.parse_directory()
        else:
            try:
                process_file(self.args.src, file_format)
            finally:
                try:
                    sys.stdout.flush()
                finally:
                    try:
                        sys.stdout.close()
                    finally:
                        try:
                            sys.stderr.flush()
                        finally:
                            sys.stderr.close()

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
            self.tool, help="tofastq help", parents=[parent_parser])
        argparser.add_argument(
            '-s', '--source', metavar="/path/to/input_file", action='store',
            help='Path to the source sequence file', required=True, dest="src")
        argparser.add_argument(
            '-d', '--directory', action='store_true', help='Process a folder '
            'of files rather than a monolith',
            required=False, dest="dir", default=False)
        argparser.add_argument(
            '-p', '--pickup', action='store', help='For directory parsing, '
            'specify filetype to pickup selectively pickup {FASTQ|BAM}',
            required=False, dest="pickup", default=None)
