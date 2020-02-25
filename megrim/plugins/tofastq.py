from megrim.environment import MegrimPlugin
from megrim.reference_genome import ReferenceGenome
from megrim.genome_geometry import BamHandler, include_flounder
import pandas as pd
import argparse
import logging
import warnings
import gzip
import bz2
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
import os


FASTQ = "FASTQ"
FQ = "FQ"
FASTA = "FASTA"
FA = "FA"
BAM = "BAM"
SAM = "SAM"
FAST5 = "FAST5"
DIRECTORY = "DIRECTORY"
candidates = {FASTQ: FASTQ,  FQ: FASTQ, FASTA: FASTA, FA: FASTA, BAM: BAM, SAM: SAM, FAST5: FAST5}


def has_compress_suffix(file):
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
    encoding = guess_type(file)[1]
    # print(f"file encoding checked == {encoding}")
    _open = open
    if encoding == "gzip":
        _open = partial(gzip.open, mode="rt")
    if encoding == "bzip2":
        _open = partial(bz2.open, mode="rt")
    with _open(file) as f:
        delim = "+"
        for record in SeqIO.parse(f, 'fastq'):
            print(record.format("fastq"))


def process_file(file, file_type):
    if file_type == FASTQ:
        process_fastq(file)


class BamStats(MegrimPlugin):
    def __init__(self):
        super().__init__()
        self.tool = "ToFastq"
        self.args = None

    def guess_format(self, file=None, iter=False):
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
        for file in os.listdir(self.args.src):
            file_type = self.guess_format(file, True)
            if file_type is not None:
                if (self.args.pickup is None) | (file_type == self.args.pickup.upper()):
                    process_file(os.path.join(self.args.src, file), file_type)

    def execute(self, args):
        self.args = args

        file_format = self.guess_format()
        if file_format == DIRECTORY:
            self.parse_directory()
        else:
            process_file(self.args.src, file_format)



    def arg_params(self, subparsers, parent_parser):
        argparser = subparsers.add_parser(self.tool, help="tofastq help", parents=[parent_parser])
        argparser.add_argument('-s', '--source', metavar="/path/to/input_file", action='store', help='Path to the source sequence file',
                               required=True, dest="src")
        argparser.add_argument('-d', '--directory', action='store_true', help='Process a folder of files rather than a monolith',
                               required=False, dest="dir", default=False)
        argparser.add_argument('-p', '--pickup', action='store', help='For directory parsing, specify filetype to pickup selectively pickup {FASTQ|BAM}',
                               required=False, dest="pickup", default=None)
