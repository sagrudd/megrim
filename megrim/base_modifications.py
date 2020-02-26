#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The intention of this module is to provide functionality to parse FAST5
files to replicate the earlier R-based functionality aiming to demystify
the processes of base-modification calling. This is based on the ont-fast5-api
and documentation that was posted

https://community.nanoporetech.com/posts/guppy-release-v3-2

Created on Mon Feb 10 15:52:14 2020

@author: srudd
"""


from ont_fast5_api.fast5_interface import get_fast5_file
from megrim.genome_geometry import BamHandler
from megrim.reference_genome import ReferenceGenome
import pandas as pd
import glob
from megrim.environment import Flounder
from megrim.cigar import cigar_rle
import concurrent.futures
from tqdm import tqdm
import multiprocessing
import sys
import hashlib
import argparse
import numpy as np
import logging
import pyranges as pr



class BaseModifications(Flounder):

    def __init__(self, fast5, bam, reference, args=None):
        Flounder.__init__(self)
        if args is not None:
            if isinstance(args, argparse.Namespace):
                self.argparse(args)
            if isinstance(args, dict):
                self.dictparse(args)

        self.fast5 = fast5
        self.reference = reference
        self.bam = bam
        self.modification = "5mC"
        self.threshold = 0.85
        self.context = "CG"
        self.index = hashlib.md5(
                    f" {self.fast5} {self.reference.fasta.filename} {self.bam.bam} {self.modification} {self.context} {self.threshold} ".encode()).hexdigest()[0:7]
        logging.info(f"workflow_index == [{self.index}]")



    def fast5s_to_basemods(self, processes=None, force=False):
        """
        Collate primitive base-modification information from defined folder.

        This method implements the fast5_to_basemods function across a whole
        directory of fast5 files - the results will, if possible, be cached within
        the Flounder framework for faster recall and recovery.

        Parameters
        ----------
        path: String
            Path to the FAST5 directory that will be processed.
        latest_basecall: String
            The ont_fast5_api parses the FAST5 file for a version of a basecall
            analysis suitable for pulling out the basemod information. This can
            be specified at the command_line if required.
        modification: String
            The base-modification to be reported. The default is "5mC". The other
            allowed value in "6mA"
        threshold: float
            The probability threshold with which to filter the returned base
            modifications. The default is 0.75.
        context: String
            The local sequence context to use. The default is "CG". If the
            base-modification does not appear within this context it will not be
            returned - there is some TODO here to expand utility.
        force: boolean
            If True any existing cached (Flounder) data will be ignored.
        processes: int
            The number of threads to use in data calculation. The default is None
            and the number of available computer cores will be used.

        Returns
        -------
        pandas.DataFrame
            The results will be returned in the format of a Pandas DataFrame
            similar to the workflow's previous incarnation in R.

        """
        result = None

        if not force:
            result = self.read_cache(
                hashlib.md5(self.fast5.encode()).hexdigest()[0:7], pd.DataFrame(), self.modification, self.context)

        if result is None:
            result = []
            if processes is None:
                processes = self.thread_processes
            files = glob.glob("{}/*.fast5".format(self.fast5))

            # {result.append(self.fast5_to_basemods(file, threshold=0)): file for file in files}
            with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as pool:
                future_fast5 = {
                    pool.submit(
                        lfast5_to_basemods,
                        file,
                        context=self.context,
                        modification=self.modification,
                        force=force): file for file in files}
                for future in tqdm(
                        concurrent.futures.as_completed(
                            future_fast5), total=len(files)):
                    result.append(future.result())

            result = pd.concat(result, sort=False)
            result.set_index("read_id", drop=False, inplace=True)
            self.write_cache(
                    hashlib.md5(self.fast5.encode()).hexdigest()[0:7], result, self.modification, self.context)
        return result


    def filter_modifications_by_prob(self, force=False):
        f_modifications = None
        if not force:
            f_modifications = self.read_cache(self.index, pd.DataFrame())

        if f_modifications is None:
            # derive results from within class - saves a long chain of maintaining variables ...
            modifications = self.fast5s_to_basemods()

            f_modifications = modifications.loc[modifications[self.modification] >= float(self.threshold)]
            if "flounder" in globals():
                flounder.write_cache(self.index, f_modifications)
        return f_modifications

    def map_methylation_signal(self, force=False):
        logging.info("Mapping methylation signals")

        if not force:
            chr_mapped_reads = self.read_cache(self.index, pd.DataFrame())

        if chr_mapped_reads is None:
            chr_mapped_reads = []

            modifications = self.filter_modifications_by_prob(force=force)

            for bam_chunk in self.bam.chunk_generator():
                chromosome = bam_chunk.iloc[0].reference_name
                block_start = bam_chunk.iloc[0].block_start
                block_end = bam_chunk.iloc[0].block_end
                fasta = list(self.reference.get_whole_sequence(chromosome))

                mapped_read_chunk = self.map_methylation_signal_chunk(bam_chunk, modifications)

                logging.debug("cleaning up mapped_methylation dataset")
                mapped_read_chunk = mapped_read_chunk.loc[
                    mapped_read_chunk.pos.notna()]
                mapped_read_chunk = mapped_read_chunk.astype(
                    dtype={"pos": "int32", "fwd": "int32",
                           "rev": "int32", "prob": "float32"})

                # since the reads that span the target region extend beyond the boundaries of
                # the target, we should clip the read set to include basemod loci that fall
                # within the canonical boundaries
                logging.debug("removing off target content")
                mapped_read_chunk = mapped_read_chunk.loc[
                    (mapped_read_chunk.pos >= block_start) &
                    (mapped_read_chunk.pos < block_end)]

                # the reference base calling was dropped earlier due to resources
                # and parallelisation - let's put back the reference bases ...
                logging.debug("adding reference sequence base")
                bases = [fasta[x] for x in mapped_read_chunk.pos.tolist()]
                mapped_read_chunk['ref_base'] = bases

                # and augment these information with some additional context
                # relating to depth of coverage and strandedness
                logging.debug("calculating pyranges")
                rle = pr.PyRanges(chromosomes=bam_chunk.reference_name,
                                  starts=bam_chunk.reference_start,
                                  ends=bam_chunk.reference_end,
                                  strands=bam_chunk.strand).to_rle(strand=True)

                mapped_read_chunk['fwd_cov'] = rle[chromosome, "+"][mapped_read_chunk.pos]
                mapped_read_chunk['rev_cov'] = rle[chromosome, "-"][mapped_read_chunk.pos]
                logging.debug(mapped_read_chunk)

                chr_mapped_reads.append(mapped_read_chunk)

            chr_mapped_reads = pd.concat(chr_mapped_reads, sort=False)
            logging.debug(chr_mapped_reads)

            self.write_cache(self.index, chr_mapped_reads)
        return chr_mapped_reads

    def map_methylation_signal_chunk(self, bam_chunk, modifications, force=False, processes=None):

        chromosome = bam_chunk.iloc[0].reference_name
        block_start = bam_chunk.iloc[0].block_start
        block_end = bam_chunk.iloc[0].block_end
        bam_hash = hashlib.md5(f"{chromosome} {block_start} {block_end}".encode()).hexdigest()[0:7]

        if not force:
            methylation_chunk = self.read_cache(self.index, pd.DataFrame(), bam_hash)

        if methylation_chunk is None:
            if processes is None:
                processes = self.thread_processes

            methylation_chunk = []
            # perform an inner join on the datasets ...
            bam_chunk = bam_chunk.join(modifications, how="inner")
            keys = bam_chunk["read_id"].unique()

            with concurrent.futures.ProcessPoolExecutor(
                    max_workers=processes) as pool:
                future_list = []
                for key in tqdm(keys):
                    read_chunk = bam_chunk.loc[key, :]
                    # mung_bam_variant(read_chunk)
                    future = pool.submit(
                        mung_bam_variant,
                        read_chunk)
                    future_list.append(future)
                # print("future list has {} members".format(len(future_list)))
                for future in tqdm(concurrent.futures.as_completed(future_list), total=len(future_list)):
                    bam_mapped = future.result()
                    for row in bam_mapped:
                        methylation_chunk.append(row)

            methylation_chunk = pd.DataFrame(
                methylation_chunk,
                columns=["read_id", "chromosome", "pos", "op", "prob", "fwd", "rev", "seq_context"])
            methylation_chunk.set_index("read_id", drop=False, inplace=True)
            self.write_cache(self.index, methylation_chunk, bam_hash)
        return methylation_chunk


def include_flounder(args):
    # setup a Flounder for this workflow ...
    global flounder
    flounder = Flounder()
    if isinstance(args, argparse.Namespace):
        flounder.argparse(args)
    if isinstance(args, dict):
        flounder.dictparse(args)



def lfast5_basemods_to_df(read, latest_basecall, modification, context):
    mod_base_df = pd.DataFrame(
        data=read.get_analysis_dataset(
            latest_basecall,
            "BaseCalled_template/ModBaseProbs"),
        columns=["A", "6mA", "C", "5mC", "G", "T"])
    mod_base_df['position'] = mod_base_df.index
    mod_base_df['read_id'] = read.read_id

    # reduce the data by modification and probability threshold
    mod_base_df[modification] = mod_base_df[modification] / 255
    mod_base_df = mod_base_df.loc[:, ["read_id", "position", modification]]

    # reduce by local sequence context
    fastq = read.get_analysis_dataset(
        latest_basecall, "BaseCalled_template/Fastq")
    sequence = fastq.splitlines()[1]

    # instead of just a lambda function, this is a little over-
    # engineered; expecting more interesting IUPAC based contexts in
    # future
    def a_context(pos, pos2):
        if offset > 0:
            return "{}[{}]{}".format(sequence[max(int(pos)-offset, 0): int(pos)],
                                     sequence[int(pos): int(pos2)],
                                     sequence[int(pos2): min(int(pos2)+offset, len(sequence)-1)])
        return sequence[max(int(pos), 0):min(int(pos2), len(sequence)-1)]

    offset = 0
    mod_base_df['contexts'] = list(map(a_context, mod_base_df.position, mod_base_df.position+len(context)))
    mod_base_df = mod_base_df.loc[mod_base_df.contexts == context]
    offset = 5
    mod_base_df['contexts'] = list(map(a_context, mod_base_df.position, mod_base_df.position + len(context), ))
    return mod_base_df





def lfast5_to_basemods(fast5file, modification, context, force=False):
    print(fast5file)
    result = None
    if (not force) & ("flounder" in globals()):
        cached = flounder.read_cache(
            fast5file, pd.DataFrame(), modification, context)

    if result is None:
        latest_basecall = None
        result = []
        with get_fast5_file(fast5file, mode="r") as f5:
            for read_id in f5.get_read_ids():
                # print(read_id)
                read = f5.get_read(read_id)
                if latest_basecall is None:
                    latest_basecall = read.get_latest_analysis("Basecall_1D")

                mod_base_df = lfast5_basemods_to_df(read, latest_basecall, modification, context)
                result.append(mod_base_df)

        result = pd.concat(result, sort=False)
        if "flounder" in globals():
            flounder.write_cache(
                fast5file, result, modification, context)
    return result




















def mung_bam_variant(rows, monly=True):
    """
    Sync a base-modification variant with FAST5 base modification data.

    Parameters
    ----------
    row : TYPE
        DESCRIPTION.
    ref : TYPE
        DESCRIPTION.
    chromosome : TYPE
        DESCRIPTION.
    modification : TYPE
        DESCRIPTION.

    Returns
    -------
    results : TYPE
        DESCRIPTION.

    """
    if isinstance(rows, pd.Series):
        rows = pd.DataFrame(rows).T
    # fix the modification column name
    modification = None
    if "5mC" in rows.columns:
        modification = "5mC"
    elif "6mA" in rows.columns:
        modification = "6mA"
    else:
        raise ValueError(
            "Suitable modification column not found")

    rle = cigar_rle(rows.iloc[0].cigar)
    chromosome = rows.iloc[0].reference_name

    def anchor_base(position, strand, reference_start, query_length, contexts, prob):
        ref_pos = -1
        operation = "?"
        if strand == "+":
            if position >= 1:
                operation, ref_pos = rle.q_to_r(position)
        else:  # maps to the reverse strand
            pos = (query_length - position - 1)
            if pos > 0:
                operation, ref_pos = rle.q_to_r(pos)
        ref_pos = ref_pos + reference_start
        return rows.iloc[0].read_id, chromosome, ref_pos, operation, prob, int(strand == "+"), int(strand == "-"), contexts

    return list(map(anchor_base, rows.position, rows.strand, rows.reference_start, rows.query_length, rows.contexts, rows[modification]))





def reduce_mapped_methylation_signal(dataframe, force=False):
    logging.info("reduce_mapped_methylation_signal")
    data_hash = hashlib.md5(
        str(
            len(dataframe.index)).join(
                dataframe.index.unique()).encode()).hexdigest()
    df = None
    if (not force) & ("flounder" in globals()):
        df = flounder.read_cache(
            data_hash, pd.DataFrame())
    if df is None:
        df = dataframe.groupby(["chromosome", "pos"]).agg(
                {"chromosome": "first", "pos": "first", "prob": np.mean,
                 "fwd": np.sum, "rev": np.sum, "seq_context": "first",
                 "ref_base": "first", "fwd_cov": "first", "rev_cov": "first"})
        if "flounder" in globals():
            flounder.write_cache(
                data_hash, df)
    return df


