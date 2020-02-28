"""
Module to enable base-modification analyses for tutorial framework.

The intention of this module is to provide functionality to parse FAST5
files to replicate the earlier R-based functionality aiming to demystify
the processes of base-modification calling. This is based on the ont-fast5-api
and documentation that was posted

https://community.nanoporetech.com/posts/guppy-release-v3-2

Created on Mon Feb 10 15:52:14 2020

@author: srudd
"""


from ont_fast5_api.fast5_interface import get_fast5_file
import pandas as pd
import glob
from megrim.environment import Flounder
from megrim.cigar import cigar_rle
import concurrent.futures
from tqdm import tqdm
import hashlib
import argparse
import numpy as np
import logging
import pyranges as pr
from Bio.Seq import Seq


class BaseModifications(Flounder):
    """A class to perform base modification extraction from FAST5 context."""

    def __init__(self, fast5, bam, reference, modification="5mC", 
                 threshold=0.85, context = "CG", args=None):
        Flounder.__init__(self)
        if args is not None:
            if isinstance(args, argparse.Namespace):
                self.argparse(args)
            if isinstance(args, dict):
                self.dictparse(args)
            include_flounder(args)
            # functions in module but outside of class may require these info

        self.fast5 = fast5
        self.reference = reference
        self.bam = bam
        self.modification = modification
        self.threshold = threshold
        self.context = context
        self.index = hashlib.md5(
                    f" {self.fast5} {self.reference.fasta.filename} "
                    "{self.bam.bam} {self.modification} {self.context} "
                    "{self.threshold} ".encode()).hexdigest()[0:7]
        logging.info(f"workflow_index == [{self.index}]")

    def fast5s_to_basemods(self, processes=None, force=False):
        """
        Collate primitive base-modification information from defined folder.

        This method implements the fast5_to_basemods function across a whole
        directory of fast5 files - the results will, if possible, be cached
        within the Flounder framework for faster recall and recovery.

        Parameters
        ----------
        force: boolean
            If True any existing cached (Flounder) data will be ignored.
        processes: int
            The number of threads to use in data calculation. The default is
            None and the number of available computer cores will be used.

        Returns
        -------
        pandas.DataFrame
            The results will be returned in the format of a Pandas DataFrame
            similar to the workflow's previous incarnation in R.

        """
        result = None
        fast5_hash = hashlib.md5(self.fast5.encode()).hexdigest()[0:7]

        if not force:
            result = self.read_cache(
                fast5_hash, pd.DataFrame(), self.modification, self.context)

        if result is None:
            result = []
            if processes is None:
                processes = self.thread_processes
            files = glob.glob("{}/*.fast5".format(self.fast5))

            # {result.append(
            #      self.fast5_to_basemods(file, threshold=0)): \
            #      file for file in files}
            with concurrent.futures.ProcessPoolExecutor(
                    max_workers=processes) as pool:
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
                    fast5_hash, result, self.modification, self.context)
        return result

    def filter_modifications_by_prob(self, force=False):
        """
        Filter the identified base-modifications by the probability threshold.

        The class initialisation defines a threshold at which identified
        base modifications can be stratified to those of interest, vs those
        that can be ignored. This method takes a full dataset of identified
        modifications and filters by this threshold. While there are
        performance advantages in filtering whilst parsing the FAST5 files,
        by caching the complete dataset (once) this provides a massive boost
        in relative performance during e.g. parameter sweeps.

        Parameters
        ----------
        force: boolean, optional
            Whether the prior steps in the workflow should be recalculated or
            if the cached data can be used. The default is False and the data
            will not be recalculated.

        Returns
        -------
        f_modifications: pd.DataFrame
            Dataframe of per fast5 sequence probability filtered basemods.

        """
        f_modifications = None
        if not force:
            f_modifications = self.read_cache(self.index, pd.DataFrame())

        if f_modifications is None:
            modifications = self.fast5s_to_basemods()
            f_modifications = modifications.loc[modifications[
                self.modification] >= float(self.threshold)]
            self.write_cache(self.index, f_modifications)
        return f_modifications

    def map_methylation_signal(self, force=False):
        """
        Map the per fast5-sequence base modifications to the provided BAM.

        The key step in the workflow; this is responsible for taking the
        raw per-read coordinates for base-modifications and mapping onto
        the reference genome thus yielding genomic positional information.

        Parameters
        ----------
        force: boolean, optional
            Whether the prior steps in the workflow should be recalculated or
            if the cached data can be used. The default is False and the data
            will not be recalculated.

        Returns
        -------
        chr_mapped_reads : pd.DataFrame
            Dataframe of per read genomic mapping coordinates.

        """
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

                mapped_read_chunk = self.map_methylation_signal_chunk(
                    bam_chunk, modifications)

                logging.debug("cleaning up mapped_methylation dataset")
                mapped_read_chunk = mapped_read_chunk.loc[
                    mapped_read_chunk.pos.notna()]
                mapped_read_chunk = mapped_read_chunk.astype(
                    dtype={"pos": "int32", "fwd": "int32",
                           "rev": "int32", "prob": "float32"})

                # the mapping results may contain basemods that are derived
                # from mapping events other than match - nuke them (initially)
                mapped_read_chunk = mapped_read_chunk.loc[
                    mapped_read_chunk.op == "M"]

                # since the reads that span the target region extend beyond
                # the boundaries of the target, we should clip the read set to
                # include basemod loci that fall within canonical boundaries
                logging.debug("removing off target content")
                mapped_read_chunk = mapped_read_chunk.loc[
                    (mapped_read_chunk.pos >= block_start) &
                    (mapped_read_chunk.pos < block_end)]

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

                mapped_read_chunk['fwd_cov'] = \
                    rle[chromosome, "+"][mapped_read_chunk.pos]
                mapped_read_chunk['rev_cov'] = \
                    rle[chromosome, "-"][mapped_read_chunk.pos]
                logging.debug(mapped_read_chunk)

                # should we also consider reference context since there are
                # a lot of base assignments where read context != reference
                # context, even with depth-of-coverage?
                def extract_reference_context(
                        position, reverse, offset=0, rev_comp=True):
                    start = int(position)
                    end = int(position)+len(self.context)
                    if reverse == 1:
                        # adjust coordinates ...
                        end = int(position) + 1
                        start = int(position) - len(self.context) + 1
                    if offset > 0:
                        ref_context = "{}[{}]{}".format(
                            "".join(fasta[max(start-offset, 0):start]),
                            "".join(fasta[start:end]),
                            "".join(fasta[end:min(end+offset, len(fasta))]))
                    else:
                        ref_context = "".join(fasta[start:end])
                    if (reverse == 1) & rev_comp:
                        # reverse complement ...
                        ref_context = str(
                            Seq(ref_context).reverse_complement())
                    return ref_context

                # ref_context = ["".join(fasta[x-5:x+len(self.context)+4])
                #               for x in mapped_read_chunk.pos.tolist()]

                ref_context = list(
                    map(extract_reference_context, mapped_read_chunk.pos,
                        mapped_read_chunk.rev))

                mapped_read_chunk['ref_context'] = ref_context
                # strip out any of the reference contexts that are not the
                # same as the target context - this means that we'll reject
                # any sequence loci where subject != reference
                mapped_read_chunk = mapped_read_chunk.loc[
                    mapped_read_chunk['ref_context'] == self.context]

                # and expand the reference context string just a little
                mapped_read_chunk['ref_context'] = list(
                    map(extract_reference_context, mapped_read_chunk.pos,
                        mapped_read_chunk.rev, 5, False))

                chr_mapped_reads.append(mapped_read_chunk)
                logging.debug(mapped_read_chunk)

            chr_mapped_reads = pd.concat(chr_mapped_reads, sort=False)
            logging.debug(chr_mapped_reads)

            self.write_cache(self.index, chr_mapped_reads)
        return chr_mapped_reads

    def map_methylation_signal_chunk(
            self, bam_chunk, modifications, force=False, processes=None):
        """
        Assign genomic coordinates to block of mapped basemod sequence context.

        Parameters
        ----------
        bam_chunk : TYPE
            DESCRIPTION.
        modifications : TYPE
            DESCRIPTION.
        force : TYPE, optional
            DESCRIPTION. The default is False.
        processes : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        methylation_chunk : TYPE
            DESCRIPTION.

        """
        chromosome = bam_chunk.iloc[0].reference_name
        block_start = bam_chunk.iloc[0].block_start
        block_end = bam_chunk.iloc[0].block_end
        bam_hash = hashlib.md5(
            f"{chromosome} {block_start} {block_end}".encode()).\
            hexdigest()[0:7]

        if not force:
            methylation_chunk = self.read_cache(
                self.index, pd.DataFrame(), bam_hash)

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
                for future in tqdm(
                        concurrent.futures.as_completed(future_list),
                        total=len(future_list)):
                    bam_mapped = future.result()
                    for row in bam_mapped:
                        methylation_chunk.append(row)

            methylation_chunk = pd.DataFrame(
                methylation_chunk,
                columns=["read_id", "chromosome", "pos", "op", "prob", "fwd",
                         "rev", "seq_context"])
            methylation_chunk.set_index("read_id", drop=False, inplace=True)
            self.write_cache(self.index, methylation_chunk, bam_hash)
        return methylation_chunk

    def reduce_mapped_methylation_signal(self, force=False):
        """
        Reduce dimensions of basemod data from per-read to per-genomic-locus.

        Parameters
        ----------
        force : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        df : TYPE
            DESCRIPTION.

        """
        logging.info("reduce_mapped_methylation_signal")
        df = None
        if not force:
            df = self.read_cache(self.index, pd.DataFrame())
        if df is None:
            dataframe = self.map_methylation_signal(force=force)

            df = dataframe.groupby(["chromosome", "pos"]).agg(
                {"chromosome": "first", "pos": "first", "prob": np.mean,
                 "fwd": np.sum, "rev": np.sum, "seq_context": "first",
                 "ref_base": "first", "fwd_cov": "first", "rev_cov": "first"})
            self.write_cache(self.index, df)
        return df


def include_flounder(args):
    """
    Set a flounder context visible in module global space for data caching.

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # setup a Flounder for this workflow ...
    global flounder
    flounder = Flounder()
    if isinstance(args, argparse.Namespace):
        flounder.argparse(args)
    if isinstance(args, dict):
        flounder.dictparse(args)


def lfast5_basemods_to_df(read, latest_basecall, modification, context):
    """
    Perform local extraction of basemods from fast5 sequence entries.

    Parameters
    ----------
    read: fast5_sequence_read
        This is an ont_fast5_api sequence read.
    latest_basecall: str
        The basecall version that was determined in parental step.
    modification: str
        The base-modification to extract.
    context: str
        The local sequence context to filter results by.

    Returns
    -------
    list-of-lists
        Lists of per-read base-modifications and probabilities.

    """
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
            return "{}[{}]{}".format(
                sequence[max(int(pos)-offset, 0): int(pos)],
                sequence[int(pos): int(pos2)],
                sequence[int(pos2): min(int(pos2)+offset, len(sequence)-1)])
        return sequence[max(int(pos), 0):min(int(pos2), len(sequence)-1)]

    offset = 0
    mod_base_df['contexts'] = list(map(
        a_context, mod_base_df.position, mod_base_df.position+len(context)))
    mod_base_df = mod_base_df.loc[mod_base_df.contexts == context]
    offset = 5
    mod_base_df['contexts'] = list(map(
        a_context, mod_base_df.position, mod_base_df.position + len(context)))
    return mod_base_df


def lfast5_to_basemods(fast5file, modification, context, force=False):
    """
    Extract base-modification data at per-read resolution from FAST5 file.

    Parameters
    ----------
    fast5file : TYPE
        DESCRIPTION.
    modification : TYPE
        DESCRIPTION.
    context : TYPE
        DESCRIPTION.
    force : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    result : TYPE
        DESCRIPTION.

    """
    result = None
    if (not force) & ("flounder" in globals()):
        result = flounder.read_cache(
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

                mod_base_df = lfast5_basemods_to_df(
                    read, latest_basecall, modification, context)
                result.append(mod_base_df)

        result = pd.concat(result, sort=False)
        if "flounder" in globals():
            flounder.write_cache(
                fast5file, result, modification, context)
    return result


def mung_bam_variant(rows):
    """
    Extract genomic coordinates from a pd.DataFrame block of mapping info.

    Parameters
    ----------
    rows : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
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

    def anchor_base(
            position, strand, reference_start, query_length, contexts, prob):
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
        return rows.iloc[0].read_id, chromosome, ref_pos, operation, prob, \
            int(strand == "+"), int(strand == "-"), contexts

    return list(map(anchor_base, rows.position, rows.strand,
                    rows.reference_start, rows.query_length, rows.contexts,
                    rows[modification]))
