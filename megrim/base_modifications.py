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
import pandas as pd
import glob
from megrim.environment import Flounder
from megrim.cigar import cigar_rle
import concurrent.futures
from tqdm import tqdm
import multiprocessing
import sys
import hashlib
import traceback
import numpy as np


def include_flounder(args):
    # setup a Flounder for this workflow ...
    global flounder
    flounder = Flounder()
    flounder.argparse(args)


def fast5_to_basemods(
        fast5file, latest_basecall=None, modification="5mC", threshold=0.75,
        context="CG", force=False):
    """
    Given a FAST5 file, return a DataFrame of filtered base modifications.

    Guppy provides a model for base-modification calling. The results from
    the analysis are written in tabular form to the FAST5 file. This method
    uses the ont_fast5_api to extract the modified bases, to filter on the
    basis of the probability threshold reported and to further reduce using
    the provided sequence context. The result should be a pandas DataFrame
    suitable for other downstream sequence analysis.

    Parameters
    ----------
    fast5file: String
        Path to the FAST5 file that will be processed.
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

    Returns
    -------
    pandas.DataFrame
        The results will be returned in the format of a Pandas DataFrame
        similar to the workflow's previous incarnation in R.

    """
    result = []
    if (not force) & ("flounder" in globals()):
        cached = flounder.read_cache(
            fast5file, pd.DataFrame(), modification, threshold, context)
        if cached is not None:
            return cached

    with get_fast5_file(fast5file, mode="r") as f5:
        for read_id in f5.get_read_ids():
            # print(read_id)
            read = f5.get_read(read_id)
            if latest_basecall is None:
                latest_basecall = read.get_latest_analysis("Basecall_1D")
                # print(latest_basecall)

            mod_base_df = fast5_basemods_to_df(
                read, latest_basecall, modification, threshold, context)
            result.append(mod_base_df)

    result = pd.concat(result, sort=False)
    if "flounder" in globals():
        flounder.write_cache(
            fast5file, result, modification, threshold, context)
    return result


def fast5_basemods_to_df(
        read, latest_basecall, modification, threshold, context):
    mod_base_df = pd.DataFrame(
        data=read.get_analysis_dataset(
            latest_basecall,
            "BaseCalled_template/ModBaseProbs"),
        columns=["A", "6mA", "C", "5mC", "G", "T"])
    mod_base_df['position'] = mod_base_df.index
    mod_base_df['read_id'] = read.read_id

    # reduce the data by modification and probability threshold
    mod_base_df[modification] = mod_base_df[modification] / 255
    mod_base_df = mod_base_df.loc[
        mod_base_df[modification] > threshold,
        ["read_id", "position", modification]]

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




def fast5s_to_basemods(
        path, latest_basecall=None, modification="5mC", threshold=0.75,
        context="CG", force=False, processes=None):
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
    result = []

    if (not force) & ("flounder" in globals()):
        cached = flounder.read_cache(
            path, pd.DataFrame(), modification, threshold, context)
        if cached is not None:
            return cached

    if processes is None:
        processes = multiprocessing.cpu_count()
    files = glob.glob("{}/*.fast5".format(path))

    with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as pool:
        future_fast5 = {
            pool.submit(
                fast5_to_basemods,
                fast5file=file,
                latest_basecall=latest_basecall,
                modification=modification,
                threshold=threshold,
                context=context,
                force=force): file for file in files}
        for future in tqdm(
                concurrent.futures.as_completed(
                    future_fast5), total=len(files)):
            basemods = future.result()
            result.append(basemods)

    result = pd.concat(result, sort=False)
    if "flounder" in globals():
        flounder.write_cache(
            path, result, modification, threshold, context)
    return result


def extract_bam_chunk(bam, chromosome, start, end, force=False):
    """
    Extract minimal bam associated annotation for given coordinates.

    This method will parse bamfile, bam, and return a pandas DataFrame
    containing key observations to facilitate a base-modification workflow.

    Parameters
    ----------
    bam: bam
        An object of class bam.
    chromosome: Str
        The chromosome object of interest.
    start: int
        start position in chromosome.
    end: int
        The chromosomal end position
    force: boolean, optional
        Whether to force the analysis, or whether to allow for a cached result
        to be returned. The default is False.

    Returns
    -------
    data: pd.DataFrame
        Pandas DataFrame containing the parsed entries.

    """
    # recover the cached data if possible
    if (not force) & ("flounder" in globals()):
        data = flounder.read_cache(
            bam.bam, pd.DataFrame(), chromosome, start, end)
    if data is None:
        data = []
        # extract reads from a bam file
        reads = bam.get_sam_core(chromosome, start, end)
        for read in reads:
            # select primary mappings only
            if not any([read.is_secondary, read.is_supplementary,
                        read.is_qcfail, read.is_duplicate]):
                # and select facets as required for basemodification workflow
                # row = pd.Series({"query_name": read.query_name,
                #                  "query_length": read.query_length,
                #                  "reference_name": read.reference_name,
                #                  "reference_start": read.reference_start,
                #                  "strand": "-" if read.is_reverse else "+",
                #                  "cigar": read.cigartuples})
                # creating just a list rather than a pd.Series could be faster???
                row = [read.query_name, read.query_length, read.reference_name,
                       read.reference_start, "-" if read.is_reverse else "+", read.cigartuples]
                data.append(row)

        # convert data into a DataFrame
        data = pd.DataFrame(
            data,
            columns=["query_name", "query_length", "reference_name",
                     "reference_start", "strand", "cigar"])

        # prettify the data
        data.set_index("query_name", drop=False, inplace=True)
        if "flounder" in globals():
            flounder.write_cache(bam.bam, data, chromosome, start, end)
    return data


def mung_bam_variant(rows, chromosome, monly=True):
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
    cigar_rles = {}

    def internal_function(row):
        if row.read_id not in cigar_rles.keys():
            cigar_rles[row.read_id] = cigar_rle(row.cigar)

        try:
            operation = None
            ref_pos = None
            if row.strand == "+":
                # print("processing (+) - pos {} ".format(row.position))
                # operation, ref_pos = cigar_q_to_r(row.position, row.cigar)
                if row.position >= 1:
                    # this is a strange edge case ...
                    operation, ref_pos = cigar_rles[row.read_id].q_to_r(
                        row.position)
                    ref_pos = ref_pos + row.reference_start
                # extract the corresponding nucleotide from the ref genome
                # print line below (commented out) for debugging
                # print(" {} -> {} {}  == {}".format(
                #    row.position, operation, ref_pos, base))
            else:  # maps to the reverse strand
                # print("processing (-) - pos {} ".format(row.position))
                # logic here is that the start of the read is upstream, need
                # to adjust start and end coordinates ...
                pos = (row.query_length - row.position - 1)
                # it is possible for reads that have a 0 position = junk!
                if pos > 0:
                    # operation, ref_pos = cigar_q_to_r(pos, row.cigar)
                    operation, ref_pos = cigar_rles[row.read_id].q_to_r(pos)
                    ref_pos = ref_pos + row.reference_start
            # package the data for return - this is in ancestral R format
            if ref_pos is not None:
                # if the operation is M ...
                if (monly and operation == "M") or not monly:
                    modification = None
                    if "5mC" in row.index:
                        modification = "5mC"
                    elif "6mA" in row.index:
                        modification = "6mA"
                    else:
                        raise ValueError(
                            "Suitable modification column not found")
                    # results = pd.Series({"read_id": row.read_id,
                    #                      "chromosome": chromosome,
                    #                      "pos": ref_pos,
                    #                      "fwd": int(row.strand == "+"),
                    #                      "rev": int(row.strand == "-"),
                    #                      "op": operation,
                    #                      "prob": row[modification],
                    #                      "seq_context": row.contexts})
                    results = [row.read_id, chromosome, ref_pos, int(row.strand == "+"), int(row.strand == "-"), operation, row[modification], row.contexts]
                    return results
        except Exception as ex:
            print(str(ex))
            print(ex)
            print(row)
            print(row.index)
            traceback.print_exc()
            sys.exit(0)

    apply_data = rows.apply(internal_function, axis=1).dropna()

    apply_data = pd.DataFrame(
        [item for item in apply_data],
        columns=["read_id", "chromosome", "pos", "fwd", "rev", "op",
                 "prob", "seq_context"])
    apply_data.set_index("read_id", drop=False, inplace=True)
    return apply_data


def map_methylation_signal_chunk(
        ref, bam, modifications, chromosome, start, end, force=False, processes=None):

    data_hash = hashlib.md5(
        str(
            len(modifications.index)).join(
                modifications.index.unique()).encode()).hexdigest()[0:7]
    idx_hash = hashlib.md5(
        "_".join([chromosome, str(start), str(end)]).encode()).hexdigest()[0:7]
    
    if (not force) & ("flounder" in globals()):
        methylation_chunk = flounder.read_cache(
            bam.bam, pd.DataFrame(), data_hash, idx_hash)

    if methylation_chunk is None:
        if processes is None:
            processes = multiprocessing.cpu_count()
        
        methylation_chunk = []
        # import the bam chunk
        bam_chunk = extract_bam_chunk(bam, chromosome, start, end, force)

        # perform an inner join on the datasets ...
        bam_chunk = bam_chunk.join(modifications, how="inner")
        keys = bam_chunk["read_id"].unique()

        with concurrent.futures.ProcessPoolExecutor(
                max_workers=processes) as pool:
            future_list = []
            # counter = 0
            for key in keys:
                read_chunk = bam_chunk.loc[key]
                if isinstance(read_chunk, pd.Series):
                    read_chunk = pd.DataFrame(read_chunk).T

                # mung_bam_variant(read_chunk, chromosome=chromosome)
                # counter += 1
                # if counter > 10:
                #     sys.exit(0)
                future = pool.submit(
                    mung_bam_variant,
                    read_chunk,
                    chromosome=chromosome)
                future_list.append(future)
            # print("future list has {} members".format(len(future_list)))
            for future in tqdm(concurrent.futures.as_completed(future_list),
                                total=len(future_list)):
                bam_mapped = future.result()
                methylation_chunk.append(bam_mapped)

        methylation_chunk = pd.concat(methylation_chunk)
        if "flounder" in globals():
            flounder.write_cache(
                bam.bam, methylation_chunk, data_hash, idx_hash)
    return methylation_chunk


def map_methylation_signal(ref, bam, modifications, force=False):
    print("Mapping methylation signals")
    
    data_hash = hashlib.md5(
        str(
            len(modifications.index)).join(
                modifications.index.unique()).encode()).hexdigest()
    
    if (not force) & ("flounder" in globals()):
        mapped_reads = flounder.read_cache(
            bam.bam, pd.DataFrame(), data_hash)
    if mapped_reads is None:
        mapped_reads = []
        for chromo in ref.get_chromosomes():
            print(chromo)
            chr_mapped_reads = []
            chromosome_tiles = ref.get_tiled_chromosome(chromo, tile_size=5000000)
            for index in chromosome_tiles.df.index:
                print(index)
                mapped_read_chunk = map_methylation_signal_chunk(
                    ref, bam, modifications,
                    chromosome=chromo,
                    start=chromosome_tiles.df.iloc[index].Start,
                    end=chromosome_tiles.df.iloc[index].End)
                chr_mapped_reads.append(mapped_read_chunk)
            chr_mapped_reads = pd.concat(chr_mapped_reads, sort=False)
            # the reference base calling was dropped earlier due to resources
            # and parallelisation - let's put back the reference bases ...
            chr_mapped_reads = chr_mapped_reads.loc[
                chr_mapped_reads.pos.notna()]
            chr_mapped_reads = chr_mapped_reads.astype(
                dtype={"pos": "int32", "fwd": "int32",
                       "rev": "int32", "prob": "float32"})
            print(chr_mapped_reads)
            
            fasta = list(ref.get_whole_sequence(chromo))
            print("picking bases")
            # bases = itemgetter(chr_mapped_reads.pos.tolist())(fasta)
            bases = [fasta[x] for x in chr_mapped_reads.pos.tolist()]
            print("There are {} bases called - ".format(len(bases)))
            print("setting bases")
            chr_mapped_reads['ref_base'] = bases
            # due to the way that the parse was performed there were duplicates
            # introduced ... remove the duplicates ...
            chr_mapped_reads.drop_duplicates(
                subset=["chromosome", "pos", "prob", "read_id"], inplace=True)
            # and replace the base enumerations
            chr_mapped_reads['A'] = 0
            chr_mapped_reads['C'] = 0
            chr_mapped_reads['G'] = 0
            chr_mapped_reads['T'] = 0
            chr_mapped_reads.loc[chr_mapped_reads.ref_base == "A", "A"] = 1
            chr_mapped_reads.loc[chr_mapped_reads.ref_base == "C", "C"] = 1
            chr_mapped_reads.loc[chr_mapped_reads.ref_base == "G", "G"] = 1
            chr_mapped_reads.loc[chr_mapped_reads.ref_base == "T", "T"] = 1
            mapped_reads.append(chr_mapped_reads)
        mapped_reads = pd.concat(mapped_reads)
        if "flounder" in globals():
            flounder.write_cache(
                bam.bam, mapped_reads, data_hash)
    return mapped_reads


def reduce_mapped_methylation_signal(dataframe, force=False):
    print("reduce_mapped_methylation_signal")
    print("...")
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
                 "ref_base": "first"})
        if "flounder" in globals():
            flounder.write_cache(
                data_hash, df)
    return df


if __name__ == '__main__':

    pd.set_option("display.max_columns", None)
    pd.set_option("display.expand_frame_repr", False)
    pd.set_option("max_colwidth", -1)

    test = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/A1-D1-PAD851010.fast5"
    with get_fast5_file(test, mode="r") as f5:

        read_id = "ff873e4e-af54-4f52-a5f2-c5ea2e1800a3"
        read = f5.get_read(read_id)
        latest_basecall = read.get_latest_analysis("Basecall_1D")

        mod_base_df = fast5_basemods_to_df(
            read, latest_basecall, "5mC", 0.75, "CG")
        print(mod_base_df)


    sys.exit(0)

    flounder = Flounder()
    f5path = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/"

    # define working parameters
    modification = "5mC"
    threshold = 0.75
    context = "CG"

    # parse modifications from available fast5 files
    modifications = fast5s_to_basemods(f5path, modification=modification,
                                       threshold=threshold, context=context)
    modifications.set_index("read_id", drop=False, inplace=True)
    # print out the modifications - quick reality checkk
    print(modifications)

    # define a reference genome and bam file
    reference = "/Volumes/Samsung_T5/MethylationPyTutorial/ReferenceData/human_g1k_v37.chr20.fasta"
    reference = ReferenceGenome(reference)
    bam = "/Volumes/Samsung_T5/MethylationPyTutorial/Analysis/minimap2/Native.bam"
    bam = BamHandler(bam)



    # associated mapped bases with the available modifications
    mapped_reads = map_methylation_signal(reference, bam, modifications)
    print(mapped_reads)
    reduced_reads = reduce_mapped_methylation_signal(mapped_reads)
    print(reduced_reads)
    
    
    
    