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
from megrim.reference_genome import ReferenceGenome
from megrim.genome_geometry import BamHandler
from megrim.cigar import cigar_q_to_r
import concurrent.futures
from tqdm import tqdm
import multiprocessing
import sys
from Bio.Seq import Seq
import hashlib


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
    result = pd.DataFrame()
    if (not force) & ("flounder" in globals()):
        cached = flounder.read_cache(
            fast5file, result, modification, threshold, context)
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
            result = result.append(mod_base_df, sort=False)

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
    def is_context(row):
        position = row.position
        localcontext = sequence[position:position+len(context)]
        if localcontext == context:
            return True
        return False

    def get_context(row, width=5):
        position = row.position
        start = position-width
        if start < 0:
            start = 0
        end = position+len(context)+width
        if end >= len(fastq):
            end = len(fastq)-1
        return sequence[start: end]

    if len(mod_base_df.index) > 0:
        contextful = mod_base_df.apply(is_context, axis=1)
        mod_base_df = mod_base_df[contextful]

    # finally augment the tabular data with actual sequence context
    if len(mod_base_df.index) > 0:
        mod_base_df["contexts"] = mod_base_df.apply(
            get_context, axis=1)
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
    result = pd.DataFrame()

    if (not force) & ("flounder" in globals()):
        cached = flounder.read_cache(
            path, result, modification, threshold, context)
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
            result = result.append(basemods, sort=False)

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
        data = pd.DataFrame()
        # extract reads from a bam file
        reads = bam.get_sam_core(chromosome, start, end)
        for read in reads:
            # select primary mappings only
            if not any([read.is_secondary, read.is_supplementary,
                        read.is_qcfail, read.is_duplicate]):
                # and select facets as required for basemodification workflow
                row = pd.Series({"query_name": read.query_name,
                                 "query_length": read.query_length,
                                 "reference_name": read.reference_name,
                                 "reference_start": read.reference_start,
                                 "strand": "-" if read.is_reverse else "+",
                                 "cigar": read.cigartuples})
                data = data.append(row, ignore_index=True)
        # prettify the data
        data.set_index("query_name", drop=False, inplace=True)
        if "flounder" in globals():
            flounder.write_cache(bam.bam, data, chromosome, start, end)
    return data


def sync_bam_variant(row, ref, chromosome, modification):
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
    # print(row)
    base = None
    try:
        if row.strand == "+":
            # print("processing (+) - pos {} ".format(row.position))
            operation, ref_pos = cigar_q_to_r(row.position, row.cigar)
            ref_pos = ref_pos + row.reference_start
            # extract the corresponding nucleotide from the ref genome
            base = str(Seq(ref.get_sequence(
                chromosome, ref_pos, ref_pos)).complement())
            # print line below (commented out) for debugging
            # print(" {} -> {} {}  == {}".format(
            #    row.position, operation, ref_pos, base))
        else:  # maps to the reverse strand
            # print("processing (-) - pos {} ".format(row.position))
            # logic here is that the start of the read is upstream, need
            # to adjust start and end coordinates ...
            pos = (row.query_length - row.position - 2)
            # it is possible for a few reads that we have a 0 position - junk!
            if pos > 0:
                operation, ref_pos = cigar_q_to_r(pos, row.cigar)
                ref_pos = ref_pos + row.reference_start
                base = ref.get_sequence(chromosome, ref_pos, ref_pos)

        # package the data for return - this is in ancestral R format
        if base is not None:
            results = pd.Series({"read_id": row.read_id,
                                 "pos": ref_pos,
                                 "fwd": int(row.strand == "+"),
                                 "rev": int(row.strand == "-"),
                                 "op": operation,
                                 "A": int(base == "A"),
                                 "C": int(base == "C"),
                                 "G": int(base == "G"),
                                 "T": int(base == "T"),
                                 "prob": row[modification],
                                 "seq_context": row.contexts})
            return results
    except Exception as ex:
        print(str(ex))
        print(row)


def layer_modifications_on_bam(
        file, ref, data, chromosome, start, latest_basecall, modification,
        threshold, context, force):
    bam_mapped = None
    if (not force) & ("flounder" in globals()):
        bam_mapped = flounder.read_cache(
            bam.bam, pd.DataFrame(), chromosome, start,
            hashlib.md5(file.encode()).hexdigest()[0:7])
    if bam_mapped is None:
        basemods = fast5_to_basemods(
            file, latest_basecall=latest_basecall, modification=modification,
            threshold=threshold, context=context, force=force)
        basemods.set_index("read_id", drop=False, inplace=True)

        # identify the overlap in the data ...
        join = data.join(basemods, how="inner")

        # and mung these data into the previous R format
        bam_mapped = pd.DataFrame(
            join.apply(
                sync_bam_variant, axis=1, ref=ref, chromosome=chromosome,
                modification=modification))
        if "flounder" in globals():
            flounder.write_cache(
                bam.bam, pd.DataFrame(), chromosome, start,
                hashlib.md5(file.encode()).hexdigest()[0:7])
    return bam_mapped


def old_map_methylation_chunk(
        ref, bam, f5path, chromosome, start, end, latest_basecall=None,
        modification="5mC", threshold=0.75, context="CG", force=False):
    methylation_chunk = None

    # let's try to cache the result sets ..
    if (not force) & ("flounder" in globals()):
        methylation_chunk = flounder.read_cache(
            bam.bam, pd.DataFrame(), "methyl_chunk", chromosome, start, end)
    if methylation_chunk is None:
        methylation_chunk = pd.DataFrame()
        # import the bam chunk
        data = extract_bam_chunk(bam, chromosome, start, end, force)
        # layer on the base-modification content
        files = glob.glob("{}/*.fast5".format(f5path))
        for file in tqdm(files):
            bam_mapped = layer_modifications_on_bam(
                file, ref, data, chromosome, start, latest_basecall,
                modification, threshold, context, force)
            methylation_chunk = methylation_chunk.append(
                bam_mapped, sort=False)
        # and cache the data ...
        if "flounder" in globals():
            flounder.write_cache(bam.bam, pd.DataFrame(), "methyl_chunk",
                                 chromosome, start, end)
    return methylation_chunk




def fast5_indexer(fast5file):
    fast5_index = pd.DataFrame()
    with get_fast5_file(fast5file, mode="r") as f5:
        for read_id in f5.get_read_ids():
            row = pd.Series({"fast5file": fast5file, "read_id": read_id})
            fast5_index = fast5_index.append(row, ignore_index=True)
    return fast5_index


def index_fast5_content(f5path, force=False, processes=8):
    files = glob.glob("{}/*.fast5".format(f5path))
    fast5_index = None
    hashstr = hashlib.md5("".join(files).encode()).hexdigest()
    if (not force) & ("flounder" in globals()):
        fast5_index = flounder.read_cache(
            f5path, pd.DataFrame(), hashstr)
    if fast5_index is None:
        print("indexing fast5 content")
        fast5_index = pd.DataFrame()

        with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as pool:
            future_fast5 = {
                pool.submit(fast5_indexer, fast5file=file):
                    file for file in files
                    }
            for future in tqdm(
                concurrent.futures.as_completed(future_fast5),
                total=len(files)):

                index = future.result()
                fast5_index = fast5_index.append(index, ignore_index=True)

        if "flounder" in globals():
            flounder.write_cache(f5path, fast5_index, hashstr)
    return fast5_index



def extract_f5_basemods(multiindex, bam_chunk, ref_filename, chromosome, modification, threshold, context):
    ref = ReferenceGenome(ref_filename)
    multiresult = None

    for index in multiindex:
        row = bam_chunk.loc[index]
        # print("processing {}".format(index))

        with get_fast5_file(row.fast5file, mode="r") as f5:
            read = f5.get_read(row.name)
            latest_basecall = read.get_latest_analysis("Basecall_1D")

            mod_base_df = fast5_basemods_to_df(
                read, latest_basecall, modification, threshold,
                context)

            mod_base_df.set_index("read_id", drop=False, inplace=True)
            mod_base_df['strand'] = row.strand
            mod_base_df['query_length'] = row.query_length
            mod_base_df['cigar'] = row.cigar
            mod_base_df['reference_start'] = row.reference_start

            # print(mod_base_df.drop(["cigar"], axis=1))
            # identify the overlap in the data ...

            # and mung these data into the previous R format
            bam_mapped = pd.DataFrame(
                mod_base_df.apply(
                    sync_bam_variant, axis=1, ref=ref,
                    chromosome=chromosome,
                    modification=modification))

            if not bam_mapped.empty:
                # print(bam_mapped)
                if multiresult is None:
                    multiresult = bam_mapped
                else:
                    multiresult = pd.concat([multiresult, bam_mapped], ignore_index=True)
                # multiresult.append(bam_mapped, sort=False)
            
    # print(multiresult)
    return multiresult


def map_methylation_chunk(
        ref, bam, f5index, chromosome, start, end, latest_basecall=None,
        modification="5mC", threshold=0.75, context="CG", force=False):
    methylation_chunk = None

    # let's try to cache the result sets ..
    if (not force) & ("flounder" in globals()):
        methylation_chunk = flounder.read_cache(
            bam.bam, pd.DataFrame(), "methyl_chunk", chromosome, start, end)
    if methylation_chunk is None:
        methylation_chunk = pd.DataFrame()
        # import the bam chunk
        bam_chunk = extract_bam_chunk(bam, chromosome, start, end, force)

        bam_chunk["fast5file"] = f5index.loc[bam_chunk.index]["fast5file"]

        # ProcessPoolExecutor appears to have limits for # of tasks ...
        items = len(bam_chunk.index)
        chunk_list = bam_chunk.index.tolist()
        chunksize = 25
        if items > 10000:
            chunksize = 250

        def chunk(lst, n):
            n = max(1, n)
            return (lst[i:i + n] for i in range(0, len(lst), n))

        chunks = chunk(chunk_list, chunksize)

        with concurrent.futures.ProcessPoolExecutor(max_workers=5) as pool:
            
            future_list = {pool.submit(
                extract_f5_basemods, chunk_, bam_chunk.loc[chunk_], 
                ref.fasta.filename, chromosome=chromosome, 
                modification=modification, threshold=threshold, 
                context=context):
                    chunk_ for chunk_ in chunks}
            
            print("future list has {} members".format(len(future_list)))
            for future in tqdm(concurrent.futures.as_completed(future_list),
                               total=len(future_list)):
                
                bam_mapped = future.result()
                methylation_chunk = methylation_chunk.append(bam_mapped, sort=False)
                
            #     
            #     methylation_chunk = methylation_chunk.append(
            #         bam_mapped, sort=False)
                
            # and harvest the items ...
                
        
        # for item in bam_chunk.index:
        #     row = bam_chunk.loc[item]
        #     f5file = f5index.loc[item]["fast5file"]
        #     extract_f5basemods(
        #         f5file, row, chromosome, ref, modification=modification,
        #         threshold=threshold, context=context)
        
        print(methylation_chunk)
        sys.exit(0)
            
        
        if "flounder" in globals():
            flounder.write_cache(bam.bam, pd.DataFrame(), "methyl_chunk",
                                 chromosome, start, end)
    return methylation_chunk



def map_methylation(ref, bam, f5path):
    print("Mapping methylation ...")
    print(ref)
    print(bam)
    print(ref.get_chromosome_lengths(["20:"])[0])

    for chromo in ref.get_chromosomes():
        print(chromo)
        chromosome_tiles = ref.get_tiled_chromosome(chromo, tile_size=1000000)
        f5index = index_fast5_content(f5path)
        f5index.set_index("read_id", drop=False, inplace=True)
        for index in chromosome_tiles.df.index:
            print(index)
            mapped_reads = map_methylation_chunk(ref, bam, f5index,
                chromosome=chromo,
                start=chromosome_tiles.df.iloc[index].Start,
                end=chromosome_tiles.df.iloc[index].End)
            print(mapped_reads)
           
                    
                


if __name__ == '__main__':
    flounder = Flounder()
    f5path = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/"
    # f5file = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/A1-D1-PAD851010.fast5"
    # print(fast5_to_basemods(f5file))
    print(fast5s_to_basemods(f5path))
    sys.exit(0)
    reference = "/Volumes/Samsung_T5/MethylationPyTutorial/ReferenceData/human_g1k_v37.chr20.fasta"
    reference = ReferenceGenome(reference)
    
    bam = "/Volumes/Samsung_T5/MethylationPyTutorial/Analysis/minimap2/Native.bam"
    bam = BamHandler(bam)
    
    pd.set_option("display.max_columns", None)
    pd.set_option("display.expand_frame_repr", False)
    pd.set_option("max_colwidth", -1)
    
    map_methylation(reference, bam, f5path)
    
    
    
    