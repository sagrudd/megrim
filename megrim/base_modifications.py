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
import concurrent.futures
from tqdm import tqdm
import multiprocessing


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

            mod_base_df = pd.DataFrame(
                data=read.get_analysis_dataset(
                    latest_basecall,
                    "BaseCalled_template/ModBaseProbs"),
                columns=["A", "6mA", "C", "5mC", "G", "T"])
            mod_base_df['position'] = mod_base_df.index
            mod_base_df['read_id'] = read_id

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

            result = result.append(mod_base_df, sort=False)

    if "flounder" in globals():
        flounder.write_cache(
            fast5file, result, modification, threshold, context)
    return result


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


if __name__ == '__main__':
    flounder = Flounder()
    f5path = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/"
    # f5file = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/A1-D1-PAD851010.fast5"
    # print(fast5_to_basemods(f5file))
    print(fast5s_to_basemods(f5path))
    