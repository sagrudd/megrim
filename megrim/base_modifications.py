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


def fast5_to_basemods(
        fast5file, latest_basecall=None, modification="5mC", threshold=0.75, 
        context="CG"):
    result = pd.DataFrame()
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
                mod_base_df["contexts"] = mod_base_df.apply(get_context, axis=1)

            result = result.append(mod_base_df)
    print(result)
            
    
    


if __name__ == '__main__':
    f5file = "/Volumes/Samsung_T5/MethylationPyTutorial/RawData/ONLL04465/fast5chr20_mods/workspace/A1-D1-PAD851010.fast5"
    fast5_to_basemods(f5file)