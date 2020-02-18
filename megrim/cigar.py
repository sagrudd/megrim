#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:02:04 2020

@author: srudd
"""

import pandas as pd
import numpy as np


tuple_map = pd.DataFrame({"op": ["M", "I", "D", "N", "S", "H", "P", "=", "X"],
                          "query": [1, 1, 0, 0, 1, 0, 0, 1, 1],
                          "reference": [1, 0, 1, 1, 0, 0, 0, 1, 1]},
                         index=[0, 1, 2, 3, 4, 5, 6, 7, 8])


def build_position_map(cigar_tuple):
    """
    Prepare a matrix of positions from provided pysam cigar_tuple.

    Parameters
    ----------
    cigar_tuple: tuple of operations and lengths
        This is expected to be produced by pysam package

    Returns
    -------
    pandas DataFrame
        With columns for downstream processing ...
    """
    if isinstance(cigar_tuple, str):
        cigar_tuple = eval(cigar_tuple)

    feature = [i[0] for i in cigar_tuple]
    repeats = [i[1] for i in cigar_tuple]

    source = np.repeat(feature, repeats)
    key = tuple_map.iloc[source]["op"]
    q_type = tuple_map.iloc[source]["query"]
    r_type = tuple_map.iloc[source]["reference"]
    q_posn = np.cumsum(q_type)
    r_posn = np.cumsum(r_type)

    return pd.DataFrame({"source": source,
                         "operation": key,
                         "q_type": q_type,
                         "r_type": r_type,
                         "q_posn": q_posn,
                         "r_posn": r_posn})


def cigar_q_to_r(poi, cigar_tuple):
    """
    Given a coordinate of interest and a cigar_tuple return reference_coord.

    Parameters
    ----------
    poi: int
        position of interest.
    cigar_tuple: list of tuple
        The tuple object as prepared by pysam.

    Raises
    ------
    ValueError
        It is possible that there are some features beyond coded - hopefully
        this will help me catch the exceptions.

    Returns
    -------
    cigar operation
        operation at the query position of interest.
    int
        the corresponding coordinate on the reference genome.

    """
    pmap = build_position_map(cigar_tuple)
    matches = pmap[pmap.q_posn == poi]
    if len(matches.index) > 1:
        m_matches = matches[matches.operation == "M"]
        if len(m_matches.index) != 1:
            print(matches)
            raise ValueError("The SAM coordinates are beyond experience")
        idx = m_matches.index[0]
        return m_matches["operation"][idx], m_matches["r_posn"][idx]
    idx = matches.index[0]
    return matches["operation"][idx], matches["r_posn"][idx]



class cigar_rle:

    def __init__(self, cigar_tuple):
        self.pmap = build_position_map(cigar_tuple)

    def q_to_r(self, poi):
        matches = self.pmap[self.pmap.q_posn == poi]
        if len(matches.index) > 1:
            m_matches = matches[matches.operation == "M"]
            # we may only have INS / DEL ???
            if len(m_matches.index) > 1:
                print(matches)
                raise ValueError("The SAM coordinates are beyond experience")
            elif len(m_matches.index) == 1:
                idx = m_matches.index[0]
                return m_matches["operation"][idx], m_matches["r_posn"][idx]
        idx = matches.index[0]
        return matches["operation"][idx], matches["r_posn"][idx]
        
        

