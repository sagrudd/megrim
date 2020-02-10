#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 20:05:43 2020

@author: stephen
"""


from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import SeqIO
from Bio.Seq import Seq
import math
import pandas as pd


class genbank:
    """
    Class for handling genbank parsing and feature extraction.

    I may be creating too many classes; but I find classes useful. Genbank
    sequences contain a load of useful; this class aims to simplify the
    process of pulling annotation and derived context from Genbank flatfiles.
    """

    INTERGENIC = "intergenic"

    def __init__(self, gbk_file):
        self.gbk_file = gbk_file

    def get_thing(self, thing):
        """
        Get the Genbank feature for a named gene's CDS.

        This method iterates through Genbank annotations and will return
        the Biopython feature object for the named gene.

        Parameters
        ----------
        thing: String
            The name of the gene to parse features for.

        Returns
        -------
        feat: Biopython feature
            Biopython feature.

        """
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                if feat.qualifiers['gene'] == thing:
                    return feat
        return None

    def get_sequence(self):
        """
        Get the Biopython format DNA sequence for the given Genbank object.

        Returns
        -------
        Sequence
            In Biopython format.

        """
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            return rec.seq
        return None

    def get_translation(self, thing):
        """
        Return the translated CDS from a named gene feature.

        Parameters
        ----------
        thing: String
            The name of the gene's CDS.

        Returns
        -------
        Translation
            Peptide sequence as Biopython sequence object.

        """
        return self.get_transcript(thing).translate()

    def get_transcript(self, thing):
        """
        Prepare the DNA transcript sequence for the named gene feature.

        Parameters
        ----------
        thing: String
            The name of the gene's CDS.

        Raises
        ------
        ValueError
            Methods have been implemented for FeatureLocation and
            CompoundLocation. I am not sure if there could be edge cases
            with a different format.

        Returns
        -------
        Transcript
            Transcript DNA sequence in Biopython Seq format.

        """
        features = self.get_thing(thing)
        if isinstance(features.location, FeatureLocation):
            return features.location.extract(self.get_sequence())
        elif isinstance(features.location, CompoundLocation):
            return features.location.extract(self.get_sequence())
        else:
            raise ValueError("The method for ?? transcripts is missing")

    def get_position_in_transcript(self, pos):
        """
        Given a genome, return the position within an underlying transcript.

        This is an accessory method. When provided with a genome coordinate
        the method will identify the CDS feature overlapping this position
        and will return the position of the given nucleotide within the
        transcript sequence.

        Parameters
        ----------
        pos: int
            The position to look for assess for CDS and CDS sequence coords.

        Raises
        ------
        ValueError
            If there is no annotated CDS at the given position, ValueError.

        Returns
        -------
        results: int
            The position of the given base within the CDS.

        """
        thing = genbank.what_is_at(pos)
        if thing == genbank.INTERGENIC:
            return None
        features = self.get_thing(thing)
        if features.type == 'CDS':
            results = []
            if features.location.strand < 0:
                raise ValueError("This hasn't been tested on minus strand")
            if isinstance(features.location, FeatureLocation):
                results.append(pos - (features.location.start))
            elif isinstance(features.location, CompoundLocation):
                # we are working over multiple intervals - maintain an offset
                offset = 0
                for part in features.location.parts:
                    if (pos >= part.start) & (pos <= part.end):
                        results.append(pos - part.start + offset)
                    offset += len(part)
            else:
                raise ValueError("The method for ?? transcripts is missing")
            return results
        raise ValueError(
            "get_position_in_transcript {!CDS} missing")

    def get_residue_in_genome(self, position):
        """
        Get the canonical base at given position from reference genome.

        Simple method; gets the given base for QC purposes mainly. Please
        note that Python is 0-based whilst genbank entries are 1-based.

        Parameters
        ----------
        position: int
            The base position (1-based) to get from the reference genome.

        Returns
        -------
        String
            Of a single letter.

        """
        return self.get_sequence()[position-1]

    def get_triplet(self, pos):
        """
        Given a position in a genome, return a codon triplet for CDS.

        If we have a genome coordinate from e.g. a VCF file this method will
        quickly pull the corresponding codon (and base position in codon) if
        the position overlaps with an annotated CDS.

        This method will return a tuple of None if the given position does
        not overlap with an annotated CDS.

        Parameters
        ----------
        pos: int
            The position within thge genome.

        Returns
        -------
        triplet
            Codon sequence as Biopython object containing the codon which
            overlaps the given CDS.
        trippos
            0-based index of given pos within the codon sequence returned.

        """
        thing = genbank.what_is_at(pos)
        if thing == genbank.INTERGENIC:
            return None, None
        features = self.get_thing(thing)
        if features.type == 'CDS':
            # print("we have CDS ...")
            transcript = self.get_transcript(thing)
            positions = self.get_position_in_transcript(pos)
            # prepare the triplets ...
            chunk_size = 3
            triplets = [
                transcript[i:i+chunk_size] for i in range(
                    0, len(transcript), chunk_size)]

            # identify the triplet to process
            triplet = triplets[math.floor((positions[0]-1)/3)]
            # and get the position in the triplet ...
            tripos = positions[0]-1 - math.floor((positions[0]-1)/3)*3
            return triplet, tripos

    def get_variant_triplet(self, pos, base):
        """
        Get the codon sequence for a modified base.

        This method pulls out the wild-type codon from the CDS overlapping
        given position pos, and mutates the appropriate codon position to
        give the variant codon.

        Parameters
        ----------
        position: int
            The genomic location of the base being considered.
        base: character
            The variant observed at position pos.

        Returns
        -------
        Codon
            Biopython Seq format codon for mutated triplet at position pos.

        """
        triplet, tripos = self.get_triplet(pos)
        triplet = list(triplet)
        triplet[tripos] = base
        return Seq("".join(triplet))

    def what_is_at(self, position):
        """
        Identify whether given genomic position contains a gene.

        This method iterates through genomic features to see whether a gene
        overlaps the given coordinate. If a CDS is found the correponding
        gene name is returned. Is the coordinate is intergenic the class
        constant INTERGENIC is returned instead.

        Parameters
        ----------
        position: int
            Position within the genome of interest.

        Returns
        -------
        Str
            Gene name of overlapping gene or constant defined by INTERGENIC.

        """
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                if (position >= feat.location.start) & (
                        position <= feat.location.end):
                    return feat.qualifiers['gene']
        return genbank.INTERGENIC

    def get_variant_effect(self, position, variant):
        """
        Get variant effect summary for genetic variant at defined position.

        This is an accessory method for reporting in human-readable form the
        probable variant effect of the observed variant, variant, at genomic
        position position. If the given position occurs within a CDS the
        variant effect will be returned as synonymous, non-synonymous or
        nonsense depending on the observed CDS and AA change.

        Parameters
        ----------
        position: int
            The location within the reference genome.
        variant: character
            The variant nucleotide observed at position position.

        Returns
        -------
        Str
            A human readable summary describing the observed effect. If the
            location is intergenic a None value will be returned.

        """
        triplet, tripos = genbank.get_triplet(position)
        if triplet is None:
            return None
        mtriplet = genbank.get_variant_triplet(position, variant)
        if triplet.translate() == mtriplet.translate():
            return "Synonymous [{}]".format(triplet.translate())
        elif mtriplet.translate() == "*":
            return "Nonsense [{}->{}]".format(
                triplet.translate(), mtriplet.translate())
        else:
            return "Non-synonymous [{}->{}]".format(
                triplet.translate(), mtriplet.translate())

    def print_mutation(self, row):
        """
        Summarise information for a SNP defined within a pd.DataFrame.

        This method is aimed to assist in the reporting of SNP meanings in
        a whole genome context. A pd.DataFrame containing columns named
        'positions' and 'variants' should be prepared. This method will
        consume a row, as within a pd.DataFrame.apply call and will return
        Pd.Series containing summary information for the variant defined
        within.

        Parameters
        ----------
        row: pd.Series
            A row from a pd.DataFrame.

        Returns
        -------
        mutation: pd.Series
            A pd.Series containing annotation for variant defined in row.

        """
        position = row.positions
        variant = row.variants
        mutation = pd.Series(
            {"position": position,
             "change": "{}->{}".format(
                 self.get_residue_in_genome(position), variant),
             "location": self.what_is_at(position),
             "effect": self.get_variant_effect(position, variant)})
        return mutation


if __name__ == "__main__":

    genbank = genbank("/Users/srudd/Downloads/sequence.gb")
    # genbank.to_ranges()

    starting_data = pd.DataFrame({"positions": [100, 19065, 22303, 26144],
                                  "variants": ["A", "C", "G", "T"]})

    print(pd.DataFrame(
        starting_data.apply(
            genbank.print_mutation, axis=1, result_type='expand')))
