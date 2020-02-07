#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 20:05:43 2020

@author: stephen
"""


from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys
import math
import pandas as pd

class genbank:
    
    INTERGENIC = "intergenic"
    
    def __init__(self, gbk_file):
        self.gbk_file = gbk_file
    
    
    
    def to_ranges(self):
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            print(rec)
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                print(feat)
                print(feat.location)
                print(feat.qualifiers['gene'])
                print("\n\n\n")
    
    
    def get_thing(self, thing):
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                if feat.qualifiers['gene'] == thing:
                    return feat
        return None
    
    def get_sequence(self):
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            return rec.seq
        return None
    
    def get_translation(self, thing):
        return self.get_transcript(thing).translate()
    
    def get_transcript(self, thing):
        features = self.get_thing(thing)
        # if isinstance(features.location, FeatureLocation):
        #     print("FeatureLocation")
        #     if features.location.strand == +1:
        #         return sequence[features.location.start:features.location.end]
        #     else:
        #         raise ValueError("The method for reverse strand transcripts is missing")
        # elif isinstance(features.location, CompoundLocation):
        #     print("CompoundLocation")
        #     for fl in features.location:
        #         print(fl)
        #     raise ValueError("The method for Compound transcripts is missing")
        # else:
        #     print(type(features.location))
        if isinstance(features.location, FeatureLocation):
            return features.location.extract(self.get_sequence())
        elif isinstance(features.location, CompoundLocation):
            return features.location.extract(self.get_sequence())
        else:
            raise ValueError("The method for ?? transcripts is missing")
    
    def get_position_in_transcript(self, pos):
        thing = genbank.what_is_at(pos)
        if thing == genbank.INTERGENIC:
            return None
        features = self.get_thing(thing)
        if features.type == 'CDS':
            results = []
            if features.location.strand < 0:
                raise ValueError("This hasn't been tested on minus strand")
            if isinstance(features.location, FeatureLocation):
                results.append(pos - features.location.start)
            elif isinstance(features.location, CompoundLocation):
                print(features.location)
                # we are working over multiple intervals - maintain an offset
                offset = 0
                for part in features.location.parts:
                    if (pos >= part.start) & (pos <= part.end):
                        print(part)
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
        thing = genbank.what_is_at(pos)
        if thing == genbank.INTERGENIC:
            return None, None
        features = self.get_thing(thing)
        if features.type == 'CDS':
            print("we have CDS ...")
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

    def get_variant_triplet(self, position, base):
        triplet, tripos = self.get_triplet(position)
        triplet = list(triplet)
        triplet[tripos] = base
        return Seq("".join(triplet))
    
    
    def extract_coords(self, stuff):
        if isinstance(stuff, CompoundLocation):
            return stuff.start, stuff.end, stuff.strand
        elif isinstance(stuff, FeatureLocation):
            return stuff.start, stuff.end, stuff.strand
        return 1, 1, "*"
    
    def what_is_at(self, position):
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                start, stop, strand = self.extract_coords(feat.location)
                if (position >= start) & (position <= stop):
                    return feat.qualifiers['gene']
        return genbank.INTERGENIC
    
    
    def get_variant_effect(self, position, variant):
        triplet, tripos = genbank.get_triplet(position)
        if triplet is None:
            return None
        mtriplet = genbank.get_variant_triplet(position, variant)
        if triplet.translate() == mtriplet.translate():
            return "Synonymous [{}]".format(triplet.translate())
        elif mtriplet.translate() == "*":
            return "Nonsense [{}->{}]".format(triplet.translate(), mtriplet.translate())
        else:
            return "Non-synonymous [{}->{}]".format(triplet.translate(), mtriplet.translate())
        
    
    def print_mutation(self, row):
        position = row.positions
        variant = row.variants
        mutation = pd.Series({"position": position,
                              "change": "{}->{}".format(self.get_residue_in_genome(position), variant),
                              "location": self.what_is_at(position),
                              "effect": self.get_variant_effect(position, variant)})
        return mutation
    
if __name__ == "__main__":
    
    genbank = genbank("/Users/srudd/Downloads/sequence.gb")
    # genbank.to_ranges()
    
    
    starting_data = pd.DataFrame({"positions":[100, 19065, 22303, 26144],
                                  "variants": ["A", "C", "G", "T"]})
    
    
    print(pd.DataFrame(starting_data.apply(genbank.print_mutation, axis=1, result_type='expand')))
    

        