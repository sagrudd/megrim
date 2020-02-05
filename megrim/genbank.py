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
                    print(feat)
                    return feat
        return None
    
    def get_sequence(self):
        recs = [rec for rec in SeqIO.parse(self.gbk_file, "genbank")]
        for rec in recs:
            return rec.seq
        return None
    
    def get_translation(self, thing):
        features = self.get_thing(thing)
        sequence = self.get_sequence()
        if isinstance(features.location, FeatureLocation):
            print("FeatureLocation")
            if features.location.strand == +1:
                print(sequence[features.location.start:features.location.end].translate())
        else:
            print(type(features))
    
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
    
    
if __name__ == "__main__":
    
    genbank = genbank("/Users/stephen/Downloads/sequence.gb")
    genbank.to_ranges()
    
    pos = 29560
    
    thing = genbank.what_is_at(pos)
    print(thing)
    if thing != genbank.INTERGENIC:
        features = genbank.get_translation(thing)
        