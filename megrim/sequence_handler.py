import os
from Bio import SeqIO
import pandas as pd
import hashlib
from math import log10
from dateutil.parser import parse

from megrim.basic_qc import SequenceSummaryHandler
from megrim.reference_genome import ReferenceGenome

class SequenceHandler:

    def __init__(self, src):
        self.src = src
        self.handle = None
        self.cache = None

    def has_next(self):
        try:
            self.cache = self.get_next_sequence()
            return True
        except StopIteration:
            return False

    def get_next_sequence(self):
        if self.cache is None:
            if not self.is_file_open() and not self.eof():
                self.open()
            if self.get_file_type() in [".fasta", ".fa", ".fastq", ".fq"]:
                return Sequence(next(self.handle))
        return self.get_cache()

    def get_cache(self):
        cache = self.cache
        self.cache = None
        return cache

    def is_file_open(self):
        return self.handle != None

    def open(self):
        if self.get_file_type() in [".fasta", ".fa"]:
            print("opening fasta")
            self.handle = SeqIO.parse(open(self.src), "fasta")
            print(self.handle)
        elif self.get_file_type() in [".fastq", ".fq"]:
            print("opening fastq")
            self.handle = SeqIO.parse(open(self.src), "fastq")
            print(self.handle)
        else:
            raise Exception('[ {} ] is not a valid sequence format'.format(self.get_file_type()))

    def eof(self):
        return False

    def get_file_type(self):
        return os.path.splitext(self.src)[1].lower()

    def get_file_md5sum(self):
        """
        Get the defined file's md5sum

        The md5 checksum is a simple way to ensure that a file is a faithful copy of the original.
        This method calculates and returns the md5 checksum.

        Returns
        -------
        str
            A string of the md5sum.
        """
        return hashlib.md5(open(self.src, "rb").read()).hexdigest()

    def fastq_report(self, limit=-1):
        if self.get_file_type() not in [".fastq", ".fq"]:
            raise Exception('[ {} ] not valid - fastq required'.format(self.get_file_type()))
        counter = 0
        results = []
        while self.has_next():
            counter += 1
            seq = self.get_next_sequence()
            results.append(seq.get_nanopore_summary())
            if counter == limit:
                break

        target_data = pd.concat(
            results, axis=1).transpose()
        target_data = target_data.astype({'passes_filtering': 'bool'})
        #return target_data
        return SequenceSummaryHandler(target_data=target_data)


class Sequence:

    def __init__(self, record):
        self.record = record
        self.annotation = None

    def __str__(self):
        return str(self.record.id)

    def get_nt_rle(self, min=4):
        """
        Get a DataFrame of RunLengthEncoded nucleotide runs over threshold.

        Homopolymer sequences are naturally occuring and largely stochastic
        repeats within the sequence space. Frequency of occurence may be
        influenced by genome size and GC richness. This simple method consumes
        a sequence to return to return homopolymer runs over a minimum
        threshold size, min.

        Parameters
        ----------
        min: int
            The minimum number of nucleotides that should be observed within
            a run for it to be reported. The default value is 4

        Returns
        -------
        pd.DataFrame
            A pandas dataframe ordered by position of observed repeat run.

        """
        offset = 0
        run = None
        run_len = 0
        run_start = None
        runs = []

        def append_run():
            if run_len >= min:
                series = pd.Series(
                    (run_start, run, run_len),
                    index=["position", "run", "length"])
                runs.append(series)

        for c in str(self.record.seq):
            offset += 1
            if run is None or run != c:
                append_run()
                run = c
                run_len = 1
                run_start = offset
            else:
                # extend
                run_len += 1
        # and the residuals
        append_run()
        return pd.concat(
            runs, axis=1, keys=[s.position for s in runs]).transpose()

    def get_longest_runs(self, n=5, min=4):
        """
        Get a summary of longest homopolymer runs observed within a sequence.

        This method subsets the homopolymer repeats reported by
        sequence_handler.Sequence.get_nt_rle to report the longest
        homopolymer repeats.

        Parameters
        ----------
        n: int
            The number of repeat items to display. The default is 5.
        min: int
            The minimum number of nucleotides that should be observed within
            a run for it to be reported. The default value is 4

        Returns
        -------
        pd.DataFrame
            A pandas dataframe ordered by position of observed repeat run.

        """
        df = self.get_nt_rle(min)
        print(df.sort_values(
            by=["length", "position"], ascending=False).head(n))

    def get_read_mean_quality(self):
        return -10 * log10((10 ** (pd.Series(self.record.letter_annotations["phred_quality"]) / -10)).mean())

    def normalise_start_time(self, time):
        #return np.datetime64(time).astype(int)
        return int(round(parse(time).timestamp()))

    def load_annotations(self):
        items = self.record.description.split(" ")
        self.annotation = {}
        for item in items:
            if "=" in item:
                key, val = item.split("=", maxsplit=1)
                if key == "start_time":
                    val = self.normalise_start_time(val)
                self.annotation[key]=val

    def get_annotation(self, key):
        if self.annotation is None:
            self.load_annotations()
        return self.annotation[key]

    def get_nanopore_summary(self):
        # returning id, sequence_length, mean_q,
        series = pd.Series(
            (self.record.id, len(str(self.record.seq)),
             self.get_read_mean_quality(), self.get_annotation("ch"),
             self.get_annotation("start_time"), True
             ),
            index=["id", "sequence_length_template", "qual", "channel", "start_time", "passes_filtering"])
        return series
