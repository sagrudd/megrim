import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

class SequenceHandler:

    def __init__(self, src):
        self.src = src
        self.handle = None

    def get_next_sequence(self):
        if not self.is_file_open() and not self.eof():
            self.open()
        if self.get_file_type() in [".fasta", ".fa"]:
            return Sequence(next(self.handle))

    def is_file_open(self):
        return self.handle != None

    def open(self):
        if self.get_file_type() in [".fasta", ".fa"]:
            print("opening fasta")
            self.handle = SeqIO.parse(open(self.src), "fasta")
            print(self.handle)
        else:
            raise Exception('[ {} ] is not a valid sequence format'.format(self.get_file_type()))

    def eof(self):
        return False

    def get_file_type(self):
        return os.path.splitext(self.src)[1].lower()


class Sequence:

    def __init__(self, record):
        self.record = record

    def __str__(self):
        return str(self.record.seq)

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


if __name__ == '__main__':
    seq_hand = SequenceHandler("/Users/srudd/SynologyDrive/nCoV/EPI_ISL_402119.fasta")
    print(seq_hand.get_file_type())

    seq = seq_hand.get_next_sequence()
    print(seq)
    seq.get_longest_runs(n=15)
    