import os
from Bio import SeqIO

class SequenceHandler:

    def __init__(self, src):
        self.src = src
        self.handle = None

    def get_next_sequence(self):
        if not self.is_file_open() and not self.eof():
            self.open()
        return 1

    def is_file_open(self):
        return False

    def open(self):
        if self.get_file_type() in [".fasta", ".fa"]:
            print("opening fasta")
            self.handle = open(self.src, "rU")
            print(self.handle)
        else:
            raise Exception('[ {} ] is not a valid sequence format'.format(self.get_file_type()))

    def eof(self):
        return False

    def get_file_type(self):
        return os.path.splitext(self.src)[1].lower()


class Sequence:

    def __init__(self):
        self.sequence = ""


if __name__ == '__main__':
    seq_hand = SequenceHandler("/home/SRudd/Downloads/EPI_ISL_402119.fasta")
    print(seq_hand.get_file_type())
    seq = seq_hand.get_next_sequence()
    print(seq)