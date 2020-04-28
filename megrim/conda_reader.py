import logging
from megrim.environment import Flounder
import argparse


class CondaGit(Flounder):

    def __init__(self, args):
        Flounder.__init__(self)

        if args is not None:
            if isinstance(args, argparse.Namespace):
                self.argparse(args)
            if isinstance(args, dict):
                self.dictparse(args)

    def lookup(self):
        logging.info("looking up a conda entity ...")