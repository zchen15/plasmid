#!/usr/bin/env python

# computing libraries
import numpy as np
import pandas as pd

# system libraries
import os
import re
import subprocess
import tempfile

from . import io


def read_fasta(fname: str) -> pd.DataFrame:
    '''
    Read polished sequences from a fasta file.

    Args:
        fname: filename of the fasta file of polished sequences.

    Returns:
        A dataframe of sequences with columns [name, sequence].
    '''
    return io.read_fasta(fname)


class Aligner:
    '''
    This class holds functions pertaining to polishing reads with raco.

    The class takes fasta and fastq filenames as input. Polished sequences are
    parsed from the polished fasta file into a pandas dataframe.

    Attributes:
        verbose: print run settings when True.
        binary: raco binary.
        config: raco configuration string.
    '''
    verbose = False
    binary = 'racon'
    config = '-w 500 -q 10'

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        Raises:
            TypeError: if kwargs contains a name that is not an Aligner setting.
        '''
        for key, value in kwargs.items():
            if not hasattr(self, key):
                continue
            if value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='racon_')

    # ------------------------------------------------------------------
    # raco polishing
    # ------------------------------------------------------------------
    def run(self, seq, paf, ref, config: str=None, out_file: str=None) -> pd.DataFrame:
        '''
        Wrapper for polishing reads against a reference alignment with raco.

        Args:
            seq: sequencing reads fasta or fastq file.
            paf: alignment file for reads against the draft reference sequence.
            ref: reference sequence fasta or fastq file.
            config: config options to pass to raco. Defaults to the raco config.
            out_file: filename of the polished fasta file.

        Returns:
            A pandas dataframe of polished sequences with columns [name, sequence].

        Raises:
            RuntimeError: if the raco command fails.
        '''
        work_dir = self.work_dir.name
        if out_file is None:
            out_file = work_dir + '/results.fa'
        if config is None:
            config = self.config

        # run raco
        if len(str(config)) > 0:
            cmd = [self.binary] + str(config).split(' ') + [seq, paf, ref]
        else:
            cmd = [self.binary, seq, paf, ref]
        ret = subprocess.run(cmd, stdout=subprocess.PIPE)
        if ret.returncode != 0:
            raise RuntimeError('racon failed. '+str(ret))

        # write polished fasta output from stdout
        with open(out_file, 'w') as f:
            f.write(ret.stdout.decode('utf-8'))

        # read the polished sequences
        data = read_fasta(out_file)
        return data
