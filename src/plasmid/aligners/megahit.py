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


def read_contigs(fname: str) -> pd.DataFrame:
    '''
    Read assembled contig fasta file into a dataframe of contigs.

    Args:
        fname: filename of the assembled contigs fasta file.

    Returns:
        A dataframe of assembled contigs with columns [name, sequence, length].
    '''
    df = io.read_fasta(fname)
    df['name'] = [i.split(' ')[0] for i in df['name']]
    df['length'] = [len(i) for i in df['sequence']]
    return df.sort_values(by=['length'], ascending=False)


class Aligner:
    '''
    This class holds functions pertaining to sequence assembly with megahit.

    The class takes fasta and fastq filenames as input. Assembly results are
    parsed from the assembled contig fasta file into a pandas dataframe.

    Attributes:
        verbose: print run settings when True.
        binary: megahit binary.
        config: megahit configuration string.
    '''
    verbose = False
    binary = 'megahit'
    config = ''

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='megahit_')

    # ------------------------------------------------------------------
    # megahit assembly
    # ------------------------------------------------------------------
    def run(self, qfile, mfile: str='', infile: str=None) -> pd.DataFrame:
        '''
        Wrapper for running metagenome assembly with megahit.

        Args:
            qfile: fastq file of reads.
            mfile: fastq file of mate reads.
            infile: fasta file of assembled contigs to read into a dataframe.

        Returns:
            A pandas dataframe of assembled contigs with columns [name, sequence, length].

        Raises:
            RuntimeError: if the megahit command fails.
        '''
        work_dir = self.work_dir.name
        if infile is None:
            # run megahit
            if mfile == '':
                cmd = [self.binary] + str(self.config).split(' ') + ['-r', qfile, '-o', work_dir + '/results/']
            else:
                cmd = [self.binary] + str(self.config).split(' ') + ['-1', qfile, '-2', mfile, '-o', work_dir + '/results/']
            ret = subprocess.run(cmd)
            if ret.returncode != 0:
                raise RuntimeError('megahit failed. '+str(ret))
            infile = work_dir + '/results/final.contigs.fa'

        # read assembled contigs
        data = read_contigs(infile)
        return data
