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
from .minimap2 import read_PAF


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with miniprot.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the PAF output of miniprot into a pandas dataframe.

    Attributes:
        verbose: print run settings when True.
        binary: miniprot binary.
        config: miniprot configuration string.
        fn_index: filename of the cached miniprot database.
    '''
    verbose = False
    binary = 'miniprot'
    config = ''
    fn_index = None

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
        self.work_dir = tempfile.TemporaryDirectory(prefix='miniprot_')

    # ------------------------------------------------------------------
    # miniprot database and alignment
    # ------------------------------------------------------------------
    def get_index(self, dbfile: str) -> str:
        '''
        This is a wrapper for miniprot to build the fn_index

        Args:
            dbfile: fasta file of amino acid sequences to index.

        Returns:
            A string containing filename of the miniprot database.
        '''
        self.fn_index = dbfile
        return dbfile

    def run(self, qfile, dbfile) -> pd.DataFrame:
        '''
        Wrapper for performing sequence alignment with miniprot.

        Args:
            qfile: amino acid sequences to search for in the database.
            dbfile: fasta file of amino acid sequences of the database.

        Returns:
            A pandas dataframe of alignment results from the PAF file.

        Raises:
            RuntimeError: if the miniprot command fails.
        '''
        work_dir = self.work_dir.name
        out_file = work_dir + '/results.paf'

        # call the aligner
        if len(str(self.config)) > 0:
            cmd = [self.binary] + str(self.config).split(' ') + [dbfile, qfile]
        else:
            cmd = [self.binary, dbfile, qfile]
        ret = subprocess.run(cmd, stdout=subprocess.PIPE)
        if ret.returncode != 0:
            raise RuntimeError('miniprot failed. '+str(ret))

        # write paf output from stdout
        with open(out_file, 'w') as f:
            f.write(ret.stdout.decode('utf-8'))

        # extract the PAF file information
        data = read_PAF(out_file)
        return data
