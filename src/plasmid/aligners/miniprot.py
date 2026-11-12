#!/usr/bin/env python

# computing libraries
import pandas as pd

# system libraries
import os
import subprocess
import tempfile

from .minimap2 import read_PAF


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with miniprot.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the PAF output of miniprot into a pandas dataframe.

    Attributes:
        verbose: print run settings when True.
        binary: miniprot binary.
        config: miniprot configuration string for alignment.
        index_config: miniprot configuration string for index building.
        fn_index: filename of the cached miniprot index.
    '''
    verbose = False
    binary = 'miniprot'
    config = '-I'
    index_config = '-t8'
    fn_index = None

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='miniprot_')

    # ------------------------------------------------------------------
    # miniprot database and alignment
    # ------------------------------------------------------------------
    def get_index(self, dbfile: str, fn_index: str=None) -> str:
        '''
        Build a miniprot index from a genome fasta file.

        Args:
            dbfile: fasta file of the genome to index.
            fn_index: filename of the index. Defaults to an index in the
                working directory.

        Returns:
            A string containing filename of the miniprot index.

        Raises:
            RuntimeError: if the miniprot index command fails.
        '''
        if fn_index is None:
            fn_index = os.path.join(self.work_dir.name, 'index.mpi')
        options = str(self.index_config).split(' ') if str(self.index_config) else []
        cmd = [self.binary] + options + ['-d', fn_index, dbfile]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('miniprot index failed. '+str(ret))
        self.fn_index = fn_index
        return fn_index

    def run(self, qfile: str, dbfile: str=None, paf_file: str=None) -> pd.DataFrame:
        '''
        Wrapper for performing sequence alignment with miniprot.

        Args:
            qfile: amino acid sequences to search for in the database.
            dbfile: genome fasta file or pre-built miniprot index. Defaults to
                the index built with get_index.
            paf_file: filename of the PAF alignment file.

        Returns:
            A pandas dataframe of alignment results from the PAF file.

        Raises:
            ValueError: if no index has been built and no database is given.
            RuntimeError: if the miniprot command fails.
        '''
        if dbfile is None:
            if self.fn_index is None:
                raise ValueError('database is required when no index has been built')
            dbfile = self.fn_index
        if paf_file is None:
            paf_file = os.path.join(self.work_dir.name, 'results.paf')

        # call the aligner
        options = str(self.config).split(' ') if str(self.config) else []
        cmd = [self.binary] + options + [dbfile, qfile]
        ret = subprocess.run(cmd, stdout=subprocess.PIPE)
        if ret.returncode != 0:
            raise RuntimeError('miniprot failed. '+str(ret))

        # write paf output from stdout
        with open(paf_file, 'w') as f:
            f.write(ret.stdout.decode('utf-8'))

        # extract the PAF file information
        data = read_PAF(paf_file)
        return data
