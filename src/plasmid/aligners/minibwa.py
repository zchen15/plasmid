#!/usr/bin/env python

# computing libraries
import pandas as pd

# system libraries
import subprocess
import tempfile

from .minimap2 import read_PAF, read_SAM


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with minibwa.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the SAM or PAF output into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        binary: minibwa binary.
        config: minibwa configuration string for alignment.
        index_config: minibwa configuration string for index building.
        config_meth: minibwa configuration string for bisulfite mode.
        config_hic: minibwa configuration string for Hi-C mode.
        meth: enable bisulfite alignment mode when True.
        hic: enable Hi-C alignment mode when True.
    '''
    verbose = False
    binary = 'minibwa'
    config = '-t8'
    index_config = '-t8'
    config_meth = '--meth'
    config_hic = '--hic'
    meth = False
    hic = False

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='minibwa_')

    # ------------------------------------------------------------------
    # minibwa index and alignment
    # ------------------------------------------------------------------
    def get_index(self, dbfile: str, fn_index: str=None) -> str:
        '''
        This is a wrapper for minibwa to build the fn_index

        Args:
            dbfile: fasta file of sequences to index.
            fn_index: index prefix. Defaults to the database filename, which is
                the minibwa default.

        Returns:
            A string containing filename of the minibwa index prefix.
        '''
        minibwa = self.binary
        if fn_index is None:
            fn_index = dbfile
        options = str(self.index_config).split(' ') if str(self.index_config) else []
        cmd = [minibwa, 'index'] + options + [dbfile, fn_index]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minibwa index failed. '+str(ret))
        self.fn_index = fn_index
        return fn_index

    def run(self, qfile, mfile: str='', dbfile: str=None,
            sam_file: str=None, paf_file: str=None) -> pd.DataFrame:
        '''
        This is a wrapper for the minibwa aligner.

        Args:
            qfile: fasta file of query sequences to search for in database.
            mfile: fasta file of the reverse pair of paired end reads.
            dbfile: fasta file or index prefix of the database sequences to
                search through. Defaults to the index built with get_index.
            sam_file: filename of the sam file to write and read results from.
            paf_file: filename of the paf file to write and read results from.

        Returns:
            a pandas dataframe of alignment results from the sam or paf file.

        Raises:
            ValueError: if no database has been built and no database is given.
        '''
        work_dir = self.work_dir.name
        if dbfile is None:
            if self.fn_index is None:
                raise ValueError('database is required when no index has been built')
            dbfile = self.fn_index
        config = str(self.config).split(' ') if str(self.config) else []
        if self.meth:
            config += str(self.config_meth).split(' ') if str(self.config_meth) else []
        if self.hic:
            config += str(self.config_hic).split(' ') if str(self.config_hic) else []
        minibwa = self.binary

        # do paired end search if needed
        if mfile == '':
            # single end or long read alignment with minibwa outputs a PAF file
            if paf_file is None:
                paf_file = work_dir + '/results.paf'
            cmd = [minibwa, 'map', '-f'] + config + [dbfile, qfile, '-o', paf_file]
            ret = subprocess.run(cmd)
            if ret.returncode != 0:
                raise RuntimeError('minibwa failed. '+str(ret))
            data = read_PAF(paf_file)
        else:
            # paired end alignment with minibwa outputs a SAM file
            if sam_file is None:
                sam_file = work_dir + '/results.sam'
            cmd = [minibwa, 'map'] + config + [dbfile, qfile, mfile, '-o', sam_file]
            ret = subprocess.run(cmd)
            if ret.returncode != 0:
                raise RuntimeError('minibwa failed. '+str(ret))
            data = read_SAM(sam_file)
        return data
