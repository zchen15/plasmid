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

class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with minimap2.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the PAF, SAM, or BAM output into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        config: minimap2 configuration string.
        index: filename of the cached minimap2 index.
        binary: minimap2 binary.
    '''
    verbose = False
    binary = 'minimap2'
    config = '-c -k 15 -w 7'
    fn_index = None

    # alignment file parsing is implemented in aligners.io
    read_PAF = staticmethod(io.read_PAF)
    read_paf = staticmethod(io.read_PAF)
    read_SAM = staticmethod(io.read_SAM)
    read_BAM = staticmethod(io.read_BAM)

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if key not in settings:
                raise TypeError('unknown Aligner setting: '+key)
            setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='minimap2_')

    # ------------------------------------------------------------------
    # minimap2 index and alignment
    # ------------------------------------------------------------------
    def get_mmi_index(self, database:str, fn_index:str=None):
        '''
        Build a minimap2 index.

        Args:
            database: filename containing the reference sequences.
            filename: name of the index file inside the working directory.

        Returns:
            Filename of the minimap2 index.
        '''
        if fn_index is None:
            fn_index = os.path.join(self.work_dir.name, "index.mmi")
        
        cmd = [self.binary,
               '-k', str(self.k),
               '-w', str(self.w)]
        if self.threads is not None:
            cmd+= ['-t', str(self.threads)]
        cmd+= ['-d', fn_index, database]
        
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minimap2 index failed. '+str(ret))

        self.fn_index = fn_index
        return fn_index

    def minimap2(self, qfile, mfile=None, dbfile=None, paf_file=None, sam_file=None):
        '''
        Wrapper for performing sequence alignment with minimap2.
        '''
        work_dir = self.work_dir.name
        if dbfile is None and self.fn_index is not None:
            dbfile = self.fn_index
        else:
            raise ValueError('database is required when no index has been built')
        
        config = self.config
        if sam_file is not None:
            config = f"-a {config}"
        elif paf_file is not None:
            outfile = paf_file
        else:
            outfile = os.path.join(self.work_dir.name, "results.paf")
        config = config.split(" ")
        
        # call the aligner
        if mfile is None:
            cmd = [self.binary]+config+[dbfile, qfile, '-o', outfile]
        else:
            cmd = [self.binary]+config+[dbfile, qfile, mfile, '-o', outfile]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minimap2 failed. '+str(ret))
        
        # parse the alignment output
        if sam_file is not None:
            data = read_SAM(outfile)
        else:
            data = read_PAF(outfile)
        
        # format output columns
        int_cols = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end',
                    'match', 'tot', 'mapq', 'cm', 's1', 's2', 'NM', 'AS',
                    'ms', 'nn', 'rl']
        for c in int_cols:
            if c in data.columns:
                data[c] = data[c].astype(int)
        return data
