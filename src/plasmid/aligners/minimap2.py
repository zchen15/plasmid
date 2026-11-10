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
from .bwa import read_SAM


def format_alignment(data):
    '''
    Format numeric columns and add match_score metrics.

    Args:
        data: dataframe of alignment results.

    Returns:
        A formatted dataframe.
    '''
    int_cols = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end',
                'match', 'tot', 'mapq', 'cm', 's1', 's2', 'NM', 'AS',
                'ms', 'nn', 'rl']
    for c in int_cols:
        if c in data.columns:
            data[c] = pd.to_numeric(data[c], errors='coerce').fillna(0).astype(int)
    if 'de' in data.columns:
        data['de'] = pd.to_numeric(data['de'], errors='coerce').fillna(0.0).astype(float)

    if 'match' in data.columns and 'tot' in data.columns:
        tot = pd.to_numeric(data['tot'], errors='coerce').replace(0, np.nan)
        data['match_score'] = (pd.to_numeric(data['match'], errors='coerce')/tot).fillna(0.0)
    elif 'match_score' not in data.columns:
        data['match_score'] = pd.Series(dtype=float)
    return data


def read_PAF(fname):
    '''
    Read a minimap2 PAF file.

    Args:
        fname: PAF filename.

    Returns:
        A pandas dataframe of alignment results.
    '''
    col1 = ['q_id', 'q_len', 'q_start', 'q_end', 'orientation',
            't_id', 't_len', 't_start', 't_end', 'match', 'tot', 'mapq']
    col2 = ['tp', 'cm', 's1', 's2', 'NM', 'AS', 'ms', 'nn', 'rl', 'cg']
    key = {col2[i]: i for i in range(len(col2))}

    out = []
    with io.decompress(fname) as f:
        for line in f:
            line = line.rstrip('\n')
            if len(line) == 0 or line[0] == '#':
                continue
            fields = line.split('\t')
            if len(fields) < len(col1):
                continue
            info = [0]*len(col2)
            for field in fields[len(col1):]:
                x = field.split(':', 2)
                if len(x) == 3 and x[0] in key:
                    info[key[x[0]]] = x[2]
            out.append(fields[:len(col1)]+info)

    data = pd.DataFrame(out, columns=col1+col2)
    data = data.rename(columns={'cg': 'CIGAR'})
    if 'CIGAR' in data.columns:
        data['CIGAR'] = [str(i).replace('\n', '') for i in data['CIGAR']]
    return format_alignment(data)


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

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='minimap2_')

    # ------------------------------------------------------------------
    # minimap2 index and alignment
    # ------------------------------------------------------------------
    def get_index(self, database: str, fn_index: str=None):
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

    def run(self, qfile, mfile=None, dbfile=None, paf_file=None, sam_file=None):
        '''
        Wrapper for performing sequence alignment with minimap2.

        Args:
            qfile: filename of the query sequences (fasta or fastq).
            mfile: filename of the paired end query sequences (fasta or fastq).
            dbfile: filename of the minimap2 index of the database. Required if no index
                has been built with get_index.
            paf_file: filename of the PAF alignment file.
            sam_file: filename of the SAM alignment file.

        Returns:
            A pandas dataframe of alignment results from the PAF or SAM file.
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

        return data
