#!/usr/bin/env python

# computing libraries
import numpy as np
import pandas as pd

# BAM file parsing
import pysam

# system libraries
import os
import re
import subprocess
import tempfile

from . import io


def get_BAM_tags(tag, cols) -> list:
    '''
    Parse BAM tags for columns

    Args:
        tag: list of tuples of tag key and value
        cols: values to extract

    Returns:
        list of values
    '''
    out = [None]*len(cols)
    for j in range(len(cols)):
        for k in tag:
            if k[0] == cols[j]:
                out[j] = k[1]
                break
    return out


def read_SAM(fname, mode='rb') -> pd.DataFrame:
    '''
    Reads a bam or sam file

    Args:
        fname: filename of bam or sam file
        mode: 'r' for SAM or 'rb' for BAM

    Returns:
        A pandas dataframe of alignment information
    '''
    samfile = pysam.AlignmentFile(str(fname), mode=mode, check_sq=False)

    out = []
    col2 = ["NM","ms","AS","nn","tp","cm","s1","s2","de","rl"]
    for read in samfile:
        q_id = read.qname
        q_len = read.infer_read_length()
        q_start = read.qstart
        q_end = read.qend

        t_id = read.reference_name
        t_len = read.reference_length
        t_start = read.reference_start
        t_end = read.reference_end

        if t_id is not None:
            orientation = read.is_forward
            cigar = read.cigarstring
            mapq = read.mapq
            tags = read.tags
            [NM, ms, AS, nn, tp, cm, s1, s2, de, rl] = get_BAM_tags(tags, col2)
            out.append([q_id, q_len, q_start, q_end, orientation,
                        t_id, t_len, t_start, t_end, cigar, mapq,
                        NM, ms, AS, nn, tp, cm, s1, s2, de, rl])
    samfile.close()
    col1 = ["q_id","q_len","q_start","q_end","orientation",
            "t_id","t_len","t_start","t_end","CIGAR","mapq"]
    cols = col1 + col2
    df = pd.DataFrame(out, columns=cols)

    c = df["orientation"].values
    df.loc[c, "orientation"] = "+"
    df.loc[~c, "orientation"] = "-"

    # compute blast like score
    df["tot"] = df["q_end"] - df["q_start"]
    df["match"] = df["tot"] - df["NM"]
    df["match_score"] = df["match"]/df["tot"]
    return df


def format_alignment(data):
    '''
    Format numeric columns and add match_score metrics.

    Args:
        data: dataframe of alignment results.

    Returns:
        A formatted dataframe.
    '''
    int_cols = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end',
                'match', 'mismatch', 'indel', 'AS', 'NM']
    for c in int_cols:
        if c in data.columns:
            data[c] = pd.to_numeric(data[c], errors='coerce').fillna(0).astype(int)
    if 'tot' in data.columns and 'match' in data.columns:
        tot = pd.to_numeric(data['tot'], errors='coerce').replace(0, np.nan)
        data['match_score'] = (pd.to_numeric(data['match'], errors='coerce')/tot).fillna(0.0)
    elif 'match_score' not in data.columns:
        data['match_score'] = pd.Series(dtype=float)
    return data


def read_PAF(fname) -> pd.DataFrame:
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
    This class holds functions pertaining to sequence alignment with minibwa.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the SAM or PAF output into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        binary: minibwa binary.
        config: minibwa configuration string.
        fn_index: filename of the cached minibwa index.
    '''
    verbose = False
    binary = 'minibwa'
    config = '-t8'
    fn_index = None

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if not hasattr(self, key):
                continue
            if value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='minibwa_')

    # ------------------------------------------------------------------
    # minibwa index and alignment
    # ------------------------------------------------------------------
    def get_index(self, dbfile: str) -> str:
        '''
        This is a wrapper for minibwa to build the fn_index

        Args:
            dbfile: fasta file of sequences to index

        Returns:
            A string containing filename of the minibwa index
        '''
        minibwa = self.binary
        cmd = [minibwa, 'index', dbfile]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minibwa index failed. '+str(ret))
        self.fn_index = dbfile
        return dbfile

    def run(self, qfile, mfile: str='', dbfile: str=None, sam_file: str=None) -> pd.DataFrame:
        '''
        This is a wrapper for the minibwa aligner.

        Args:
            qfile: fasta file of query sequences to search for in database
            mfile: fasta file of the reverse pair of paired end reads
            dbfile: fasta file of the database sequences to search through
            sam_file: filename of the sam file to write and read results from

        Returns:
            a pandas dataframe of alignment results from the sam file
        '''
        work_dir = self.work_dir.name
        if sam_file is None:
            sam_file = work_dir + '/results.sam'
        config = self.config
        minibwa = self.binary

        # do paired end search if needed
        if mfile == '':
            # single end or long read alignment with minibwa outputs a PAF file
            cmd = [minibwa, 'map', '-ft8', dbfile, qfile]
            paf_file = work_dir + '/results.paf'
            with open(paf_file, 'w') as f:
                ret = subprocess.run(cmd, stdout=f)
            if ret.returncode != 0:
                raise RuntimeError('minibwa failed. '+str(ret))
            data = read_PAF(paf_file)
        else:
            # paired end alignment with minibwa outputs a SAM file
            cmd = [minibwa, 'map', '-t8', dbfile, qfile, mfile]
            with open(sam_file, 'w') as f:
                ret = subprocess.run(cmd, stdout=f)
            if ret.returncode != 0:
                raise RuntimeError('minibwa failed. '+str(ret))
            data = read_SAM(sam_file, 'r')
        return data
