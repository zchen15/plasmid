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


def trim_spaces(text) -> str:
    '''
    Leading and previous spaces

    Args:
        text: text string

    Returns:
        string with spaces left and right removed
    '''
    while len(text) > 0 and text[0] == ' ':
        text = text[1:]
    while len(text) > 0 and text[-1] == ' ':
        text = text[:-1]
    return text


def get_infernal_line(row: str, width) -> list:
    '''
    Get each line of infernal results file

    Args:
        row: each row
        width: list of number representing column width

    Returns:
        list of infernal output for each column
    '''
    s1 = 0
    out = []
    for w in width[:-1]:
        val = row[s1:s1+w]
        val = trim_spaces(val)
        out.append(val)
        s1 += w+1
    val = row[s1:]
    val = trim_spaces(val)
    out.append(val)
    return out


def read_infernal(fname: str) -> pd.DataFrame:
    '''
    Function for parsing infernal output results

    Args:
        fname: filename of infernal alignment file

    Returns:
        pandas dataframe of alignments
    '''
    with open(fname, 'r') as f:
        text = f.read()
    text = text.split('\n')
    cols = text[0]
    widths = text[1]
    text = text[2:]

    # compute column widths
    widths = widths.split(' ')
    widths = [len(i) for i in widths]

    # format column names
    cols = get_infernal_line(cols, widths)
    cols[2] = 't_accession'
    cols[4] = 'q_accession'

    # parse lines
    out = []
    for line in text:
        if len(line) > 0 and line[0] != '#':
            val = get_infernal_line(line, widths)
            out.append(val)
    df = pd.DataFrame(out, columns=cols)

    # rename columns
    cols = {"query name": "q_id",
            "target name": "t_id", "E-value": "evalue",
            "seq to": "q_end", "seq from": "q_start", "strand": "orientation",
            "mdl from": "t_start", "mdl to": "t_end",
            "description of target": "description",
            "score": "AS"}
    df = df.rename(columns=cols)
    df = df.drop(columns=['#idx'])

    # fix t_start and t_end
    c = df["orientation"] == "-"
    q_start = df.loc[c, "q_start"]
    q_end = df.loc[c, "q_end"]
    df.loc[c, "q_start"] = q_end
    df.loc[c, "q_end"] = q_start

    # format data types
    cols = ["t_start","t_end","q_start","q_end"]
    for c in cols:
        df[c] = df[c].astype(int)
    df["evalue"] = df["evalue"].astype(float)
    df["AS"] = df["AS"].astype(float)
    df["match_score"] = 1 - df["evalue"]

    # correct the start position by 1bp
    df["q_start"] -= 1
    df["t_start"] -= 1
    return df


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with infernal.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the infernal results file into a pandas dataframe.

    Attributes:
        verbose: print run settings when True.
        binary: infernal binary.
        config: infernal configuration string.
        fn_index: filename of the cached infernal database.
    '''
    verbose = False
    binary = 'infernal'
    config = '--rfam --noali --nohmmonly'
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
        self.work_dir = tempfile.TemporaryDirectory(prefix='infernal_')

    def run(self, qfile, dbfile) -> pd.DataFrame:
        '''
        Wrapper for running sequence alignment with infernal.

        Args:
            qfile: input fasta file of sequences to search against the infernal database.
            dbfile: filename of the infernal database file.

        Returns:
            A pandas dataframe of alignment results from the infernal results file.

        Raises:
            RuntimeError: if the infernal command fails.
        '''
        work_dir = self.work_dir.name

        # run infernal cmscan
        cmd = [self.binary]
        cmd += ['--fmt', '2']
        cmd += str(self.config).split(' ')
        cmd += ['--tblout', work_dir + '/results.txt']
        cmd += [dbfile, qfile]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('infernal failed. '+str(ret))

        # parse results file
        data = read_infernal(work_dir + '/results.txt')
        return data
