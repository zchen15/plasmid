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

import parasail


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
    if 'match_score' not in data.columns:
        if 'match' in data.columns and 'q_len' in data.columns:
            match_score = pd.to_numeric(data['match'], errors='coerce')/pd.to_numeric(data['q_len'], errors='coerce')
            data['match_score'] = match_score.fillna(0.0)
    return data


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with parasail.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        config: parasail configuration string.
        fn_index: filename of the cached parasail database.
        binary: parasail binary.
    '''
    verbose = False
    binary = 'parasail'
    config = 'sg_scan_16 -g1=1 -g2=1 -mat=1'
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
        self.work_dir = tempfile.TemporaryDirectory(prefix='parasail_')

    # ------------------------------------------------------------------
    # parasail alignment
    # ------------------------------------------------------------------
    def get_index(self, dbfile: str) -> str:
        '''
        Build a parasail database.

        Args:
            dbfile: filename containing the reference sequences.

        Returns:
            Filename of the parasail database.
        '''
        self.fn_index = dbfile
        return dbfile

    def run(self, qfile, mfile=None, dbfile=None, paf_file=None, sam_file=None) -> pd.DataFrame:
        '''
        Wrapper for performing sequence alignment with parasail.

        Args:
            qfile: filename of the query sequences (fasta or fastq).
            mfile: filename of the paired end query sequences (fasta or fastq).
            dbfile: filename of the database sequences to align against.

        Returns:
            A pandas dataframe of alignment results.

        Raises:
            ValueError: if dbfile is not provided.
        '''
        work_dir = self.work_dir.name
        if dbfile is None:
            raise ValueError('database is required')

        # read in the reference sequences
        ref = io.read_fasta(dbfile)
        # read in the query sequences
        qry = io.read_fasta(qfile)

        out = []
        c = str(self.config).split(' ')
        method = c[0]
        g1 = 1
        g2 = 1
        mat = 1
        for p in c[1:]:
            k, v = str(p).split('=')
            k = k.lstrip('-')
            if k == 'g1':
                g1 = int(v)
            elif k == 'g2':
                g2 = int(v)
            elif k == 'mat':
                mat = int(v)

        # align the queries against the references and their reverse complements
        for rname, rseq in ref[['name', 'sequence']].values:
            for qname, qseq in qry[['name', 'sequence']].values:
                qseqR = io.revcomp(qseq)
                res = getattr(parasail, method)(qseq, rseq, g1, g2, mat)
                resR = getattr(parasail, method)(qseqR, rseq, g1, g2, mat)
                cigarF = ''
                cigarR = ''
                if 'trace' in str(method):
                    cigarF = res.cigar.decode.decode('UTF-8')
                    cigarR = resR.cigar.decode.decode('UTF-8')
                sF = io.cigar_to_stats(cigarF, qseq, rseq)
                sF = [qname, rname, 1, res.score]+sF+[cigarF]
                sR = io.cigar_to_stats(cigarR, qseq, rseq)
                sR = [qname, rname, -1, resR.score]+sR+[cigarR]
                out.append(sF)
                out.append(sR)
        col = ['query_id','ref_id','strand','AS',
               'q_len','q_start','q_end',
               't_len','t_start','t_end',
               'match','mismatch','indel','CIGAR']
        df = pd.DataFrame(out, columns=col)
        df['NM'] = df['mismatch']+df['indel']
        return format_alignment(df)
