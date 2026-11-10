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
    This class holds functions pertaining to sequence alignment with DIAMOND blast.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the DIAMOND output into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        binary: DIAMOND binary.
        config: DIAMOND configuration string.
        fn_db: filename of the cached DIAMOND database.
    '''
    verbose = False
    binary = 'diamond'
    config = 'blastn'
    fn_db = None

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
        self.work_dir = tempfile.TemporaryDirectory(prefix='diamond_')

    # ------------------------------------------------------------------
    # diamond database and alignment
    # ------------------------------------------------------------------
    def get_index(self, database: str, fn_db: str=None) -> str:
        '''
        Build a DIAMOND database.

        Args:
            database: filename containing the reference sequences.
            fn_db: filename of the database file inside the working directory.

        Returns:
            Filename of the DIAMOND database.
        '''
        if fn_db is None:
            fn_db = os.path.join(self.work_dir.name, "db")
        cmd = [self.binary, 'makedb', '--in', database, '-d', fn_db]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('diamond makedb failed. '+str(ret))
        self.fn_db = fn_db
        return fn_db

    def run(self, qfile, dbfile: str=None, out_file: str=None) -> pd.DataFrame:
        '''
        Wrapper for performing sequence alignment with DIAMOND blast.

        Args:
            qfile: filename of the query sequences (fasta or fastq).
            dbfile: filename of the DIAMOND database. Required if no database
                has been built with get_index.
            out_file: filename of the DIAMOND output file.

        Returns:
            A pandas dataframe of alignment results.

        Raises:
            ValueError: if no database has been built and no database is given.
            RuntimeError: if the DIAMOND blast command fails.
        '''
        work_dir = self.work_dir.name
        if dbfile is None and self.fn_db is not None:
            dbfile = self.fn_db
        else:
            raise ValueError('database is required when no diamond database has been built')

        if out_file is None:
            out_file = os.path.join(self.work_dir.name, "results.out")

        cols = ['qseqid', 'qlen', 'qstart', 'qend', 'qframe',
                'sseqid', 'slen', 'sstart', 'send',
                'length', 'pident', 'nident', 'mismatch', 'score', 'evalue']
        cmd = [self.binary] + self.config.split(' ') + ['-d', dbfile, '-q', qfile, '-o', out_file, '-f', '6'] + cols
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('diamond blast failed. '+str(ret))

        # parse the alignment output
        data = self.parse_diamond(out_file, ' '.join(cols))

        # format output columns
        int_cols = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end',
                    'match', 'tot', 'NM', 'AS']
        for c in int_cols:
            if c in data.columns:
                data[c] = data[c].astype(int)
        return data

    def parse_diamond(self, fname: str, cols: str=None) -> pd.DataFrame:
        '''
        Parse a DIAMOND tabular output file into a pandas dataframe.

        Args:
            fname: filename of the DIAMOND output file.
            cols: whitespace delimited column names of the DIAMOND output format.

        Returns:
            A pandas dataframe of alignment results.

        Notes:
            DIAMOND tabular output columns map onto the alignment statistics
            columns used by the other aligners, such as ``qseqid``, ``qlen``,
            ``qstart``, ``qend``, ``sseqid``, ``slen``, ``sstart``, and ``send``.
        '''
        empty_cols = ["query_id", "q_len", "q_start", "q_end", "orientation",
                      "database_id", "t_len", "t_start", "t_end",
                      "match", "tot", "NM", "AS", "match_score"]
        if cols is None:
            cols = ('qseqid qlen qstart qend qframe '
                    'sseqid slen sstart send '
                    'length pident nident mismatch score evalue')
        if os.path.getsize(fname) == 0:
            return pd.DataFrame(columns=empty_cols)
        df = pd.read_csv(fname, delimiter="\t", header=None)
        if len(df) == 0:
            return pd.DataFrame(columns=empty_cols)
        df.columns = cols.split(' ')
        key = {"qseqid": "query_id", "qlen": "q_len",
               "qstart": "q_start", "qend": "q_end",
               "sseqid": "database_id", "slen": "t_len",
               "sstart": "t_start", "send": "t_end",
               "nident": "match", "mismatch": "NM", "score": "AS"}
        df = df.rename(columns=key)

        # get the orientation
        c = df["qframe"] > 0
        df.loc[c, "orientation"] = "+"
        df.loc[~c, "orientation"] = "-"

        # fix the coordinates of alignments on the reverse strand
        c = df["orientation"] == "-"
        q_start = df.loc[c, "q_start"]
        q_end = df.loc[c, "q_end"]
        df.loc[c, "q_start"] = q_end
        df.loc[c, "q_end"] = q_start

        # correct the start position by 1bp
        df["q_start"] -= 1
        df["t_start"] -= 1

        # compute alignment statistics
        df['tot'] = df['q_end'] - df['q_start']
        df['match'] = df['tot'] - df['NM']
        tot = pd.to_numeric(df['tot'], errors='coerce').replace(0, np.nan)
        df['match_score'] = (pd.to_numeric(df['match'], errors='coerce')/tot).fillna(0.0)

        return df
