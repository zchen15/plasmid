#!/usr/bin/env python

# computing libraries
import numpy as np
import pandas as pd

# system libraries
import multiprocessing
import tempfile

from . import io

import parasail


COLUMNS = ['query_id', 'ref_id', 'strand', 'AS',
           'q_len', 'q_start', 'q_end',
           't_len', 't_start', 't_end',
           'match', 'mismatch', 'indel', 'CIGAR']


def cigar_to_stats(cig: str, qseq: str, rseq: str) -> list:
    '''
    Get alignment statistics associated with parasail output.

    Args:
        cig: cigar string
        qseq: query sequence
        rseq: reference sequence

    Returns:
        list with columns [q_len, q_start, q_end, t_len, t_start, t_end, M, Xm, indel]
    '''
    cig = io.cigar_decode(cig)

    # compute match boundaries
    q_start = 0
    q_end = len(qseq)
    q_len = len(qseq)
    t_start = 0
    t_end = len(rseq)
    t_len = len(rseq)
    if len(cig) > 0 and cig[0][1] == 'I':
        t_start = cig[0][0]
    elif len(cig) > 0 and cig[0][1] == 'D':
        q_start = cig[0][0]
    if len(cig) > 0 and cig[-1][1] == 'I':
        t_end -= cig[-1][0]
    elif len(cig) > 0 and cig[-1][1] == 'D':
        q_end -= cig[-1][0]

    x = np.array([cig[i][0] for i in range(len(cig))])
    y = np.array([cig[i][1] for i in range(len(cig))])

    # get number of matches
    M = np.sum(x[y == '='])
    M += np.sum(x[y == 'M'])

    # get number of mismatches
    Xm = np.sum(x[y == 'X'])

    # get indel count
    indel = np.sum(x[y == 'I'])
    indel += np.sum(x[y == 'D'])
    return [q_len, q_start, q_end, t_len, t_start, t_end, M, Xm, indel]


def _matrix(mat):
    '''
    Resolve a parasail substitution matrix.

    Args:
        mat: parasail matrix object or matrix name.

    Returns:
        A parasail matrix object.
    '''
    if isinstance(mat, str):
        return getattr(parasail, mat)
    return mat


def sg_trace_scan_16(task: list) -> list:
    '''
    Align one query against one reference in the forward orientation with
    parasail.sg_trace_scan_16. Used by the multiprocessing pool.

    Args:
        task: list of [qname, qseq, rname, rseq, g1, g2, mat].

    Returns:
        list containing one alignment result row.
    '''
    qname, qseq, rname, rseq, g1, g2, mat = task
    mat = _matrix(mat)

    res = parasail.sg_trace_scan_16(qseq, rseq, g1, g2, mat)
    cigar = res.cigar.decode.decode('UTF-8')

    stats = cigar_to_stats(cigar, qseq, rseq)
    row = [qname, rname, 1, res.score]+stats+[cigar]
    return [row]


def _empty_results() -> pd.DataFrame:
    '''
    Create an empty, typed alignment dataframe.

    Returns:
        An empty pandas dataframe with the parasail alignment schema.
    '''
    int_cols = ['strand', 'AS', 'q_len', 'q_start', 'q_end',
                't_len', 't_start', 't_end', 'match', 'mismatch',
                'indel', 'NM']
    data = {}
    for c in int_cols:
        data[c] = pd.Series(dtype='int64')
    for c in ['query_id', 'ref_id', 'CIGAR']:
        data[c] = pd.Series(dtype='object')
    data['match_score'] = pd.Series(dtype='float64')
    return pd.DataFrame(data)[COLUMNS+['NM', 'match_score']]


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with parasail.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        g1: parasail gap open score.
        g2: parasail gap extension score.
        mat: parasail substitution matrix.
        threads: number of threads to run the alignment task.
        batch_size: number of query sequences to align per reference sequence.
        chunksize: number of tasks sent to each worker at a time.
    '''
    verbose = False
    g1 = 5
    g2 = 1
    mat = parasail.nuc44
    threads = 4
    batch_size = 1000
    chunksize = 100

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='parasail_')

    # ------------------------------------------------------------------
    # parasail alignment
    # ------------------------------------------------------------------
    def run(self, qfile: str, dbfile: str) -> pd.DataFrame:
        '''
        Wrapper for performing sequence alignment with parasail.

        Args:
            qfile: filename of the query sequences (fasta or fastq).
            dbfile: filename of the database sequences to align against.

        Returns:
            A pandas dataframe of alignment results.

        Raises:
            ValueError: if dbfile is not provided.
        '''
        if dbfile is None:
            raise ValueError('database is required')

        # read in the reference sequences
        try:
            ref = io.read_fasta(dbfile)
        except ValueError:
            return _empty_results()
        # read in the query sequences
        try:
            qry = io.read_fasta(qfile)
        except ValueError:
            return _empty_results()

        # resolve the matrix by name so it can be sent to worker processes
        mat = self.mat
        if not isinstance(mat, str):
            name = getattr(mat, 'name', None)
            if isinstance(name, bytes):
                name = name.decode('utf-8')
            if name:
                mat = name

        queries = [[qname, qseq]
                   for qname, qseq in qry[['name', 'sequence']].values]

        out = []
        with multiprocessing.Pool(processes=self.threads) as pool:
            for rname, rseq in ref[['name', 'sequence']].values:
                # Batch the query sequences ~1000 per 1 reference sequence
                n_batches = (len(queries) + self.batch_size - 1)//self.batch_size
                for i in range(0, len(queries), self.batch_size):
                    batch = queries[i:i+self.batch_size]
                    tasks = [[qname, qseq, rname, rseq,
                              self.g1, self.g2, mat]
                             for qname, qseq in batch]
                    if self.verbose:
                        batch_id = i//self.batch_size + 1
                        print(f'parasail: aligning batch {batch_id}/{n_batches} '
                              f'({len(batch)} queries against {rname})')
                    results = pool.map(sg_trace_scan_16, tasks,
                                       chunksize=self.chunksize)
                    for r in results:
                        out += r

        if len(out) == 0:
            return _empty_results()

        df = pd.DataFrame(out, columns=COLUMNS)
        df['NM'] = df['mismatch']+df['indel']

        # format numeric columns and add match_score metrics
        int_cols = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end',
                    'match', 'mismatch', 'indel', 'AS', 'NM']
        for c in int_cols:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0).astype(int)
        df['match_score'] = (
            pd.to_numeric(df['match'], errors='coerce')/
            pd.to_numeric(df['q_len'], errors='coerce')
        ).fillna(0.0)
        return df
