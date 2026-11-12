#!/usr/bin/env python

# computing libraries
import numpy as np

# multi-sequence alignment libraries
import spoa

# system libraries
from functools import partial
import multiprocessing
import os
import tempfile

# import custom libraries
from ..genbank.terminal import color_by_array

colormaps = {'nucleic': {'A': 'green',
                         'T': 'blue',
                         'U': 'cyan',
                         'C': 'red',
                         'G': 'orange',
                         '-': 'white'}}


def format_msa(msa: list, width: int = 100) -> str:
    '''
    Format multi-sequence alignment string to a certain width

    Args:
        msa: list of sequences
        width: characters per line

    Returns:
        A string output
    '''
    x = np.arange(0, len(msa[0]), width)
    if x[-1] < len(msa[0]):
        x = np.append(x, len(msa[0]))
    output = ''
    for i in range(len(x)-1):
        for j in range(len(msa)):
            output+= msa[j][x[i]:x[i+1]]+' '+str(x[i+1])+' seq'+str(j)+'\n'
        # find mismatches
        match = np.array([True]*(x[i+1]-x[i]))
        for j in range(len(msa)-1):
            seq1 = np.array([k for k in msa[j][x[i]:x[i+1]]])
            seq2 = np.array([k for k in msa[j+1][x[i]:x[i+1]]])
            match*= (seq1==seq2)
        mchar = np.array([' ']*(x[i+1]-x[i]))
        mchar[~match] = '.'
        output+= ''.join(mchar)+'\n'
    return output


def get_msa_stats(msa: list) -> list:
    '''
    Compute counts for each character in the multi-alignment

    Args:
        msa: list sequences

    Returns:
        list with columns [characters, frequency]
    '''
    seq = np.array([[i for i in m] for m in msa])
    seq = seq.T
    hist = [np.unique(s, return_counts=True) for s in seq]
    # sort the characters by frequency
    out = []
    for i in range(len(hist)):
        freq = hist[i][1]
        chars = hist[i][0]
        s = np.argsort(freq)
        s = s[::-1]
        out.append([chars[s], freq[s]])
    return out


def colorize_msa(msa: list, width: int = 100, cmap: dict = None) -> str:
    '''
    Formats a multi-sequence alignment with colors and to a certain width

    Args:
        msa: list of sequences
        width: characters per line
        cmap: colormap from characters to colors

    Returns:
        A string output
    '''
    # use default nucleic acid colors
    if cmap is None:
        cmap = colormaps['nucleic']

    # get character frequencies
    stats = get_msa_stats(msa)

    # colorize characters in the minor population
    seq = np.array([[i for i in m] for m in msa])
    dim = seq.shape
    carr = np.array([['']*dim[1]]*dim[0])

    for i in range(len(msa[0])):
        chars = stats[i][0]
        freq = stats[i][1]
        b = 1
        # catch cases where at least two characters are frequent
        if len(chars) > 1 and freq[0] == freq[1]:
            b = 0
        for j in range(b, len(chars)):
            c = seq[:, i]==chars[j]
            carr[c, i] = chars[j]

    # format lines
    x = np.arange(0, len(msa[0]), width)
    if x[-1] < len(msa[0]):
        x = np.append(x, len(msa[0]))

    # format to colors
    output = ''
    for i in range(len(x)-1):
        for j in range(len(msa)):
            start = x[i]
            stop = x[i+1]
            text = msa[j][start:stop]
            fgarr = carr[j, start:stop]
            text = color_by_array(text, cmap, bgarr=fgarr)
            output+= text+' '+str(x[i+1])+' seq'+str(j)+'\n'
    return output


class Aligner:
    '''
    This class holds functions pertaining to multi-sequence alignment with spoa.

    Attributes:
        verbose: print run settings when True.
        algorithm: spoa algorithm to use for multi-sequence alignment.
        genmsa: generate multiple sequence alignment.
        m: spoa parameter m.
        n: spoa parameter n.
        g: spoa parameter g.
        e: spoa parameter e.
        q: spoa parameter q.
        c: spoa parameter c.
        threads: number of clusters to align in parallel.
        batch_size: number of clusters submitted to the pool at once.
    '''
    verbose = False
    algorithm = 0
    genmsa = True
    m = 5
    n = -4
    g = -8
    e = -6
    q = -10
    c = -4
    threads = os.cpu_count()
    batch_size = 1000

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='spoa_')

    # ------------------------------------------------------------------
    # spoa multi-sequence alignment
    # ------------------------------------------------------------------
    def run(self, seq: list) -> list:
        '''
        Runs multi-sequence alignment on clusters of sequences with spoa.
        Multiple clusters are aligned in parallel with multiprocessing.

        Args:
            seq: list of lists of sequences, one list per cluster.

        Returns:
            A list of dictionaries, one per cluster, in the same order as the
            input. Each dictionary contains
            {'cluster_id':index, 'consensus':sequence, 'msa':[sequences]}.

        Raises:
            ValueError: if no clusters are given or if a cluster is empty.
        '''
        if len(seq) == 0:
            raise ValueError('at least one cluster is required')
        for cluster_id, cluster in enumerate(seq):
            if len(cluster) == 0:
                raise ValueError(f'cluster {cluster_id} is empty')

        poa = partial(spoa.poa, genmsa=self.genmsa, algorithm=self.algorithm,
                      m=self.m, n=self.n, g=self.g, e=self.e, q=self.q,
                      c=self.c)

        if self.verbose:
            print(f'spoa: aligning {len(seq)} cluster(s) with {self.threads} threads')

        if self.threads > 1 and len(seq) > 1:
            results = []
            with multiprocessing.Pool(processes=self.threads) as pool:
                for i in range(0, len(seq), self.batch_size):
                    batch = seq[i:i+self.batch_size]
                    results += pool.map(poa, batch)
        else:
            results = [poa(cluster) for cluster in seq]

        return [{'cluster_id': cluster_id,
                 'consensus': consensus,
                 'msa': msa}
                for cluster_id, (consensus, msa) in enumerate(results)]
