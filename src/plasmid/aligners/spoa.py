#!/usr/bin/env python

# computing libraries
import numpy as np
import pandas as pd

# multi-sequence alignment libraries
import spoa

# system libraries
import tempfile

# import custom libraries
from ..genbank.terminal import color_by_array

colormaps = {'nucleic': {'A': 'green',
                         'T': 'blue',
                         'U': 'cyan',
                         'C': 'red',
                         'G': 'orange',
                         '-': 'white'}}

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
        if self.verbose:
            print({'algorithm': self.algorithm, 'genmsa': self.genmsa, 'm': self.m,
                   'n': self.n, 'g': self.g, 'e': self.e, 'q': self.q, 'c': self.c})
        self.work_dir = tempfile.TemporaryDirectory(prefix='spoa_')

    # ------------------------------------------------------------------
    # spoa multi-sequence alignment
    # ------------------------------------------------------------------
    def run(self, seq: list) -> dict:
        '''
        Runs multi-sequence alignment on provided sequences with spoa

        Args:
            seq: list of sequences

        Returns:
            A dictionary containing {'consensus':sequence,
                                     'msa':[sequences]}
        '''
        consensus, msa = spoa.poa(sequences=seq, genmsa=self.genmsa, algorithm=self.algorithm, m=self.m, n=self.n, g=self.g, e=self.e, q=self.q, c=self.c)
        out = {'consensus':consensus,
               'msa':msa,
               'params': {'algorithm':self.algorithm,
                          'genmsa':self.genmsa,
                          'm':self.m,
                          'n':self.n,
                          'g':self.g,
                          'e':self.e,
                          'q':self.q,
                          'c':self.c}}
        return out

    # ------------------------------------------------------------------
    # multi-sequence alignment formatting and visualization
    # ------------------------------------------------------------------
    @staticmethod
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

    @staticmethod
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

    @staticmethod
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
        stats = Aligner.get_msa_stats(msa)

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
