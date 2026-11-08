#!/usr/bin/env python
"""Class functions for multi-sequence alignment with spoa."""

# computing libraries
import numpy as np
import pandas as pd

# multi-sequence alignment libraries
import spoa

# system and time
import tempfile

__all__ = ['Aligner']


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
        work_dir: temporary working directory for intermediate files
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
    work_dir = None

    def __init__(self, args=None, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            args: Command line arguments or dictionary of settings. Defaults to None.
            kwargs: settings that override the class defaults.

        Raises:
            TypeError: if kwargs contains a name that is not an Aligner setting.
        '''
        settings = {k: v for k, v in vars(type(self)).items()
                    if not k.startswith('__') and not callable(v)}
        if args is not None:
            values = args if isinstance(args, dict) else vars(args)
            for key, value in values.items():
                if value is not None and key in settings:
                    setattr(self, key, value)
        for key, value in kwargs.items():
            if key not in settings:
                raise TypeError('unknown Aligner setting: '+key)
            setattr(self, key, value)
        if self.verbose:
            print({'algorithm': self.algorithm, 'genmsa': self.genmsa, 'm': self.m,
                   'n': self.n, 'g': self.g, 'e': self.e, 'q': self.q, 'c': self.c})
        self.create_tempdir()

    def __del__(self):
        '''
        Destructor class to clean up temporary files
        '''
        try:
            self.work_dir.cleanup()
        except Exception:
            pass

    def __enter__(self):
        '''
        Entry function for use in with statements
        '''
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        '''
        Exit function for use in with statements
        '''
        self.cleanup()
        return False

    def create_tempdir(self):
        '''
        Create a temporary working directory for the aligner. Name of the temporary working directory is stored in self.work_dir
        '''
        self.work_dir = tempfile.TemporaryDirectory(prefix='spoa_')

    def cleanup(self):
        '''
        Remove the temporary working directory.
        '''
        if self.work_dir is not None:
            self.work_dir.cleanup()
            self.work_dir = None

    def spoa(self, seq: list) -> dict:
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
