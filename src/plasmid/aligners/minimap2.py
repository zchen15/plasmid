#!/usr/bin/env python
"""Class functions for sequence alignment with minimap2.

This module consolidates the minimap2 functionality that was spread across
``aligners/aligner.py`` and ``tbd/search.py``. The wrapper keeps the
``Aligner`` API used by the rest of the package (``minimap2`` and
``minimap2_get_index``) and the same output columns, while using the file
based interface and subprocess error handling from ``search.py``.

Alignment settings are class attributes. They can be overridden per instance
through keyword arguments or an argparse namespace. Inputs are fasta or fastq
filenames. The module only depends on the core ``numpy``/``pandas`` stack;
``pysam`` is imported lazily and is only needed to read BAM output.
"""

# computing libraries
import numpy as np
import pandas as pd

# system libraries
import os
import re
import subprocess
import tempfile

from . import io

__all__ = ['Aligner']


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with minimap2.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the PAF, SAM, or BAM output into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        work_dir: TemporaryDirectory used for intermediate files.
        k: minimap2 k-mer size.
        w: minimap2 minimizer window size.
        minimap2_config: minimap2 configuration string.
        cigar: add -c to PAF output when True.
        index: filename of the cached minimap2 index.
        threads: number of minimap2 threads.
        minimap2_bin: minimap2 binary.
    '''
    verbose = False
    work_dir = None
    k = 7
    w = 1
    minimap2_config = '-x map-ont -P'
    cigar = True
    index = None
    threads = None
    minimap2_bin = 'minimap2'

    # alignment file parsing is implemented in aligners.io
    read_PAF = staticmethod(io.read_PAF)
    read_paf = staticmethod(io.read_PAF)
    read_SAM = staticmethod(io.read_SAM)
    read_BAM = staticmethod(io.read_BAM)

    def __init__(self, args=None, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            args: argparse namespace or dictionary of settings. Defaults to None.
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
            print({'k': self.k, 'w': self.w, 'config': self.minimap2_config,
                   'threads': self.threads})
        self.work_dir = tempfile.TemporaryDirectory(prefix='minimap2_')

    def __del__(self):
        '''
        Destructor class to clean up temporary files.
        '''
        try:
            self.cleanup()
        except Exception:
            pass

    def __enter__(self):
        '''
        Entry function for use in with statements.
        '''
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        '''
        Exit function for use in with statements.
        '''
        self.cleanup()
        return False

    def cleanup(self):
        '''
        Remove the temporary working directory.
        '''
        if self.work_dir is not None:
            self.work_dir.cleanup()
            self.work_dir = None

    # ------------------------------------------------------------------
    # minimap2 index and alignment
    # ------------------------------------------------------------------
    def minimap2_get_index(self, database, filename='index.mmi'):
        '''
        Build a minimap2 index.

        Args:
            database: filename containing the reference sequences.
            filename: name of the index file inside the working directory.

        Returns:
            Filename of the minimap2 index.
        '''
        work_dir = self.work_dir.name
        fn_index = os.path.join(work_dir, filename)
        cmd = [self.minimap2_bin,
               '-k', str(self.k),
               '-w', str(self.w)]
        if self.threads is not None:
            cmd+= ['-t', str(self.threads)]
        cmd+= ['-d', fn_index, str(database)]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minimap2 index failed. '+str(ret))

        self.index = fn_index
        return fn_index

    def minimap2(self, query, query2=None, database=None, paf_file=None,
                 sam_file=None, config=None):
        '''
        Wrapper for performing sequence alignment with minimap2.

        Args:
            query: fasta or fastq filename containing the query sequences.
            query2: fasta or fastq filename containing the reverse pair of
                paired end reads. Defaults to None.
            database: fasta or fastq filename, a minimap2 index filename, or
                None to use the index built by minimap2_get_index.
            paf_file: optional output PAF filename. Defaults to a file in the
                working directory.
            sam_file: optional output SAM filename. Takes precedence over
                paf_file.
            config: optional minimap2 configuration string. Overrides
                minimap2_config.

        Returns:
            A pandas dataframe of alignment results.
        '''
        work_dir = self.work_dir.name
        qfile = str(query)
        mfile = None if query2 is None else str(query2)

        if database is None:
            if self.index is None:
                raise ValueError('database is required when no index has been built')
            dbfile = str(self.index)
        else:
            dbfile = str(database)

        # assemble the minimap2 configuration
        if config is None:
            config = self.minimap2_config
        if isinstance(config, str):
            tokens = config.split(' ')
        else:
            tokens = list(config)

        # the -a flag switches minimap2 to SAM output
        sam_output = sam_file is not None or '-a' in tokens
        if sam_output:
            tokens = [i for i in tokens if i != '-c']
            if '-a' not in tokens:
                tokens = ['-a']+tokens
        elif self.cigar and '-c' not in tokens:
            tokens+= ['-c']

        if '-k' not in tokens:
            tokens+= ['-k', str(self.k)]
        if '-w' not in tokens:
            tokens+= ['-w', str(self.w)]
        if self.threads is not None and '-t' not in tokens:
            tokens+= ['-t', str(self.threads)]

        # define the output file
        if sam_file is not None:
            outfile = str(sam_file)
        elif paf_file is not None:
            outfile = str(paf_file)
        elif sam_output:
            outfile = work_dir+'/results.sam'
        else:
            outfile = work_dir+'/results.paf'

        # call the aligner
        if mfile is None:
            cmd = [self.minimap2_bin]+tokens+[dbfile, qfile, '-o', outfile]
        else:
            cmd = [self.minimap2_bin]+tokens+[dbfile, qfile, mfile, '-o', outfile]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minimap2 failed. '+str(ret))

        # parse the alignment output
        if sam_output:
            data = Aligner.read_SAM(outfile)
        else:
            data = Aligner.read_PAF(outfile)
        return Aligner.format_alignment(data)

    # ------------------------------------------------------------------
    # alignment file parsing
    # ------------------------------------------------------------------
    @staticmethod
    def cigar_query_bounds(cigar, q_len):
        '''
        Compute the aligned query start and end from a CIGAR string.

        Args:
            cigar: CIGAR string.
            q_len: query length.

        Returns:
            list containing [q_start, q_end].
        '''
        ops = re.findall('([0-9]+)([MIDNSHP=X])', cigar)
        q_start = 0
        q_end = q_len
        if len(ops) > 0 and ops[0][1] == 'S':
            q_start = int(ops[0][0])
        if len(ops) > 0 and ops[-1][1] == 'S':
            q_end = q_len-int(ops[-1][0])
        return [q_start, q_end]

    @staticmethod
    def cigar_ref_span(cigar):
        '''
        Compute the reference span of a CIGAR string.

        Args:
            cigar: CIGAR string.

        Returns:
            Number of reference bases spanned by the alignment.
        '''
        ops = re.findall('([0-9]+)([MIDNSHP=X])', cigar)
        span = 0
        for n, op in ops:
            if op in 'MDN=X':
                span+= int(n)
        return span

    # ------------------------------------------------------------------
    # dataframe utilities
    # ------------------------------------------------------------------
    @staticmethod
    def format_alignment(data):
        '''
        Format numeric columns and add match_score and similarity metrics.

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

        # add similarity metric
        if all(i in data.columns for i in ['q_len', 't_len', 'NM']):
            if len(data) == 0:
                data['similarity'] = pd.Series(dtype=float)
            else:
                q_len = pd.to_numeric(data['q_len'], errors='coerce').fillna(0).values
                t_len = pd.to_numeric(data['t_len'], errors='coerce').fillna(0).values
                nm = pd.to_numeric(data['NM'], errors='coerce').fillna(0).values
                v1 = np.min([q_len, t_len], axis=0)
                v1 = np.where(v1 == 0, 1, v1)
                data['similarity'] = 1-nm/v1
                data.loc[data['similarity'] < 0, 'similarity'] = 0
        elif 'similarity' not in data.columns:
            data['similarity'] = pd.Series(dtype=float)
        return data
