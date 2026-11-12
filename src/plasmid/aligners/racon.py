#!/usr/bin/env python

# computing libraries
import pandas as pd

# system libraries
import subprocess
import tempfile


def parse_fasta(text: str) -> pd.DataFrame:
    '''
    Parse fasta formatted text into a dataframe.

    Args:
        text: fasta formatted text, such as racon output from stdout.

    Returns:
        A dataframe of sequences with columns [name, sequence].
    '''
    data = []
    name = None
    seq = []
    for line in text.splitlines():
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if name is not None:
                data.append([name, ''.join(seq)])
            name = line[1:].split(' ')[0]
            seq = []
        elif name is not None:
            seq.append(line)
    if name is not None:
        data.append([name, ''.join(seq)])
    return pd.DataFrame(data, columns=['name', 'sequence'])


class Aligner:
    '''
    This class holds functions pertaining to polishing reads with raco.

    The class takes fasta and fastq filenames as input. Polished sequences are
    parsed from the polished fasta file into a pandas dataframe.

    Attributes:
        verbose: print run settings when True.
        binary: raco binary.
        config: raco configuration string.
    '''
    verbose = False
    binary = 'racon'
    config = '-w 500 -q 10'

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='racon_')

    # ------------------------------------------------------------------
    # raco polishing
    # ------------------------------------------------------------------
    def run(self, seq, paf, ref, config: str=None, out_file: str=None) -> pd.DataFrame:
        '''
        Wrapper for polishing reads against a reference alignment with raco.

        Args:
            seq: sequencing reads fasta or fastq file.
            paf: alignment file for reads against the draft reference sequence.
            ref: reference sequence fasta or fastq file.
            config: config options to pass to raco. Defaults to the raco config.
            out_file: filename of the polished fasta file.

        Returns:
            A pandas dataframe of polished sequences with columns [name, sequence].

        Raises:
            RuntimeError: if the raco command fails.
        '''
        if config is None:
            config = self.config

        # run raco
        options = str(config).split(' ') if str(config) else []
        cmd = [self.binary] + options + [seq, paf, ref]
        ret = subprocess.run(cmd, stdout=subprocess.PIPE)
        if ret.returncode != 0:
            raise RuntimeError('racon failed. '+str(ret))

        # parse the polished fasta output directly from stdout
        text = ret.stdout.decode('utf-8')
        if out_file is not None:
            with open(out_file, 'w') as f:
                f.write(text)
        data = parse_fasta(text)
        return data
