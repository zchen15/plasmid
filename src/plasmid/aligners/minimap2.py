#!/usr/bin/env python

# computing libraries
import numpy as np
import pandas as pd

# BAM file parsing
import pysam

# system libraries
import os
import subprocess
import tempfile


def get_BAM_tags(tag, cols) -> list:
    '''
    Parse BAM tags for columns.

    Args:
        tag: list of tuples of tag key and value.
        cols: values to extract.

    Returns:
        list of values.
    '''
    out = [None]*len(cols)
    for j in range(len(cols)):
        for k in tag:
            if k[0] == cols[j]:
                out[j] = k[1]
                break
    return out


def read_SAM(fname: str, mode: str='r') -> pd.DataFrame:
    '''
    Read a SAM or BAM file.

    Args:
        fname: filename of the SAM or BAM file.
        mode: 'r' for SAM or 'rb' for BAM.

    Returns:
        A pandas dataframe of alignment information.
    '''
    samfile = pysam.AlignmentFile(str(fname), mode=mode, check_sq=False)

    out = []
    col2 = ['NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 's2', 'de', 'rl']
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
    col1 = ['q_id', 'q_len', 'q_start', 'q_end', 'orientation',
            't_id', 't_len', 't_start', 't_end', 'CIGAR', 'mapq']
    cols = col1 + col2
    df = pd.DataFrame(out, columns=cols)

    df['orientation'] = df['orientation'].replace({True: '+', False: '-'})

    # compute blast like score
    df['tot'] = df['q_end'] - df['q_start']
    df['match'] = df['tot'] - df['NM']
    tot = pd.to_numeric(df['tot'], errors='coerce').replace(0, np.nan)
    df['match_score'] = (pd.to_numeric(df['match'], errors='coerce')/tot).fillna(0.0)
    return df


def read_PAF(fname: str) -> pd.DataFrame:
    '''
    Read a minimap2 PAF file.

    Args:
        fname: PAF filename.

    Returns:
        A pandas dataframe of alignment results.
    '''
    col1 = ['q_id', 'q_len', 'q_start', 'q_end', 'orientation',
            't_id', 't_len', 't_start', 't_end', 'match', 'tot', 'mapq']
    col2 = ['tp', 'cm', 's1', 's2', 'NM', 'AS', 'ms', 'nn', 'rl', 'CIGAR']
    key = {name: i for i, name in enumerate(col2)}

    out = []
    with open(str(fname), 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line[0] == '#':
                continue
            fields = line.split('\t')
            if len(fields) < len(col1):
                continue
            info = ['']*len(col2)
            for field in fields[len(col1):]:
                x = field.split(':', 2)
                if len(x) != 3:
                    continue
                tag, value = x[0], x[2]
                if tag == 'cg':
                    tag = 'CIGAR'
                if tag in key:
                    info[key[tag]] = value
            out.append(fields[:len(col1)]+info)

    data = pd.DataFrame(out, columns=col1+col2)

    # format numeric columns
    int_cols = ['q_len', 'q_start', 'q_end', 't_len', 't_start', 't_end',
                'match', 'tot', 'mapq', 'cm', 's1', 's2', 'NM', 'AS',
                'ms', 'nn', 'rl']
    for c in int_cols:
        data[c] = pd.to_numeric(data[c], errors='coerce').fillna(0).astype(int)

    # compute blast like score
    tot = data['tot'].replace(0, np.nan)
    data['match_score'] = (data['match']/tot).fillna(0.0)
    return data


class Aligner:
    '''
    This class holds functions pertaining to sequence alignment with minimap2.

    The class takes fasta and fastq filenames as input. Alignment results are
    parsed from the PAF, SAM, or BAM output into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        binary: minimap2 binary.
        config: minimap2 configuration string for alignment.
        index_config: minimap2 configuration string for index building.
    '''
    verbose = False
    binary = 'minimap2'
    config = '-c -k 15 -w 7'
    index_config = '-k 15 -w 7'

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='minimap2_')

    # ------------------------------------------------------------------
    # minimap2 index and alignment
    # ------------------------------------------------------------------
    def get_index(self, database: str, fn_index: str=None):
        '''
        Build a minimap2 index.

        Args:
            database: filename containing the reference sequences.
            fn_index: name of the index file inside the working directory.

        Returns:
            Filename of the minimap2 index.
        '''
        if fn_index is None:
            fn_index = os.path.join(self.work_dir.name, "index.mmi")

        options = str(self.index_config).split(' ') if str(self.index_config) else []
        cmd = [self.binary] + options + ['-d', fn_index, database]

        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minimap2 index failed. '+str(ret))

        return fn_index

    def run(self, qfile: str, mfile: str=None, dbfile: str=None,
            paf_file: str=None, sam_file: str=None) -> pd.DataFrame:
        '''
        Wrapper for performing sequence alignment with minimap2.

        Args:
            qfile: filename of the query sequences (fasta or fastq).
            mfile: filename of the paired end query sequences (fasta or fastq).
            dbfile: filename of the minimap2 index of the database. Required if no index
                has been built with get_index.
            paf_file: filename of the PAF alignment file.
            sam_file: filename of the SAM alignment file.

        Returns:
            A pandas dataframe of alignment results from the PAF or SAM file.

        Raises:
            ValueError: if no index has been built and no database is given.
        '''
        if dbfile is None:
            raise ValueError('database is required')

        config = str(self.config).split(' ')
        if sam_file is not None:
            config = ['-a'] + config
            outfile = sam_file
        elif paf_file is not None:
            outfile = paf_file
        else:
            outfile = os.path.join(self.work_dir.name, 'results.paf')

        # call the aligner
        if mfile is None:
            cmd = [self.binary]+config+[dbfile, qfile, '-o', outfile]
        else:
            cmd = [self.binary]+config+[dbfile, qfile, mfile, '-o', outfile]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('minimap2 failed. '+str(ret))

        # parse the alignment output
        if sam_file is not None:
            data = read_SAM(outfile)
        else:
            data = read_PAF(outfile)

        return data
