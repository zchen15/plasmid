#!/usr/bin/env python
"""Functions for reading and writing sequence and alignment files."""

# computing libraries
import numpy as np
import pandas as pd

# compression libraries
import bz2
import gzip
import Bio
import Bio.SeqIO

# BAM file parsing
import pysam

# system and time
import os
import re
import sys
import logging


def decompress(fname):
    '''
    Opens a compressed file
    fname = name of compressed file with .bz2 or .gz
    returns filehandle passed through the decompression
    '''
    fname = str(fname)
    if fname.endswith('.bz2'):
        return bz2.open(fname, mode='rt')
    elif fname.endswith('.gz'):
        return gzip.open(fname, mode='rt')
    else:
        return open(fname, 'r')


def compress(fname):
    '''
    Checks if a filename will be bz2 or gz compressed.
    fname = filename with .bz2 or .gz or no extension
    return a file handle
    '''
    fname = str(fname)
    if fname.endswith('.bz2'):
        return bz2.open(fname, mode='wt')
    elif fname.endswith('.gz'):
        return gzip.open(fname, mode='wt')
    else:
        return open(fname, 'w')


def add_seq(df):
    '''
    Adds sequence and quality information back to dataframe
    df = dataframe containing at least the columns:
        [filename, id, seek1, rlen] for fasta files
        [filename, id, seek1, seek2, rlen] for fastq files
    return dataframe with columns [id, sequence, quality] or [id, sequence]
    to do testing
    '''
    # if data is already there, dont do anything
    if 'sequence' in df.keys():
        return df
    # check if relevant columns are there
    cols = ['filename', 'name', 'seek1', 'rlen']
    for c in cols:
        if (c in df.keys()) == False:
            logging.warning('columns '+c+' not found in dataframe')
            return df
    # sort and iterate
    df = df.sort_values(by=['filename', 'seek1']).copy()
    data = []
    if 'seek2' in df.keys():
        # work on fastq
        x = df[['filename', 'name', 'seek1', 'seek2', 'rlen']].values
        files = np.unique(x[:, 0])
        for fname in files:
            y = x[x[:, 0] == fname]
            # decompress file if needed
            with decompress(fname) as f:
                for i in range(0, len(y)):
                    f.seek(y[i, 2])
                    seq = f.read(y[i, 4])
                    f.seek(y[i, 3])
                    q = f.read(y[i, 4])
                    data.append([seq, q])
        df['sequence'] = [i[0] for i in data]
        df['quality'] = [i[1] for i in data]
        cols = ['seek1', 'seek2', 'rlen', 'filename']
        df = df.drop(columns=cols)
    else:
        # work on fasta
        x = df[['filename', 'name', 'seek1', 'rlen']].values
        files = np.unique(x[:, 0])
        for fname in files:
            y = x[x[:, 0] == fname]
            # decompress the file
            with decompress(fname) as f:
                for i in range(0, len(y)):
                    f.seek(y[i, 2])
                    seq = f.read(y[i, 3])
                    seq = seq.replace('\n', '')
                    data.append(seq)
        df['sequence'] = data
        cols = ['seek1', 'rlen', 'filename']
        df = df.drop(columns=cols)
    return df


def write_fastq(fname, data):
    '''
    Write list of sequences into a fastq file format
    data = dataframe with columns [name, sequence, quality]
    fname = file to write the information.
            File extension can be .bz2 or .gz if compression is desired
    '''
    data = data[['name', 'sequence', 'quality']].values
    with compress(fname) as f:
        for i in data:
            text = '@' + str(i[0]) + '\n'
            text += str(i[1]) + '\n'
            text += '+\n'
            text += str(i[2]) + '\n'
            f.write(text)


def read_fastq(fname):
    '''
    Reads list of sequences from a fastq file format.
    fname = fname to read the information

    returns list of names, sequence, quality
    '''
    with decompress(fname) as f2:
        rec = Bio.SeqIO.parse(f2, 'fastq')
        # iterate through the record
        data = []
        for i in rec:
            x = i.format('fastq').split('\n')
            data.append([i.description, str(i.seq), x[-2]])

    # raise error if empty
    if len(data) == 0:
        raise ValueError('Failed to read '+str(fname))

    # format to pandas dataframe
    col = ['name', 'sequence', 'quality']
    data = pd.DataFrame(data, columns=col)
    return data


def read_fasta(fname):
    '''
    Reads list of sequences from a fasta file format
    fname = file to read the information
    returns list of sequence and names
    '''
    data = []
    with decompress(fname) as f2:
        rec = Bio.SeqIO.parse(f2, 'fasta')
        # iterate through the record
        data = [[i.description, str(i.seq)] for i in rec]

    # raise error if empty
    if len(data) == 0:
        raise ValueError('Failed to read '+str(fname))

    # format to pandas dataframe
    col = ['name', 'sequence']
    data = pd.DataFrame(data, columns=col)
    return data


def write_fasta(fname, data):
    '''
    Write list of sequences into a fasta file format.
    data =  dataframe columns [name, sequence]
    fname = file to write the information
            .bz2 and .gz extension can be added to filename if compression is desired
    '''
    data = data[['name', 'sequence']].values
    with compress(fname) as f:
        for i in data:
            text = '>' + str(i[0]) + '\n'
            text += str(i[1]) + '\n'
            f.write(text)


def read_fastx(fname):
    '''
    Reads data from fasta or fastq
    returns filename
    '''
    # write data to fastq file
    try:
        return read_fastq(fname)
    except Exception:
        logging.info('failed to read '+str(fname)+' as fastq')

    try:
        return read_fasta(fname)
    except Exception:
        logging.info('failed to read '+str(fname)+' as fasta')
    # throws an error to do
    sys.exit(1)


def write_fastx(prefix, data):
    '''
    Writes data to fasta or fastq
    returns filename
    '''
    # write data to fastq file
    try:
        fname = prefix+'.fq'
        write_fastq(fname, data)
        return fname
    except Exception:
        logging.info('failed to write fastq')

    # write data to fasta file
    try:
        fname = prefix+'.fa'
        write_fasta(fname, data)
        return fname
    except Exception:
        logging.info('failed to write fasta')

    # throws an error to do
    sys.exit(1)


def parse_SAM(fname):
    '''
    This function parses the information from a SAM file and returns a dataframe
    fname = filename of .sam file
    returns dataframe with columns
    ['query_id','database_id','orientation','flag','t_start','t_end',
    'mapq','CIGAR','AS','XS','XN','XM','XO','XG','NM']
    '''
    fsize = os.path.getsize(fname)
    record = False
    col1 = ['query_id', 'database_id', 'orientation', 'flag', 't_start', 't_end', 'mapq', 'CIGAR']
    col2 = ['AS', 'XS', 'XN', 'XM', 'XO', 'XG', 'NM']
    key = {col2[i]: i for i in range(0, len(col2))}
    data = []
    with open(fname, 'r') as f:
        while fsize > f.tell():
            text = f.readline()
            # dont record anything until we read past the header
            if ~record and text[:3] == '@PG':
                record = True
                text = f.readline()
            # start recording data
            if record:
                out, info = get_SAM_info(text, key)
                if out[3][-2] == '1':
                    text = f.readline()
                    out2, info2 = get_SAM_info(text, key)
                    out[0] += ' '+out2[0]
                    out[7] += ' '+out2[7]
                    info = info+info2
                data.append(out+[j for j in info])
    return pd.DataFrame(data, columns=col1+col2)


def get_SAM_info(read, key):
    '''
    Function to parse useful info from a line of the SAM file
    read = each line of sam file
    key = columns of data to look for
    '''
    read = read.split('\t')
    qname = read[0]
    flag = '{0:10b}'.format(int(read[1]))  # decode the sam flag via binary to see if reverse complement is aligned
    rname = read[2]
    if flag[-5] == '1':
        orientation = '-'
        pos2 = int(read[3])
        pos1 = pos2 + int(read[8])
    else:
        orientation = '+'
        pos1 = int(read[3])
        pos2 = pos1 + int(read[8])
    mapq = int(read[4])
    cigar = read[5]
    # get another sam info
    info = np.array([0]*len(key))
    if rname != '*':
        for i in range(11, len(read)):
            x = read[i].split(':')  # parse it backwards
            if x[0] in key.keys():
                info[key[x[0]]] = x[-1]
    data = [qname, rname, orientation, flag, pos1, pos2, mapq, cigar]
    return data, info


def read_PAF(fname):
    '''
    Read a minimap2 PAF file.

    Args:
        fname: PAF filename.

    Returns:
        A pandas dataframe of alignment results.
    '''
    col1 = ['query_id', 'q_len', 'q_start', 'q_end', 'orientation',
            'database_id', 't_len', 't_start', 't_end', 'match', 'tot', 'mapq']
    col2 = ['tp', 'cm', 's1', 's2', 'NM', 'AS', 'ms', 'nn', 'rl', 'cg']
    key = {col2[i]: i for i in range(len(col2))}

    out = []
    with decompress(fname) as f:
        for line in f:
            line = line.rstrip('\n')
            if len(line) == 0 or line[0] == '#':
                continue
            fields = line.split('\t')
            if len(fields) < len(col1):
                continue
            info = [0]*len(col2)
            for field in fields[len(col1):]:
                x = field.split(':', 2)
                if len(x) == 3 and x[0] in key:
                    info[key[x[0]]] = x[2]
            out.append(fields[:len(col1)]+info)

    data = pd.DataFrame(out, columns=col1+col2)
    data = data.rename(columns={'cg': 'CIGAR'})
    if 'CIGAR' in data.columns:
        data['CIGAR'] = [str(i).replace('\n', '') for i in data['CIGAR']]
    return format_alignment(data)


# backwards-compatible alias
read_paf = read_PAF


def read_SAM(fname, mode='r'):
    '''
    Read a SAM or BAM file.

    Args:
        fname: SAM or BAM filename.
        mode: 'r' for SAM or 'rb' for BAM. BAM files are detected by
            filename when mode is 'r'.

    Returns:
        A pandas dataframe of alignment results.

    Notes:
        Reading BAM files requires pysam. For SAM input the reference
        length (t_len) is the aligned reference span, since a SAM record
        does not contain the full reference length.
    '''
    if mode == 'rb' or str(fname).lower().endswith('.bam'):
        return read_BAM(fname)

    col1 = ['query_id', 'q_len', 'q_start', 'q_end', 'orientation',
            'database_id', 't_len', 't_start', 't_end', 'CIGAR', 'mapq']
    col2 = ['NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 's2', 'de', 'rl']
    key = {col2[i]: i for i in range(len(col2))}

    out = []
    with decompress(fname) as f:
        for line in f:
            if line[0] == '@':
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11 or fields[2] == '*':
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]
            q_len = len(seq) if seq != '*' else 0
            q_start, q_end = cigar_query_bounds(cigar, q_len)
            t_start = max(pos-1, 0)
            t_len = cigar_ref_span(cigar)
            t_end = t_start+t_len
            orientation = '-' if flag & 16 else '+'

            info = [0]*len(col2)
            for field in fields[11:]:
                x = field.split(':', 2)
                if len(x) == 3 and x[0] in key:
                    info[key[x[0]]] = x[2]
            out.append([qname, q_len, q_start, q_end, orientation,
                        rname, t_len, t_start, t_end, cigar, mapq]+info)

    data = pd.DataFrame(out, columns=col1+col2)
    data = add_match_stats(data)
    return format_alignment(data)


def read_BAM(fname):
    '''
    Read a BAM file with pysam.

    Args:
        fname: BAM filename.

    Returns:
        A pandas dataframe of alignment results.
    '''
    col1 = ['query_id', 'q_len', 'q_start', 'q_end', 'orientation',
            'database_id', 't_len', 't_start', 't_end', 'CIGAR', 'mapq']
    col2 = ['NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 's2', 'de', 'rl']

    out = []
    samfile = pysam.AlignmentFile(str(fname), 'rb', check_sq=False)
    try:
        for read in samfile:
            if read.reference_name is None:
                continue
            tags = dict(read.tags)
            q_len = read.infer_read_length()
            if q_len is None:
                q_len = read.query_length if read.query_length is not None else 0
            t_start = read.reference_start
            t_end = read.reference_end
            if t_end is None:
                t_end = t_start
            t_len = read.reference_length
            if t_len is None:
                t_len = t_end-t_start
            out.append([read.qname, q_len, read.qstart, read.qend,
                        '+' if read.is_forward else '-',
                        read.reference_name, t_len, t_start, t_end,
                        read.cigarstring if read.cigarstring is not None else '',
                        read.mapq]+[tags.get(i, 0) for i in col2])
    finally:
        samfile.close()

    data = pd.DataFrame(out, columns=col1+col2)
    data = add_match_stats(data)
    return format_alignment(data)


def add_match_stats(data):
    '''
    Add match and tot columns to a SAM derived alignment dataframe.

    Args:
        data: dataframe with columns [q_start, q_end, NM].

    Returns:
        A dataframe with columns [tot, match] added.
    '''
    if len(data) == 0:
        data['tot'] = pd.Series(dtype=int)
        data['match'] = pd.Series(dtype=int)
        return data
    nm = pd.to_numeric(data['NM'], errors='coerce').fillna(0).astype(int)
    q_start = pd.to_numeric(data['q_start'], errors='coerce').fillna(0).astype(int)
    q_end = pd.to_numeric(data['q_end'], errors='coerce').fillna(0).astype(int)
    data['tot'] = q_end-q_start
    data['match'] = data['tot']-nm
    return data


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
            span += int(n)
    return span


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

    if all(i in data.columns for i in ['q_len', 't_len', 'NM']):
        data = get_similarity(data)
    elif 'similarity' not in data.columns:
        data['similarity'] = pd.Series(dtype=float)
    return data


def get_similarity(data):
    '''
    Add a similarity metric to an alignment dataframe.

    Args:
        data: dataframe with columns [NM, q_len, t_len].

    Returns:
        A dataframe with similarity added to the columns.
    '''
    if len(data) == 0:
        data['similarity'] = pd.Series(dtype=float)
        return data
    q_len = pd.to_numeric(data['q_len'], errors='coerce').fillna(0).values
    t_len = pd.to_numeric(data['t_len'], errors='coerce').fillna(0).values
    nm = pd.to_numeric(data['NM'], errors='coerce').fillna(0).values
    v1 = np.min([q_len, t_len], axis=0)
    v1 = np.where(v1 == 0, 1, v1)
    data['similarity'] = 1-nm/v1
    data.loc[data['similarity'] < 0, 'similarity'] = 0
    return data
