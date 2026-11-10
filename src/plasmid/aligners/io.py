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

# system and time
import os
import re
import sys
import logging
import multiprocessing


def revcomp(seq: str) -> str:
    '''
    Get the reverse complement of a sequence.

    Args:
        seq: sequence string.

    Returns:
        Reverse complement of the sequence.
    '''
    seq = str(seq).upper()
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


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


def get_best(df, col: list, value: str = "score", opt: str = "idxmax") -> pd.DataFrame:
    '''
    Get idxmax or idxmin over a set of groups

    Args:
        col: column to group by
        value: column to sort
        opt: statistic to compute such as idxmax

    Returns:
        A pandas DataFrame of filtered values
    '''
    df = df.reset_index(drop=True)
    idx = df.groupby(by=col).agg({value: opt}).reset_index()
    return df.iloc[idx[value].values].reset_index()


def remove_overlapping_reads_N(df, metric: str = "AS", thresh: float = 0.1, N: int = 2) -> pd.DataFrame:
    '''
    Remove overlapping alignments for each concatemeric read

    Args:
        df: pandas dataframe of read alignments
        metric: metric to use for sorting
        thresh: threshold overlap allowed between adjacent alignments
        N: number of N overlaps

    Returns:
        A pandas dataframe of filtered reads
    '''
    # sort and remove overlapping reads
    cols = ["q_id", metric, "q_start"]
    df = df.sort_values(by=cols, ascending=False).reset_index(drop=True)
    q_start = df["q_start"].values
    q_end = df["q_end"].values
    match = df["match"].values

    keep = np.array([True]*len(df))
    for i in range(N):
        for j in range(i+1,N):
            minL = np.min([match[i::N], match[j::N]], axis=0)
            c1 = q_end[i::N] < (q_start[j::N] + thresh*minL)
            c2 = q_end[j::N] < (q_start[i::N] + thresh*minL)
            c = c1 | c2 | ~keep[i::N]
            keep[j::N] &= c
    return df[keep]


def remove_overlapping_reads(df, metric: str = "AS", thresh: float = 0.1) -> pd.DataFrame:
    '''
    Remove overlapping alignments for each concatemeric read

    Args:
        df: pandas dataframe from read_SAM or read_PAF
        metric: metric to use for sorting
        thresh: threshold overlap allowed between adjacent alignments

    Returns:
        A pandas dataframe of filtered reads
    '''
    # work on things with only one alignment
    x = df.groupby(by=["q_id"]).agg({"q_id":"count"})
    out = []
    for N in np.unique(x["q_id"]):
        c = df["q_id"].isin(x[x["q_id"]==N].index)
        if N == 1:
            y = df[c]
        else:
            y = remove_overlapping_reads_N(df[c], metric, thresh, int(N))
        out.append(y)
    # return results
    df = pd.concat(out)
    return df


def get_indel_list(df, key: str = "([0-9]*)([DIM=SX])"):
    '''
    Add indel info to each read

    Args:
        df: pandas dataframe of alignment information
        key: cigar string keys

    Returns:
        pandas dataframe
    '''
    # iterate and get start and end positions of mutations
    out = []
    for cigar in df["CIGAR"].drop_duplicates().values:
        x = cigar_decode(cigar, key)
        cols = ["cg_len", "cg_value"]
        x = pd.DataFrame(x, columns=cols)
        x["cg_step"] = x["cg_len"]
        c = x["cg_value"] == "I"
        x.loc[c, "cg_step"] = 0
        x["cg_start"] = np.cumsum(x["cg_step"]) - x["cg_step"]
        x.loc[c, "cg_step"] = 1
        x["cg_end"] = x["cg_start"] + x["cg_step"]
        x["CIGAR"] = cigar
        out.append(x)
    df_mut = pd.concat(out)

    # merge info
    df = df.merge(df_mut, on="CIGAR", how="left")
    df["cg_start"] += df["t_start"]
    df["cg_end"] += df["t_start"]

    # correct for circular alignments
    df["cg_start"] = df["cg_start"] % df["t_len"]
    df["cg_end"] = df["cg_end"] % df["t_len"]
    c = df["cg_end"] == 0
    df.loc[c, "cg_end"] = df.loc[c, "t_len"]

    # break around origin and merge
    c = df["cg_end"] < df["cg_start"]
    df1 = df[c].copy()
    df2 = df[c].copy()
    df1["cg_start"] = 0
    df2["cg_end"] = df2["t_len"]
    df = pd.concat([df[~c], df1, df2])
    return df


def get_indel_hist(df, wcol: str):
    '''
    Get indel array for plotting

    Args:
        df: pandas dataframe of alignment info
        wcol: weights column for each alignment

    Returns:
        pandas dataframe
    '''
    cols = ["cg_start", "cg_end", "cg_value", wcol]
    df = df.groupby(by=cols[:-1]).agg({wcol: "sum"}).reset_index()

    # get position bounds
    t_start = int(np.min(df["cg_start"]))
    t_end = int(np.max(df["cg_end"]))
    x = np.zeros([t_end - t_start, 5])
    x[:, 0] = t_start + np.arange(len(x))

    # sum over positions
    mut_dict = {"M": 4, "=": 4, "X": 3, "D": 2, "I": 1}
    for s1, s2, mut, w in df[cols].values:
        c = x[:, 0] >= s1
        c &= x[:, 0] < s2
        x[c, mut_dict[mut]] += w

    # format results into a dataframe
    cols = ["position", "insertion", "deletion", "mismatch", "match"]
    df = pd.DataFrame(x, columns=cols)
    return df


def get_count(df, idx: int, block_size: int, wcol: str, key: str):
    '''
    Calls get_indel_list and get_indel_hist for multiprocessing

    Args:
        df: pandas dataframe of alignment info
        idx: index of the dataframe block
        block_size: size of the dataframe block
        wcol: weights column for each alignment
        key: cigar string keys

    Returns:
        pandas dataframe
    '''
    if idx % 100 == 0:
        logging.info(f"processing {idx}/{block_size}")
    df = get_indel_list(df, key)
    df = get_indel_hist(df, wcol)
    if idx % 100 == 0:
        logging.info(f"finished {idx}/{block_size}")
    return df


def get_indel_freq(data, wcol: str = "weights", key: str = "([0-9]*)([DIM=SX])", block_size: int = 10000, threads: int = 4):
    '''
    Get indel frequency array for plotting

    Args:
        data: pandas dataframe of alignment info
        wcol: weights column for each alignment
        key: key used to decode the cigar string. Default="([0-9]*)([DIM=SX])"
        threads: number of parallel threads to run

    Returns:
        pandas dataframe
    '''
    df = data.copy()
    df = df.sort_values(by=["t_id", "t_start", "q_id", "q_start"]).reset_index(drop=True)

    out = []
    for t_id in df["t_id"].drop_duplicates().values:
        logging.info(f"processing {t_id}")

        # prepare data for parallel processing
        c = df["t_id"] == t_id
        x = df[c].copy()
        L = max(int(len(x) / block_size), 1)
        chunk = np.array_split(x, L)
        data_chunks = [[chunk[i], i, len(chunk), wcol, key] for i in range(len(chunk))]

        # process data in parallel using Pool
        logging.info(f"Processing {len(data_chunks)} datasets using {threads} threads")
        with multiprocessing.Pool(processes=threads) as pool:
            results = pool.starmap(get_count, data_chunks)

        # aggregate results
        results = pd.concat(results)
        cols = ["insertion", "deletion", "mismatch", "match"]
        stats = {k: "sum" for k in cols}
        results = results.groupby(["position"]).agg(stats).reset_index()
        results["t_id"] = t_id

        # compute indel frequency
        results["total_reads"] = 0
        cols = ["insertion", "deletion", "mismatch", "match"]
        for k in cols:
            results["total_reads"] += results[k]
        for k in cols:
            results[k + "_freq"] = results[k] / results["total_reads"]
        out.append(results)

    # aggregate final results
    df = pd.concat(out)
    return df


def cigar_decode(cigar: str, key: str = '([0-9]*)([DIM=SX])') -> list:
    '''
    Decode cigar strings into counts of indels and mismatches

    Args:
        cigar: cigar string
        key: key used to decrypt the cigar string

    Returns:
        list containing [[key value, counts],...]
    '''
    pattern = re.compile(key)
    return [[int(i), j] for i, j in pattern.findall(cigar)]


def cigar_to_stats(cig: str, qseq: str, rseq: str) -> list:
    '''
    Get alignment statistics associated with parasail output

    Args:
        cig: cigar string
        qseq: query sequence
        rseq: reference sequence

    Returns:
        list with columns [q_len, q_start, q_end, t_len, t_start, t_end, M, Xm, indel]
    '''
    cig = cigar_decode(cig)

    # compute match boundaries
    q_start = 0
    q_end = len(qseq)
    q_len = len(qseq)
    t_start = 0
    t_end = len(rseq)
    t_len = len(rseq)
    if cig[0][1] == 'I':
        t_start = cig[0][0]
    elif cig[0][1] == 'D':
        q_start = cig[0][0]
    if cig[-1][1] == 'I':
        t_end -= cig[-1][0]
    elif cig[-1][1] == 'D':
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


def cigar_to_alignment(query: str, ref: str, cigar: str) -> list:
    '''
    Function to align query and reference strings based on CIGAR string

    Args:
        query: query sequence. This must be relative to how it was used in aligner
        ref: reference sequence
        cigar: CIGAR string from aligner

    Returns:
        list with values [aligned query, aligned reference, alignment ticks]
    '''
    pattern = re.compile('([0-9]*)([DIMSX])')
    query = query.upper()
    ref = ref.upper()
    q_aligned = ''
    ref_aligned = ''
    align = ''
    i1 = 0
    i2 = 0
    for n, c in pattern.findall(cigar):
        span = int(n)
        if c == 'M':
            q_aligned += query[i1:i1 + span]
            ref_aligned += ref[i2:i2 + span]
            align += ''.join(['|'] * span)
            i1 += span
            i2 += span
        elif c == 'X':
            q_aligned += query[i1:i1 + span]
            ref_aligned += ref[i2:i2 + span]
            align += ''.join([' '] * span)
            i1 += span
            i2 += span
        elif c == 'S' or c == 'H':
            q_aligned += query[i1:i1 + span].lower()
            ref_aligned += ref[i2:i2 + span]
            align += ''.join([' '] * span)
            i1 += span
            i2 += span
        elif c == 'I':
            q_aligned += query[i1:i1 + span]
            ref_aligned += ''.join(['-'] * span)
            align += ''.join([' '] * span)
            i1 += span
        elif c == 'D':
            q_aligned += ''.join(['-'] * span)
            ref_aligned += ref[i2:i2 + span]
            align += ''.join([' '] * span)
            i2 += span
    return [q_aligned, ref_aligned, align]


def get_fastq_statistics(query: pd.DataFrame) -> pd.DataFrame:
    '''
    Compute quality statistics about each read

    Args:
        query: dataframe with columns [quality]

    Returns:
        original dataframe with ['Q_min','Q_lower','Q_median','Q_upper','Q_max'] appended
    '''
    out = []
    for v in query['quality'].values:
        v = phred_to_int(v)
        q = [0, 0.25, 0.5, 0.75, 1]
        y = np.quantile(v, q)
        out.append(y)

    col = ['Q_min', 'Q_lower', 'Q_median', 'Q_upper', 'Q_max']
    df = pd.DataFrame(out, columns=col)
    for c in col:
        query[c] = df[c]
    return query


def phred_to_int(x) -> list:
    '''
    Translate phred ASCII to quality integers

    Args:
        x: quality strings

    Returns:
        list of integers
    '''
    return [ord(i) - 33 for i in x]
