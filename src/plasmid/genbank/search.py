#!/usr/bin/env python
"""Sequence search and ORF-finding helpers.

This module contains the sequence-search functions formerly defined in
``plasmid.py``: DNA/amino-acid matching, substring search, and open
reading frame detection.
"""

# science libs
import Bio
import Bio.Seq
import pandas as pd


def find_all_substrings(qry: str, ref: str) -> list[int]:
    """Find all (possibly overlapping) occurrences of ``qry`` in ``ref``.

    Args:
        qry: Query substring.
        ref: Reference string.

    Returns:
        List of start indices.
    """
    out: list[int] = []
    i = ref.find(qry)
    while i != -1:
        out.append(i)
        i = ref.find(qry, i + 1)
    return out


def match_DNA(
    qry: str,
    ref: str,
    fwd_only: bool = False,
    circular: bool = False,
    exact: bool = False,
) -> pd.DataFrame:
    """Find all instances of the query sequence in the reference.

    Args:
        qry: DNA subsequence to find.
        ref: DNA sequence containing the query subsequence.
        fwd_only: Search in forward orientation 5'->3' only.
        circular: Search string as though it's a circular genome.
        exact: Perform exact substring search.

    Returns:
        DataFrame with columns ``t_start``, ``t_end``, ``orientation``.
    """
    qry = str(qry).upper().replace("-", "").replace(" ", "")
    ref = str(ref).upper()
    if not qry:
        return pd.DataFrame(columns=["t_start", "t_end", "orientation"])

    seq_len = len(ref)

    def _forward_hits(seq: str, is_circular: bool) -> list[list[int]]:
        """Collect forward-orientation hits as ``[start, end, orientation]``."""
        if is_circular:
            seq = seq + seq
        hits = []
        for i in find_all_substrings(qry, seq):
            if is_circular and i >= seq_len:
                continue
            end = i + len(qry)
            if is_circular:
                end = end % seq_len
                if end == 0:
                    end = seq_len
            hits.append([i, end, 1])
        return hits

    fwd_out = pd.DataFrame(_forward_hits(ref, circular),
                           columns=["t_start", "t_end", "orientation"])
    if fwd_only:
        return fwd_out

    rev_hits = _forward_hits(str(Bio.Seq.Seq(ref).reverse_complement()), circular)
    rev_out = pd.DataFrame(
        [[seq_len - end, seq_len - start, -1] for start, end, _ in rev_hits],
        columns=["t_start", "t_end", "orientation"],
    )
    return pd.concat([fwd_out, rev_out], ignore_index=True)


def match_AA(
    qry: str,
    ref: str,
    fwd_only: bool = False,
    circular: bool = False,
) -> pd.DataFrame:
    """Find all instances of the query (amino acid) sequence in the reference.

    Args:
        qry: Amino acid subsequence to find.
        ref: DNA sequence containing the query subsequence.
        fwd_only: Search in forward orientation only.
        circular: Search string as though it's a circular genome.

    Returns:
        DataFrame with columns ``t_start``, ``t_end``, ``orientation``.
    """
    qry = str(qry).upper().replace("-", "").replace(" ", "")
    ref = str(ref).upper()

    if fwd_only:
        if circular:
            out = match_AA(qry, ref + ref, fwd_only=True, circular=False)
            L = len(ref)
            out = out[out["t_start"] < L].copy()
            out["t_end"] = out["t_end"] % L
            out.loc[out["t_end"] == 0, "t_end"] = L
            return out

        seq = Bio.Seq.Seq(ref)
        frames = []
        for i in range(3):
            aa = str(seq[i:].translate())
            df = match_DNA(qry, aa, fwd_only=True, circular=False, exact=True)
            if len(df) == 0:
                continue
            df["t_start"] = df["t_start"] * 3 + i
            df["t_end"] = df["t_start"] + len(qry) * 3
            frames.append(df)
        if not frames:
            return pd.DataFrame(columns=["t_start", "t_end", "orientation"])
        return pd.concat(frames, ignore_index=True)

    fwd_out = match_AA(qry, ref, fwd_only=True, circular=circular)
    rev_ref = str(Bio.Seq.Seq(ref).reverse_complement())
    rev_out = match_AA(qry, rev_ref, fwd_only=True, circular=circular).copy()
    rev_out["orientation"] = -1
    rev_out["t_start"] = len(ref) - rev_out["t_start"]
    rev_out["t_end"] = len(ref) - rev_out["t_end"]
    rev_out = rev_out.rename(columns={"t_end": "t_start", "t_start": "t_end"})
    return pd.concat([fwd_out, rev_out], ignore_index=True)


def get_ORF(seq: str, start: str = "M", stop: str = "*") -> pd.DataFrame:
    """Get list of ORFs based on start and stop codons.

    Args:
        seq: Amino-acid sequence (already translated).
        start: Start character (default ``'M'``).
        stop: Stop character (default ``'*'``).

    Returns:
        DataFrame with columns ``start``, ``stop``, ``sequence``.
    """
    seq = str(seq)
    start_positions = [i for i, c in enumerate(seq) if c == start]
    stop_positions = [i for i, c in enumerate(seq) if c == stop]

    out: list[list] = []
    j = 0
    n_stops = len(stop_positions)
    for s1 in start_positions:
        while j < n_stops and stop_positions[j] <= s1:
            j += 1
        if j >= n_stops:
            break
        s2 = stop_positions[j]
        out.append([s1, s2 + 1, seq[s1 : s2 + 1]])
    return pd.DataFrame(out, columns=["start", "stop", "sequence"])


def search_ORF(seq: str, table: str = "Standard", fwd_only: bool = False) -> pd.DataFrame:
    """Generate a dataframe of open reading frames for a nucleotide sequence.

    Args:
        seq: Input nucleotide sequence.
        table: Codon table to use.
        fwd_only: If True, search only in forward orientation.

    Returns:
        DataFrame of ORFs with start/stop in nucleotide coordinates.
    """
    seq_str = str(seq)
    seq_obj = Bio.Seq.Seq(seq_str)

    out = []
    for i in range(3):
        aa = str(seq_obj[i:].translate(table=table))
        df = get_ORF(aa)
        if len(df) == 0:
            continue
        df["start"] = df["start"] * 3 + i
        df["stop"] = df["stop"] * 3 + i
        df["frame"] = i
        df["length"] = [len(s) for s in df["sequence"]]
        out.append(df)

    if not out:
        cols = ["start", "stop", "sequence", "frame", "length", "orientation"]
        df = pd.DataFrame(columns=cols)
    else:
        df = pd.concat(out, ignore_index=True)
    df["orientation"] = 1
    if fwd_only:
        return df

    rev_seq = str(seq_obj.reverse_complement())
    df2 = search_ORF(rev_seq, table=table, fwd_only=True)
    df2["orientation"] = -1
    stop_vals = df2["stop"].values
    df2["stop"] = len(seq_str) - df2["start"].values
    df2["start"] = len(seq_str) - stop_vals
    return pd.concat([df, df2], ignore_index=True)
