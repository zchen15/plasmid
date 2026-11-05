#!/usr/bin/env python

import Bio
import Bio.Seq
import pandas as pd
import re

# type hint imports
from typing import Union

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
    comp = str.maketrans("ACGTRYSWKMBVDHN", "TGCAYRSWMKVBHDN")
    
    hits = []
    # Lookahead regex (?=...) finds overlapping matches without consuming characters
    pattern = re.compile(f"(?={qry})")
    
    for seq, orient in [(ref, 1), (ref.translate(comp)[::-1], -1)]:
        if orient == -1 and fwd_only:
            break
            
        search_seq = seq + seq if circular else seq
        
        # Iterate through all regex matches
        for m in pattern.finditer(search_seq):
            idx = m.start()
            
            # Stop searching if we've wrapped around the circular genome
            if circular and idx >= seq_len:
                break
            
            end = idx + len(qry)
            if circular:
                end = end % seq_len or seq_len 
                
            # Adjust coordinates for the reverse strand
            s, e = (idx, end) if orient == 1 else (seq_len - end, seq_len - idx)
            hits.append([s, e, orient])

    return pd.DataFrame(hits, columns=["t_start", "t_end", "orientation"])

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
        circular: Treat the reference as a circular genome.

    Returns:
        DataFrame with columns ``t_start``, ``t_end``, ``orientation``
        (1 = forward, -1 = reverse complement), sorted by ``t_start``.
    """
    empty = pd.DataFrame(columns=["t_start", "t_end", "orientation"])
    
    qry = str(qry).upper().replace("-", "").replace(" ", "")
    ref = str(ref).upper().replace("-", "").replace(" ", "")
    
    if not qry or not ref:
        return empty
    
    seq_len = len(ref)
    results = []
    
    for strand, seq in [(1, ref), (-1, str(Seq(ref).reverse_complement()))]:
        if strand == -1 and fwd_only:
            break
            
        for frame in range(3):
            search_seq = (seq + seq)[frame:] if circular else seq[frame:]
            if not search_seq:
                continue
                
            translated = str(Seq(search_seq).translate())
            if not translated:
                continue
                
            hits = match_DNA(qry, translated, fwd_only=True, circular=False, exact=True)
            if hits.empty:
                continue
            
            # Convert to DNA coordinates
            hits["t_start"] = hits["t_start"] * 3 + frame
            hits["t_end"] = hits["t_end"] * 3 + frame
            
            # Handle circular
            if circular:
                hits = hits[hits["t_start"] < (seq_len - frame) // 3 * 3 + frame]
                hits["t_end"] = hits["t_end"] % seq_len
                hits.loc[hits["t_end"] == 0, "t_end"] = seq_len
            
            # Reverse strand
            if strand == -1:
                hits["t_start"] = (seq_len - hits["t_end"]) % seq_len
                hits["t_end"] = (seq_len - hits["t_start"]) % seq_len
                hits.loc[hits["t_start"] == 0, "t_start"] = seq_len
                hits.loc[hits["t_end"] == 0, "t_end"] = seq_len
            
            hits["orientation"] = strand
            results.append(hits)
    
    if not results:
        return empty
    
    return pd.concat(results, ignore_index=True).sort_values("t_start", kind="stable").reset_index(drop=True)


def get_ORF(seq: str, start: str = "M", stop: str = "*") -> pd.DataFrame:
    """Get list of ORFs based on start and stop codons.

    Args:
        seq: Amino-acid sequence (already translated).
        start: Start character (default ``'M'``).
        stop: Stop character (default ``'*'``).

    Returns:
        DataFrame with columns ``start``, ``stop``, ``sequence``.
    """
    pattern = re.escape(start) + ".*?" + re.escape(stop)
    matches = [
        [m.start(), m.end(), m.group()]
        for m in re.finditer(pattern, str(seq))
    ]
    return pd.DataFrame(matches, columns=["start", "stop", "sequence"])

def search_ORF(
    seq: str,
    table: Union[str, int] = "Standard",
    fwd_only: bool = False,
) -> pd.DataFrame:
    """Generate a dataframe of open reading frames for a nucleotide sequence.

    Args:
        seq: Input nucleotide sequence.
        table: Codon table to use (name or NCBI index).
        fwd_only: If True, search only in forward orientation.

    Returns:
        DataFrame of ORFs with start/stop in nucleotide coordinates.
    """
    cols = ["start", "stop", "sequence", "frame", "length", "orientation"]
    seq_str = str(seq)
    seq_len = len(seq_str)

    if seq_len < 3:
        return pd.DataFrame(columns=cols)

    seq_obj = Seq(seq_str)
    strands = [(1, seq_obj)]
    if not fwd_only:
        strands.append((-1, seq_obj.reverse_complement()))

    results: list[tuple] = []

    for orientation, strand_seq in strands:
        for frame in range(3):
            aa = str(strand_seq[frame:].translate(table=table))
            orf_df = get_ORF(aa)

            for m in orf_df.itertuples(index=False):
                start_aa, stop_aa, sequence = m
                start = start_aa * 3 + frame
                stop = stop_aa * 3 + frame

                if orientation == -1:
                    start, stop = seq_len - stop, seq_len - start

                results.append((start, stop, sequence, frame, stop - start, orientation))

    return pd.DataFrame(results, columns=cols)
