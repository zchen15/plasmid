#!/usr/bin/env python
# tests for fileIO.py

import plasmid
import glob

def test_minimap2():
    # Read a fasta and fastq file
    f1 = 'data/COI/CO1_0.fq.bz2'
    f2 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f2 = plasmid.fileIO.read_fastq(f2)
    f1 = f1.iloc[:10]
    f2 = f2.iloc[:10]
    # perform cross alignment
    x = plasmid.Aligner().minimap2(f1, f2)
    assert type(x).__name__ == 'DataFrame'
    col = ['query_id','q_len','q_start','q_end','ref_id','t_len','t_start','t_end','strand','read_num',
           'blen','mlen','match_score','NM','mapq','is_primary','cigar']
    for c in col:
        assert c in x.columns
    return x

def test_parasail():
    f1 = 'data/COI/CO1_0.fq.bz2'
    f2 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f2 = plasmid.fileIO.read_fastq(f2)
    f1 = f1.iloc[:10]
    f2 = f2.iloc[:10]
    # perform cross alignment
    x = plasmid.Aligner().parasail(f1, f2)
    assert type(x).__name__ == 'DataFrame'
    col = ['query_id','q_len','ref_id','t_len','strand','AS','NM','match','mismatch','indel','cigar']
    for c in col:
        assert c in x.columns
    return x

def test_pyspoa():
    f1 = 'data/mRNA.fa'
    # perform cross alignment
    x = plasmid.Aligner().spoa(f1, f2)
    assert type(x).__name__ == 'DataFrame'
    col = ['query_id','q_len','ref_id','t_len','strand','AS','NM','match','mismatch','indel','cigar']
    for c in col:
        assert c in x.columns
    return x
