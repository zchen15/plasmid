#!/usr/bin/env python
# tests for aligner.py

import plasmid
import glob

def test_read_fasta():
    # Read a fasta and fastq file
    fname = 'data/COI/CO1_0.fq.bz2'
    x = plasmid.fileIO.read_fastq(fname)
    assert type(x).__name__ == 'DataFrame'
    col = ['name','sequence','quality']
    for c in col:
        assert c in x.columns
    
    # Read a fasta
    fname = 'data/mRNA.fa'
    x = plasmid.fileIO.read_fasta(fname)
    assert type(x).__name__ == 'DataFrame'
    col = ['name','sequence']
    for c in col:
        assert c in x.columns

