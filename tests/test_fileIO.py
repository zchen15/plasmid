#!/usr/bin/env python
# tests for aligner.py

from plasmid.aligners import io
import glob

def test_read_fastq():
    # Read a fastq file
    fname = 'data/coi.fq.bz2'
    x = io.read_fastq(fname)
    assert type(x).__name__ == 'DataFrame'
    col = ['name','sequence','quality']
    for c in col:
        assert c in x.columns
    assert len(x) > 0

def test_read_fasta():
    # Read a fasta
    fname = 'data/mRNA.fa'
    x = io.read_fasta(fname)
    assert type(x).__name__ == 'DataFrame'
    col = ['name','sequence']
    for c in col:
        assert c in x.columns

def test_read_fastx():
    # Read a fasta
    fname = 'data/mRNA.fa'
    x = io.read_fasta(fname)
    assert type(x).__name__ == 'DataFrame'
    col = ['name','sequence']
    for c in col:
        assert c in x.columns

    # Read a fastq file
    fname = 'data/coi.fq.bz2'
    x = io.read_fastq(fname)
    assert type(x).__name__ == 'DataFrame'
    col = ['name','sequence','quality']
    for c in col:
        assert c in x.columns
    assert len(x) > 0

