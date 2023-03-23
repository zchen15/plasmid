#!/usr/bin/env python
# tests for aligner.py

import plasmid
import glob

def test_read_genbank():
    # Read a genbank file
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)        
    assert type(x).__name__ == 'Plasmid'
    return x

   
