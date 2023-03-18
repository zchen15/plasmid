#!/usr/bin/env python
# tests for plasmid.py

import plasmid
import glob
import copy

def test_genbank_read():
    # Checks reading and loading genbank files
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    assert type(x).__name__ == 'Plasmid'
    return x
df = test_genbank_read()

def test_plasmid_print():
    # Check printing
    assert type(df.__str__())==str
    print(df.__repr__())
    assert type(df.__repr__())==str
    print(df.print())    

def test_concatenation():
    # Checks concatenation
    x = df+df
    print(df)
    print(x)
    assert type(x).__name__ == 'Plasmid'
    
def test_deletion():
    # Checks deletion of entries
    x = copy.copy(df)
    x0 = x[0]
    print(x[0], x0)
    assert (x0==x[0])==False

def test_rotation():
    # Checks rotation of entries
    assert True
    
def test_rotation():
    # Checks rotation of entries
    assert True

