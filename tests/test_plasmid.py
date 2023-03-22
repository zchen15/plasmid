#!/usr/bin/env python
# tests for plasmid.py

import plasmid
import glob
import copy

def test_read_genbank():
    # Read a genbank file
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)        
    assert type(x).__name__ == 'Plasmid'
    return x
df = test_read_genbank()

def test_read_fasta():
    # Read csv and fasta
    fname = glob.glob('data/*.csv')[0]
    y = plasmid.plasmid.read_fasta(fname)
    fname = glob.glob('data/*.fa')[0]
    y+= plasmid.plasmid.read_fasta(fname)
    assert type(y) == list
    for x in y:
        assert type(x).__name__ == 'Plasmid'

def test_plasmid_print():
    # Check printing
    assert type(df.__str__())==str
    assert type(df.print())==str
    assert type(df.__repr__())==str

def test_rotation():
    # Checks rotation of entries
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    ori = x[x['locus_tag'].str.contains('TetR')]['start'].values[0]
    x = x.set_origin(ori)
    assert 'TetR' in x['locus_tag'].values[0]
    
def test_del():
    # Checks deletion of entries
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    f1 = x[0]
    f2 = x[1]
    del(x[0])
    # debug output
    print('f1', f1.__repr__())
    print('f2', f2.__repr__())
    print('x', x.__repr__())
    c1 = f2['locus_tag'].values[0]
    c2 = x['locus_tag'].values[0]
    print('==', c1==c2)
    assert c1==c2

def test_add():
    # Checks concatenation
    x = df+df
    assert type(x).__name__ == 'Plasmid'
    N = len(df)
    f1 = x[0]
    f2 = x[N]
    print(f1.__repr__())
    print(f2.__repr__())
    c1 = f1['locus_tag'].values[0]
    c2 = f2['locus_tag'].values[0]
    assert c1==c2

def test_set():
    # Checks deletion of entries
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    assert c1==c2
    
def test_slice():
    # Checks deletion of entries
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    assert c1==c2

def test_slice():
    # Checks deletion of entries
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    assert c1==c2

    