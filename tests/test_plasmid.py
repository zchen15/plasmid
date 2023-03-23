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
    
    # index set
    val = 'hello'
    col = 'locus_tag'
    i = 1
    x[col, i] = val
    assert x[col].values[i] == val 
    
    val = 'hello'
    col = 'color'
    i = 1
    x[i, col] = val
    assert x[col].values[i] == val 

    # bool setting
    val = 'blue'
    col = 'color'
    idx = x['type'] == 'terminator'
    x[idx, col] = val
    x = x['color'].values[idx]
    for y in x:
        assert y==val
    
def test_splice():
    # Checks splicing of sequences
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    i = 1
    name1, s1, s2 = x[['locus_tag','start','end']].values[i]
    seq1 = x[s1:s2].splice().__str__()
    seq2 = x.__str__()[s1:s2]

    assert seq1==seq2

def test_replace():
    # Checks sequence replacement
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    i = 2
    val = 'TTTTTTAAAAAA'
    
    # check initial values
    s1 = x[i]['start'].values[0]
    s2 = s1 + len(val)
    print(s1, s2)

    # set to new value
    x[i] = val

    # check new values
    seq = x.__str__()[s1:s2]
    print(val, seq)
    assert val == seq
    
def test_merge():
    # Checks merger of records todo
    fname = 'data/dcas9_LacI.gb'
    x = plasmid.read_genbank(fname)
    assert True


    
