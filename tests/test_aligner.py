#!/usr/bin/env python
# tests for fileIO.py

import plasmid
import pandas as pd
import glob
import os

def test_mappy():
    # Read a fasta and fastq file
    f1 = 'data/COI/CO1_0.fq.bz2'
    f2 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f2 = plasmid.fileIO.read_fastq(f2)
    f1 = f1.iloc[:10]
    f2 = f2.iloc[:10]
    # perform cross alignment
    x = plasmid.Aligner().mappy(f1, f2)
    assert type(x).__name__ == 'DataFrame'
    col = ['query_id','q_len','q_start','q_end','ref_id','t_len','t_start','t_end','strand','read_num',
           'blen','mlen','match_score','NM','mapq','is_primary','CIGAR']
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
    col = ['query_id','ref_id','strand','AS',
           'q_len','q_start','q_end',
           't_len','t_start','t_end',
           'match','mismatch','indel','CIGAR']
    for c in col:
        assert c in x.columns
    return x

def test_pyspoa():
    f1 = 'data/mRNA.fa'
    # perform cross alignment
    seq = ['atccggttcc',
           'atcccccc',
           'ataattcc']
    x = plasmid.Aligner().spoa(seq)
    assert type(x['consensus']) == str
    assert type(x['msa']) == list
    for i in x['msa']:
        assert type(i) == str
    return x

def test_minimap2():
    f1 = 'data/COI/CO1_0.fq.bz2'
    f2 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f2 = plasmid.fileIO.read_fastq(f2)
    f1 = f1.iloc[:10]
    f2 = f2.iloc[:10]
    
    # perform cross alignment
    aligner = plasmid.Aligner()
    #aligner.params['minimap']['config'] = '-x map-ont -P'
    x = aligner.minimap2(query=f1, database=f2)
    assert type(x).__name__ == 'DataFrame'
    
    # test index building
    fn_index = aligner.minimap2_get_index(f2)
    assert os.path.exists(fn_index)    
    x = aligner.minimap2(query=f1)
    
    col1 = ['query_id','q_len','q_start','q_end','orientation',
            'database_id','t_len','t_start','t_end','match','tot','mapq']
    col2 = ['tp','cm','s1','s2','NM','AS','ms','nn','rl','CIGAR']
    col = col1+col2
    for c in col:
        assert (c in x.columns)
    return x

def test_bwa():
    f1 = 'data/COI/CO1_0.fq.bz2'
    f2 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f2 = plasmid.fileIO.read_fastq(f2)
    f1 = f1.iloc[:10]
    f2 = f2.iloc[:10]
    
    # test index building
    aligner = plasmid.Aligner()
    #fn_index = aligner.bwa_get_index(f2)
    x = aligner.bwa(query=f1, database=f2)
    
    col1 = ['query_id','database_id','orientation','flag','t_start','t_end','mapq','CIGAR']
    col2 = ['AS','XS','XN','XM','XO','XG','NM']
    col = col1+col2
    for c in col:
        assert (c in x.columns)
    return x

def test_bowtie2():
    f1 = 'data/COI/CO1_0.fq.bz2'
    f2 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f2 = plasmid.fileIO.read_fastq(f2)
    f1 = f1.iloc[:10]
    f2 = f2.iloc[:10]
    
    # test index building
    aligner = plasmid.Aligner()
    #fn_index = aligner.bwa_get_index(f2)
    x = aligner.bowtie2(query=f1, database=f2)
    
    col1 = ['query_id','database_id','orientation','flag','t_start','t_end','mapq','CIGAR']
    col2 = ['AS','XS','XN','XM','XO','XG','NM']
    col = col1+col2
    for c in col:
        assert (c in x.columns)
    return x

def test_get_fastq_statistics():
    f1 = 'data/COI/CO1_0.fq.bz2'
    f1 = plasmid.fileIO.read_fastq(f1)
    f1 = f1.iloc[:10]
    x = plasmid.Aligner.get_fastq_statistics(f1)
    col = ['Q_min','Q_lower','Q_median','Q_upper','Q_max']
    for c in col:
        assert (c in x.columns)
    return x

def test_search_ORF():
    x = 'ATGGGAACCATGGGATAA'
    seq = ''+x
    seq+= 'A'+x
    seq+= 'A'+x
    
    # check in fwd direction
    df = plasmid.Aligner.search_ORF(seq)
    check_ORF(df, seq)
    
    # check reverse direction
    rev_seq = plasmid.misc.revcomp(seq)
    df = plasmid.Aligner.search_ORF(rev_seq)
    check_ORF(df, rev_seq)
    
def check_ORF(df, seq):
    # print debugging output
    print(df)
    for i in range(3):
        print('frame',i,plasmid.misc.translate(seq[i:]))
    print(seq)
    
    for s1, s2, AA, strand in df[['start','stop','sequence','orientation']].values:
        y = seq[s1:s2]
        print('sequence=',y)
        if strand==1:
            y = plasmid.misc.translate(y)
        else:
            y_rev = plasmid.misc.revcomp(y)
            y = plasmid.misc.translate(y_rev)
        print('AA seq_true',y)
        print('AA seq',AA)
        assert y == AA

def test_filter_idx():
    x = [['HAP01','HAP02',1],
         ['HAP01','HAP01',2],
         ['HAP02','HAP01',3],
         ['HAP02','HAP02',4]]
    col = ['q','d','val']
    x = pd.DataFrame(x, columns=col)
    y = plasmid.Aligner.filter_idx(x, ['q'], 'val', 'idxmax')
    print('x',x)
    print('y',y)
    c = y['q']!=y['d']
    assert sum(c)==0
