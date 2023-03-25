#!/usr/bin/env python
# tests for designer.py

import plasmid as pge
import glob

def get_parts01():
    # slice out the RFP gene
    RFP = pge.read_genbank('data/dcas9_RFP.gb')
    RFP = RFP[RFP['locus_tag'].str.contains('mRFP')].splice()
    # slice out the ribosome binding site
    RBS = pge.read_genbank('data/xRFP.gb')
    RBS = RBS[RBS['locus_tag'].str.contains('BBa_B0034')].splice()
    # slice out the promoter
    pLac = pge.read_genbank('data/xRFP.gb')
    pLac = pLac[pLac['locus_tag'].str.contains('pLac')].splice()

    # assemble the promoter, rbs, and mRFP
    insert = pLac + 'gagacc' + RBS + 'ggtctc'
    return insert, RFP


def test_xtPCR():
    pcr = pge.Designer()
    pcr.params['xtPCR']['Tm'] = 55         # target annealing temperature for xtPCR
    pcr.params['xtPCR']['len'] = [15, 60]  # defines the [min, max] primer lengths

    insert, RFP = get_parts01()
    res = pcr.xtPCR(insert, RFP, ' ')
    print(res)
    print(res.values)
    assert type(res).__name__ == 'DataFrame'

def get_parts02():
    # slice out the LacI gene
    LacI = pge.read_genbank('data/xRFP.gb')
    LacI = LacI[LacI['locus_tag'].str.contains('LacI')].splice()

    # slice out the RFP gene
    RFP = pge.read_genbank('data/dcas9_RFP.gb')
    RFP = RFP[RFP['locus_tag'].str.contains('mRFP')].splice()

    # slice out the origin of replication
    df = pge.read_genbank('data/xRFP.gb')
    vec = df[df['locus_tag'].str.contains('pSC101')]
    start = vec['start'][0]
    stop = vec['end'][0]
    vec = df[start:stop]

    seq = []
    seq+= [[' ',LacI,'AAAActttt']]
    seq+= [[' ',RFP,'CGCCctttt']]
    seq+= [[' ',vec,'GGGGctttt']]
    # return data
    return seq

def test_Gibson():
    seq = get_parts02()

    pcr = pge.Designer()
    pcr.params['gibson']['Tm'] = 50     # target annealing temperature of gibson fragments    
    pcr.params['gibson']['window'] = 30 # +/i window in bp around frag edges to look for gibson overlap
    pcr.params['gibson']['len'] = 20    # length of gibson overlap

    pcr.params['xtPCR']['Tm'] = 55         # target annealing temperature for xtPCR
    pcr.params['xtPCR']['len'] = [15, 60]  # defines the [min, max] primer lengths
    pcr.params['xtPCR']['nM'] = [20, 500] # defines the [seed, finisher] primer conc in nM

    res = pcr.Gibson(seq)
    print(res)
    print(res.values)
    assert type(res).__name__ == 'DataFrame'

def test_GoldenGate():
    seq = get_parts02()

    pcr = pge.Designer()
    pcr.params['goldengate']['window'] = 20 # +/i window in bp around frag edges to look for overlap
    pcr.params['goldengate']['ggN'] = 4     # length of golden gate overlap
    pcr.params['goldengate']['ggsite'] = 'GGTCTCc'     # golden gate enzyme site
    pcr.params['goldengate']['padding'] = 'atatatatgg' # padding around the golden gate site
    pcr.params['xtPCR']['len'] = [15, 60]  # defines the [min, max] primer lengths
    pcr.params['xtPCR']['nM'] = [20, 500] # defines the [seed, finisher] primer conc in nM
    pcr.params['xtPCR']['Tm'] = 55 # defines the [seed, finisher] primer conc in nM

    res = pcr.GoldenGate(seq)
    print(res)
    print(res.values)
    assert type(res).__name__ == 'DataFrame'

