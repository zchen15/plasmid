#!/usr/bin/env python
# Class functions to nucleic acid design

# computing libraries
import numpy as np
import scipy as sp
import scipy.optimize
import pandas as pd

# reading and writing genbank files
import Bio
import Bio.SeqUtils

# system libraries
import time

class Nucleic:
    '''
    Class object holding information for primer design
    '''
    params = {
          'verbose':False,
          'max_opt':3600*12,
          'fstop':0.05,
          'trials':1,
          'data':'strands',
          'ckpt_int':600,
          'cache':2,
          'material':'dna',
          'wobble':True,
          'temperature':37,
          'excludes':[],
          'terminator':{'BBa_B1010':'cgccgcaaaccccgcccctgacagggcggggtttcgccgc',
                        'BBa_B1006':'aaaaaaaaaccccgcccctgacagggcggggtttttttt',
                        'BBa_K864601':'ctgtaacagagcattagcgcaaggtgatttttgtcttcttgcgctaatttttt',
                        'BBa_B1002':'cgcaaaaaaccccgcttcggcggggttttttcgc'},
          'enzymes':{'BsaI':'GGTCTC',
               'Esp3I':'CGTCTC',
               'BspMI':'ACCTGC',
               'PaqCI':'CACCTGC',
               'TelN':'TATCAGCACACAATTGCCCATTATACGCGCGTATAATGGACTATTGTGTGCTGATA'}}
    '''
           'Nt.BsmAI':'GTCTC',
           'Nb.BssSI':'CACGAG',
           'Nb.BsmI':'GAATGC',
           'Nt.BbvCI':'CCTCAGC',
           'Nt.AlwI':'GGATCNNNNN',
           'Nb.BtsI':'GCAGTGNN',
           'Nb.BsrDI':'GCAATGNN',
           'Nt.BstNBI':'GAGTCNNNNN',
           'Nb.BspQI':'GCTCTTCN',
    '''
    params['xtPCR'] = {'Tm':50, 'len':[15,60], 'nM':[20,500], 'Na':50}
    params['gibson'] = {'Tm':50, 'nM':500, 'len':30, 'window':40}
    params['goldengate'] = {'Tm':50,'nM':20,
                            'padding':'atatatatgg',
                            'ggsite':'GGTCTCn','ggN':4,
                            'window':20}
    params['LAMP'] = {'outer':{'Tm':60, 'nM':500, 'len':[20,40], 'gap':[0,20]},
                      'inner':{'Tm':62, 'nM':500, 'len':[20,40,60], 'gap':[0,20]},
                      'stem':{'Tm':65, 'nM':500, 'len':[20,40], 'gap':[0,50]}}
    params['parasail'] = {'gap_penalty':[1,1]}
    params['method'] = 'differential_evolution'
        
    def __init__(self, args=None):
        if self.params['verbose']:
            print(params['verbose'])
        
    def hairpin_primers(self):
        self.info = '''
        Primers for hairpin assembly
        '''
        # fixed domains
        df = ftpbase.read_csv(self.params['infile'])
        N = len(df)
        name = df['name'].values
        seq = df['seq'].values
        idx = df['idx'].values
        dN = [nupack.Domain(seq[i], name=['seq',i]) for i in range(N)]
        dI = [nupack.Domain(idx[i], name=['idx',i]) for i in range(N)]
 
        # variable domains
        [A,B,C,D] = self.params['dim']
        dA = [nupack.Domain('N'*A, name=['A',i]) for i in range(N)]
        dB = [nupack.Domain('N'*B, name=['B',i]) for i in range(N)]
        dC = [nupack.Domain('N'*C, name=['C',i]) for i in range(N)]
        dD = [nupack.Domain('N'*D, name=['D',i]) for i in range(N)]

        # strands
        s1 = [nupack.TargetStrand([dA[i], dD[i], dN[i]], name=[name[i]+'_s1a',i]) for i in range(N)]
        s2 = [nupack.TargetStrand([~dA[i], ~dB[i], ~dC[i], dI[i], dB[i]], name=[name[i]+'_s1b',i]) for i in range(N)]
        s3 = [nupack.TargetStrand([~dD[i], ~dA[i], ~dC[i]], name=[name[i]+'_s2',i]) for i in range(N)]
        s4 = [nupack.TargetStrand([~dD[i], ~dA[i], ~dB[i], ~dI[i], dC[i], dA[i], dD[i]], name=[name[i]+'_hp', i]) for i in range(N)]
        
        # target complexes
        Cs1 = [nupack.TargetComplex([s1[i]],
               structure='U'+str(D+A+dN[i].nt()),
               name=['Cs1',i]) for i in range(N)]
        Cs2 = [nupack.TargetComplex([s2[i]],
               structure='U'+str(A)+'D'+str(B)+'(U'+str(dI[i].nt()+C)+')',
               name=['Cs2', i]) for i in range(N)]
        Cs1s2 = [nupack.TargetComplex([s2[i], s1[i]],
                  structure='D'+str(A)+'(D'+str(B)+'(U'+str(dI[i].nt()+C)+')+)U'+str(D+dN[i].nt()),
                  name=['Cs1s2',i]) for i in range(N)]
        Cs3 = [nupack.TargetComplex([s3[i]],
               structure='U'+str(A+C+D),
               name=['Cs3',i]) for i in range(N)]
        Cs4 = [nupack.TargetComplex([s4[i]],
               structure='D'+str(D+A)+'(U'+str(B+C+dI[i].nt())+')',
               name=['Cs4',i]) for i in range(N)]
        Cs3s4 = [nupack.TargetComplex([s4[i], s3[i]],
               structure='U'+str(D+A+B+dI[i].nt())+'D'+str(D+C+A)+'(+)',
               name=['Cs3s4',i]) for i in range(N)]

        # reaction test tubes
        # s1 + s2 -> s1s2
        step0 = [nupack.TargetTube({Cs1[i]:500e-9, Cs2[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[Cs1s2[i]]), name=['step0',i]) for i in range(N)]
        # s1s2 -> s1s2
        step1 = [nupack.TargetTube({Cs1s2[i]:500e-9}, nupack.SetSpec(max_size=2), name=['step1',i]) for i in range(N)]
        # s3 + s4 -> s3s4
        step2 = [nupack.TargetTube({Cs3[i]:500e-9, Cs4[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[Cs3s4[i]]), name=['step2',i]) for i in range(N)]
        # s3s4 -> s3s4
        step3 = [nupack.TargetTube({Cs3s4[i]:500e-9}, nupack.SetSpec(max_size=2), name=['step3',i]) for i in range(N)]
        self.tubes = step0 + step1 + step2 + step3
     
        # crosstalk tube to prevent primer dimers
        x1 = {}
        ex1 = [] 
        x2 = {}
        ex2 = []
        for i in range(N):
            x1.update({Cs1[i]:500e-9, Cs2[i]:500e-9})
            ex1+=[Cs1s2[i]]
            x2.update({Cs3[i]:500e-9, Cs4[i]:500e-9})
            ex2+=[Cs3s4[i]]
        cross1 = nupack.TargetTube(x1, nupack.SetSpec(max_size=2, exclude=ex1), name='cross1')
        cross2 = nupack.TargetTube(x2, nupack.SetSpec(max_size=2, exclude=ex2), name='cross2')
        self.tubes+= [cross1, cross2]
     
        # weights 
        weights = nupack.Weights(self.tubes)
        for i in range(N):
            weights[dN[i], :, :, :] = 0
            weights[dI[i], :, :, :] = 0
            weights[dA[i], :, :, :] = 2
            weights[dC[i], :, :, :] = 2
        # weight primer dimer tube
        weights[:, :, :, cross2]*= N
        self.weights = weights     
        
        # constraints
        self.scopes = [[~dA[i], ~dB[i], ~dC[i]] for i in range(N)]
        self.scopes+= [[~dD[i], ~dA[i], ~dC[i]] for i in range(N)]

    def pcr_overlaps(self):
        self.info = '''
        Orthogonal PCR overlap design
        '''
        # variable domains
        [A,nA,B,nB] = self.params['dim']
        dA = [nupack.Domain('N'*A, name=['A',i]) for i in range(nA)]
        dB = [nupack.Domain('N'*B, name=['B',i]) for i in range(nB)]
        self.scopes = [[dA[i]] for i in range(nA)]
        self.scopes+= [[dB[i]] for i in range(nB)]

        # strands
        s1 = [nupack.TargetStrand([dA[i]], name=['sAF_'+str(i)]) for i in range(nA)]
        s2 = [nupack.TargetStrand([~dA[i]], name=['sAR_'+str(i)]) for i in range(nA)]
        s3 = [nupack.TargetStrand([dB[i]], name=['sBF_'+str(i)]) for i in range(nB)]
        s4 = [nupack.TargetStrand([~dB[i]], name=['sBR_'+str(i)]) for i in range(nB)]

        # target complexes
        c1 = [nupack.TargetComplex([s1[i]], structure='U'+str(dA[i].nt()), name=['c1',i]) for i in range(nA)]
        c2 = [nupack.TargetComplex([s2[i]], structure='U'+str(dA[i].nt()), name=['c2',i]) for i in range(nA)]
        c12 = [nupack.TargetComplex([s1[i], s2[i]], structure='D'+str(dA[i].nt())+'(+)', name=['c12',i]) for i in range(nA)]
        c3 = [nupack.TargetComplex([s3[i]], structure='U'+str(dB[i].nt()), name=['c3',i]) for i in range(nB)]
        c4 = [nupack.TargetComplex([s4[i]], structure='U'+str(dB[i].nt()), name=['c4',i]) for i in range(nB)]
        c34 = [nupack.TargetComplex([s3[i], s4[i]], structure='D'+str(dB[i].nt())+'(+)', name=['c34',i]) for i in range(nB)]
 
        # s1 + s2 -> s1s2
        step0 = [nupack.TargetTube({c1[i]:500e-9, c2[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[c12[i]]), name=['step0',i]) for i in range(nA)]
        step1 = [nupack.TargetTube({c3[i]:500e-9, c4[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[c34[i]]), name=['step1',i]) for i in range(nB)]
        # s1s2 -> s1s2
        self.tubes = step0 + step1
        # crosstalk
        tubes = []
        for i in range(nA):
            cross = [nupack.TargetTube({c1[i]:500e-9, c3[j]:500e-9}, nupack.SetSpec(max_size=2), name=['crossF_R',i,j]) for j in range(nB)]
            cross+= [nupack.TargetTube({c1[i]:500e-9, c4[j]:500e-9}, nupack.SetSpec(max_size=2), name=['crossF_Rc',i,j]) for j in range(nB)]
            tb = [nupack.TargetTube({c1[i]:500e-9, c2[j]:500e-9}, nupack.SetSpec(max_size=2), name=['crossF_Fc',i,j]) for j in range(nA)]
            del(tb[i])
            cross+= tb
            tubes+= cross

        for i in range(nB):
            cross = [nupack.TargetTube({c3[i]:500e-9, c2[j]:500e-9}, nupack.SetSpec(max_size=2), name=['crossR_Fc',i,j]) for j in range(nA)]
            tb = [nupack.TargetTube({c3[i]:500e-9, c4[j]:500e-9}, nupack.SetSpec(max_size=2), name=['crossR_Rc',i,j]) for j in range(nB)]
            del(tb[i])
            cross+= tb
            tubes+= cross
        self.tubes = tubes
        # weights 
        weights = nupack.Weights(self.tubes)
        for i in range(nA):
            weights[dA[i], :, :, :] = 1
        for i in range(nB):
            weights[dB[i], :, :, :] = 1
        # weight primer dimer tube
        #weights[:, :, :, cross1]*= N
        self.weights = weights

    def rna_trig(self):
        self.info = '''
        Orthogonal rna design
        '''
        # variable domains
        [A,nA,N] = self.params['dim']
        dA = [[nupack.Domain('N'*A, name=['A',i,j]) for i in range(nA)] for j in range(N)]
        L = 3
        dL = nupack.Domain('T'*L, name=['L'])
        self.scopes = [dA[j] for j in range(N)]

        # strands
        s1 = [nupack.TargetStrand(dA[j], name=['trig',j]) for j in range(N)]
        s2 = [[nupack.TargetStrand([dA[j][i+1], dA[j][i+2], dL, ~dA[j][i+2], ~dA[j][i+1], ~dA[j][i]], name=['hp5',i,j]) for i in range(nA-2)] for j in range(N)]
        s3 = [[nupack.TargetStrand([~dA[j][i+2], ~dA[j][i+1], ~dA[j][i], dL, dA[j][i], dA[j][i+1]], name=['hp3',i,j]) for i in range(nA-2)] for j in range(N)]

        # target complexes
        c1 = [nupack.TargetComplex([s1[j]], structure='U'+str(A*nA), name=['c1',j]) for j in range(N)]
        c2 = [[nupack.TargetComplex([s2[j][i]], structure='D'+str(A*2)+'(U'+str(L)+')U'+str(A), name=['c2',i,j]) for i in range(nA-2)] for j in range(N)]
        c3 = [[nupack.TargetComplex([s3[j][i]], structure='U'+str(A)+'D'+str(A*2)+'(U'+str(L)+')', name=['c3',i,j]) for i in range(nA-2)] for j in range(N)]
        c12 = [[nupack.TargetComplex([s2[j][i], s1[j]], structure='U'+str(A*2+L)+'D'+str(A*3)+'(+U'+str(A*i)+')U'+str(A*(nA-3-i)), name=['c12',i,j]) for i in range(nA-2)] for j in range(N)]
        c13 = [[nupack.TargetComplex([s3[j][i], s1[j]], structure='D'+str(A*3)+'(U'+str(L+A*2)+'+U'+str(A*i)+')U'+str(A*(nA-3-i)), name=['c13',i,j]) for i in range(nA-2)] for j in range(N)]

        self.tubes = [] 
        # crosstalk
        for i in range(N):
            # hp5 tubes on target
            tube5 = [nupack.TargetTube({c1[i]:500e-9, c2[i][j]:500e-9}, nupack.SetSpec(max_size=2, exclude=[c12[i][j]]), name=['step00',i,j]) for j in range(nA-2)]
            # hp3 tubes on target
            tube3 = [nupack.TargetTube({c1[i]:500e-9, c3[i][j]:500e-9}, nupack.SetSpec(max_size=2, exclude=[c13[i][j]]), name=['step01',i,j]) for j in range(nA-2)]
            self.tubes+= tube5 + tube3
            # off target tubes
            for k in range(N):
                if k!=i:
                    tube5 = [nupack.TargetTube({c1[i]:500e-9, c2[k][j]:500e-9}, nupack.SetSpec(max_size=2), name=['cross0',i,k,j]) for j in range(nA-2)]
                    tube3 = [nupack.TargetTube({c1[i]:500e-9, c3[k][j]:500e-9}, nupack.SetSpec(max_size=2), name=['cross1',i,k,j]) for j in range(nA-2)]
                    self.tubes+= tube5+tube3 

        # weights 
        weights = nupack.Weights(self.tubes)
        weights[dL, :, :, :] = 0
        self.weights = weights

    def golden_gate_primers(self):
        self.info = '''
        golden gate primers
        '''
        # fixed domains
        df = ftpbase.read_csv(self.params['infile'])
        N = len(df)
        name = df['name'].values
        seq = df['seq'].values
        idx = df['idx'].values
        dN = [nupack.Domain(seq[i], name=['seq',i]) for i in range(N)]
        dI = [nupack.Domain(idx[i], name=['idx',i]) for i in range(N)]
        
        # variable domains
        [A,B,C,D] = self.params['dim']
        dA = [nupack.Domain('N'*A, name=['A',i]) for i in range(N)]
        dB = [nupack.Domain('N'*B, name=['B',i]) for i in range(N)]
        dC = [nupack.Domain('N'*C, name=['C',i]) for i in range(N)]
        dD = [nupack.Domain('N'*D, name=['D',i]) for i in range(N)]
        dE = [nupack.Domain('N3', name=['E',i]) for i in range(N)]

        # strands
        s0a = [nupack.TargetStrand([dB[i], dA[i], dN[i]], name=['s0a_'+name[i], i]) for i in range(N)]
        s0b = [nupack.TargetStrand([dB[i], dA[i]], name=['s0b_'+name[i], i]) for i in range(N)]
        s1 = [nupack.TargetStrand([dB[i], dA[i], dN[i], dE[i], ~dN[i], ~dA[i]], name=['hp1_'+name[i], i]) for i in range(N)]
        s2 = [nupack.TargetStrand([~dB[i], ~dC[i], dI[i], dD[i], dC[i]], name=['s2_'+name[i], i]) for i in range(N)]
        s3 = [nupack.TargetStrand([~dA[i], ~dD[i]], name=['s3_'+name[i], i]) for i in range(N)]
        s4 = [nupack.TargetStrand([~dA[i], ~dB[i], ~dC[i], dI[i], dD[i], dA[i]], name=['hp4_'+name[i], i]) for i in range(N)]
        
        # target complexes
        Cs0a = [nupack.TargetComplex([s0a[i]],
               structure='U'+str(B+A+dN[i].nt()),
               name=['Cs0a',i]) for i in range(N)]
        Cs0b = [nupack.TargetComplex([s0b[i]],
               structure='U'+str(B+A),
               name=['Cs0b',i]) for i in range(N)]
        Cs1 = [nupack.TargetComplex([s1[i]],
               structure='U'+str(B)+'D'+str(A+dN[i].nt())+'(U'+str(dE[i].nt())+')',
               name=['Chp1',i]) for i in range(N)]
        Cs2 = [nupack.TargetComplex([s2[i]],
               structure='U'+str(B)+'D'+str(C)+'(U'+str(dI[i].nt()+D)+')',
               name=['Cs2', i]) for i in range(N)]
        Cs1s2 = [nupack.TargetComplex([s1[i], s2[i]],
               structure='D'+str(B)+'(D'+str(A+dN[i].nt())+'(U'+str(dE[i].nt())+')+)D'+str(C)+'(U'+str(dI[i].nt()+D)+')',
               name=['Cs1s2',i]) for i in range(N)]
        Cs3 = [nupack.TargetComplex([s3[i]],
               structure='U'+str(A+D),
               name=['Cs3',i]) for i in range(N)]
        Cs4 = [nupack.TargetComplex([s4[i]],
               structure='D'+str(A)+'(U'+str(B+C+dI[i].nt()+D)+')',
               name=['Chp4',i]) for i in range(N)]
        Cs3s4 = [nupack.TargetComplex([s4[i], s3[i]],
               structure='U'+str(A+B+C+dI[i].nt())+'D'+str(D+A)+'(+)',
               name=['Cs3s4',i]) for i in range(N)]

        # reaction test tubes
        self.tubes = []
        # s1 + s2 -> s1s2
        self.tubes+= [nupack.TargetTube({Cs1[i]:500e-9, Cs2[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[Cs1s2[i]]), name=['step1',i]) for i in range(N)]
        # s1s2 -> s1s2
        self.tubes+= [nupack.TargetTube({Cs1s2[i]:500e-9}, nupack.SetSpec(max_size=2), name=['step2',i]) for i in range(N)]
        # s3 + s4 -> s3s4
        self.tubes+= [nupack.TargetTube({Cs3[i]:500e-9, Cs4[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[Cs3s4[i]]), name=['step3',i]) for i in range(N)]
        # s3s4 -> s3s4
        self.tubes+= [nupack.TargetTube({Cs3s4[i]:500e-9}, nupack.SetSpec(max_size=2), name=['step4',i]) for i in range(N)]
     
        # crosstalk tube
        x0 = {}
        x1 = {}
        x2 = {}
        ex0 = []
        ex1 = []
        ex2 = []
        for i in range(N):
            x0.update({Cs0a[i]:500e-9, Cs0b[i]:500e-9})
            x1.update({Cs1[i]:500e-9, Cs2[i]:500e-9})
            ex1+=[Cs1s2[i]]
            x2.update({Cs3[i]:500e-9, Cs4[i]:500e-9})
            ex2+=[Cs3s4[i]]
        cross0 = nupack.TargetTube(x0, nupack.SetSpec(max_size=2), name='cross0')
        cross1 = nupack.TargetTube(x1, nupack.SetSpec(max_size=2, exclude=ex1), name='cross1')
        cross2 = nupack.TargetTube(x2, nupack.SetSpec(max_size=2, exclude=ex2), name='cross2')
        self.tubes+= [cross0, cross1, cross2]
     
        # weights 
        weights = nupack.Weights(self.tubes)
        for i in range(N):
            weights[dN[i], :, :, :] = 0
            weights[dB[i], :, :, :] = 2
            weights[dD[i], :, :, :] = 2
        # weight primer dimer tube
        weights[:, :, :, cross0]*= N
        weights[:, :, :, cross1]*= N
        weights[:, :, :, cross2]*= N
        self.weights = weights     
        
        # constraints
        self.scopes = [[~dA[i], ~dB[i], ~dC[i]] for i in range(N)]
        self.scopes+= [[dD[i], dA[i]] for i in range(N)]

    def malbec_primers(self):
        self.info = '''
        malbec primers
        '''
        df = ftpbase.read_csv(self.params['infile']) 
        # fixed domains
        name = df['name'].values[::2]
        seq1 = df['seq'].values[::2]
        seq2 = df['seq'].values[1::2]
        idxF = df['idx'].values[::2]
        idxR = df['idx'].values[1::2]

        N = len(seq1)
        FdN = [nupack.Domain(seq1[i], name=['seqF',i]) for i in range(N)]
        RdN = [nupack.Domain(seq2[i], name=['seqR',i]) for i in range(N)]
        FdI = [nupack.Domain(idxF[i], name=['idxF',i]) for i in range(N)]
        RdI = [nupack.Domain(idxR[i], name=['idxR',i]) for i in range(N)]
       
        # variable domains
        [A,B,C,D,E] = self.params['dim']
        if 'dA' in df.columns:
            seqA = df['dA'].values
            dA = [nupack.Domain(seqA[i], name=['A',i]) for i in range(N)]
        else:
            dA = [nupack.Domain('N'*A, name=['A',i]) for i in range(N)]
        dB = [nupack.Domain('N'*B, name=['B',i]) for i in range(N)]
        dC = [nupack.Domain('N'*C, name=['C',i]) for i in range(N)]
        dD = [nupack.Domain('N'*D, name=['D',i]) for i in range(N)]
        dE = [nupack.Domain('N'*E, name=['E',i]) for i in range(N)]

        # strands
        Fs1 = [nupack.TargetStrand([dB[i], dA[i], FdI[i], FdN[i]], name=[name[i]+'_s1F', i]) for i in range(N)]
        Rs1 = [nupack.TargetStrand([dB[i], dA[i], RdI[i], RdN[i]], name=[name[i]+'_s1R', i]) for i in range(N)]
        s1 = [nupack.TargetStrand([dB[i], dA[i]], name=[name[i]+'_s1', i]) for i in range(N)]
        s2 = [nupack.TargetStrand([dB[i], dA[i], FdI[i], FdN[i], ~RdN[i], ~RdI[i], ~dA[i]], name=[name[i]+'_hp2', i]) for i in range(N)]
        s3 = [nupack.TargetStrand([~dB[i], ~dC[i], ~dE[i], dD[i], dC[i]], name=[name[i]+'_s3', i]) for i in range(N)]
        Fs4 = [nupack.TargetStrand([~dA[i], ~dD[i]], name=[name[i]+'_s4F', i]) for i in range(N)]
        Rs4 = [nupack.TargetStrand([~dA[i], ~dE[i]], name=[name[i]+'_s4R', i]) for i in range(N)]
        Fs5 = [nupack.TargetStrand([~dA[i], ~dB[i], ~dC[i], ~dE[i], dD[i], dA[i]], name=[name[i]+'_hp5F', i]) for i in range(N)]
        Rs5 = [nupack.TargetStrand([~dA[i], ~dB[i], ~dC[i], ~dD[i], dE[i], dA[i]], name=[name[i]+'_hp5R', i]) for i in range(N)]
        
        # target complexes
        CFs1 = [nupack.TargetComplex([Fs1[i]],
               structure='U'+str(B+dA[i].nt()+FdI[i].nt()+FdN[i].nt()),
               name=['CFs1',i]) for i in range(N)]
        CRs1 = [nupack.TargetComplex([Rs1[i]],
               structure='U'+str(B+dA[i].nt()+RdI[i].nt()+RdN[i].nt()),
               name=['CRs1', i]) for i in range(N)]
        Cs1 = [nupack.TargetComplex([s1[i]],
               structure='U'+str(B+dA[i].nt()),
               name=['Cs1', i]) for i in range(N)]
        
        Cs2 = [nupack.TargetComplex([s2[i]],
               structure='U'+str(B)+'D'+str(dA[i].nt())+'(U'+str(FdI[i].nt()+FdN[i].nt()+RdN[i].nt()+RdI[i].nt())+')',
               name=['Cs2',i]) for i in range(N)]
        Cs3 = [nupack.TargetComplex([s3[i]],
               structure='U'+str(B)+'D'+str(C)+'(U'+str(D+E)+')',
               name=['Cs3',i]) for i in range(N)]
        Cs2s3 = [nupack.TargetComplex([s2[i], s3[i]],
               structure='D'+str(B)+'(D'+str(dA[i].nt())+'(U'+str(FdI[i].nt()+FdN[i].nt()+RdN[i].nt()+RdI[i].nt())+')+)D'+str(C)+'(U'+str(E+D)+')',
               name=['Cs2s3',i]) for i in range(N)]
        
        CFs4 = [nupack.TargetComplex([Fs4[i]],
               structure='U'+str(D+dA[i].nt()),
               name=['CFs4',i]) for i in range(N)]
        CRs4 = [nupack.TargetComplex([Rs4[i]],
               structure='U'+str(E+dA[i].nt()),
               name=['CRs4',i]) for i in range(N)]
        CFs5 = [nupack.TargetComplex([Fs5[i]],
               structure='D'+str(dA[i].nt())+'(U'+str(B+C+E+D)+')',
               name=['CFs5',i]) for i in range(N)]
        CRs5 = [nupack.TargetComplex([Rs5[i]],
               structure='D'+str(dA[i].nt())+'(U'+str(B+C+E+D)+')',
               name=['CRs5',i]) for i in range(N)]
        CFs4Fs5 = [nupack.TargetComplex([Fs5[i], Fs4[i]],
               structure='U'+str(dA[i].nt()+B+C+E)+'D'+str(D+dA[i].nt())+'(+)',
               name=['CFs4Fs5',i]) for i in range(N)]
        CRs4Rs5 = [nupack.TargetComplex([Rs5[i], Rs4[i]],
               structure='U'+str(dA[i].nt()+B+C+D)+'D'+str(E+dA[i].nt())+'(+)',
               name=['CRs4Rs5',i]) for i in range(N)]
        
        # reaction test tubes
        # s2 + s3 -> s2s3
        step1 = [nupack.TargetTube({Cs2[i]:500e-9, Cs3[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[Cs2s3[i]]), name=['step1',i]) for i in range(N)]
        step2 = [nupack.TargetTube({Cs2s3[i]:500e-9}, nupack.SetSpec(max_size=2), name=['step2',i]) for i in range(N)]
        # s4 + s5 -> s4s5
        step3 = [nupack.TargetTube({CFs4[i]:500e-9, CFs5[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[CFs4Fs5[i]]), name=['Fstep3',i]) for i in range(N)]
        step3+= [nupack.TargetTube({CRs4[i]:500e-9, CRs5[i]:500e-9}, nupack.SetSpec(max_size=2, exclude=[CRs4Rs5[i]]), name=['Rstep3',i]) for i in range(N)]
        step4 = [nupack.TargetTube({CFs4Fs5[i]:500e-9}, nupack.SetSpec(max_size=2), name=['Fstep4',i]) for i in range(N)]
        step4+= [nupack.TargetTube({CRs4Rs5[i]:500e-9}, nupack.SetSpec(max_size=2), name=['Rstep4',i]) for i in range(N)]
        self.tubes = step1 + step2 + step3 + step4
     
        # crosstalk tube
        x0 = {}
        x1 = {}
        x2 = {}
        ex0 = []
        ex1 = []
        ex2 = []
        for i in range(N):
            x0.update({CFs1[i]:500e-9, CRs1[i]:500e-9, Cs1[i]:500e-9})
            x1.update({Cs2[i]:500e-9, Cs3[i]:500e-9})
            ex1+=[Cs2s3[i]]
            x2.update({CFs4[i]:500e-9, CRs4[i]:500e-9})
 
        cross0 = nupack.TargetTube(x0, nupack.SetSpec(max_size=2), name='cross0')
        cross1 = nupack.TargetTube(x1, nupack.SetSpec(max_size=2, exclude=ex1), name='cross1')
        cross2 = nupack.TargetTube(x2, nupack.SetSpec(max_size=2), name='cross2')
        self.tubes+= [cross0, cross1, cross2]
     
        # weights 
        weights = nupack.Weights(self.tubes)
        for i in range(N):
            weights[dA[i], :, :, :] = 3
            weights[dB[i], :, :, :] = 2
            weights[dC[i], :, :, :] = 3
            weights[dD[i], :, :, :] = 2
            weights[dE[i], :, :, :] = 2
        # weight primer dimer tube
        weights[:, :, :, cross0]*= N
        weights[:, :, :, cross1]*= N
        weights[:, :, :, cross2]*= N
        self.weights = weights     
        
        # constraints
        self.scopes = [[~dA[i], ~dB[i], ~dC[i], ~dE[i], dD[i], dC[i], dB[i], dA[i]] for i in range(N)]

    def get_constraints(self):
        # define soft constraints
        self.soft = []
        self.hard = []
        # setting up enzyme constraints
        Rsites = [self.params['enzymes'][i] for i in self.params['enzymes'].keys()]
        Rsites+= [str(Bio.Seq.Seq(i).reverse_complement()) for i in Rsites]
        Rsites = np.array(Rsites)
        Rlen = np.array([len(i) for i in Rsites])
        if 'enzyme' in self.params['constraints']:
            logging.info('adding enzyme constraints')
            for j in range(len(self.scopes)):
                L = sum([d.nt() for d in self.scopes[j]])
                self.hard+= [nupack.Pattern(Rsites[Rlen <= L], scope=self.scopes[j])]
        if 'diversity' in self.params['constraints']:
            logging.info('adding diversity constraints')
            self.hard+= [nupack.Diversity(word=8, types=4, scope=self.scopes[j]) for j in range(len(self.scopes))]
        # hu6 poly T forbid
        if 'polyT' in self.params['constraints']:
            logging.info('adding polyT constraints')
            self.hard+= [nupack.Pattern(['T'*4], scope=self.scopes[j]) for j in range(len(self.scopes))]
            self.hard+= [nupack.Pattern(['A'*4], scope=self.scopes[j]) for j in range(len(self.scopes))]

        # apply window constraint
        if 'window_constraint' in self.params.keys():
            wlib = ftpbase.read_csv(self.params['window_constraint'])
            for j in range(len(self.scopes)):
                name, seq = wlib[['Strand','Sequence']].values[j]
                seq = seq.upper()
                trigRNA = self.scopes[j]
                N = np.sum([d.nt() for d in trigRNA])
                # check if sequence is amino acids
                try:
                    if self.isRNA(seq):
                        logging.info('RNA constraint for trig'+str(j)+str(t))
                        self.hard+=[nupack.Window(trigRNA, sources=[seq])]
                    elif self.isAminoAcid(seq):
                        logging.info('Amino Acid constraint for trig'+str(j))
                        self.hard+=[nupack.Library(trigRNA, self.get_AAlib(seq))]
                except:
                    logging.error('sequence constraint '+name+' is invalid')
                    logging.error('sequence = '+seq)
                    logging.error('trigger length = '+str(N))
                    logging.error('mRNA length = '+str(len(seq)))
                    logging.error(str(self.params))
                    sys.exit(1)
                    
    def run_nupack(self):
        import nupack
        # change block cache size
        nupack.config.cache = self.params['cache']
        # specify model parameters
        self.model = nupack.Model(material=self.params['material'], celsius=self.params['temperature'])
        self.options = nupack.DesignOptions(f_stop=self.params['fstop'], wobble_mutations=self.params['wobble'])
        # figure out which test tube to use 
        tube_func = getattr(self, self.params['tube_spec'])
        tube_func()
        # print design info
        if self.params['verbose']:
            logging.info(self.info)
            logging.info(str(self.params))
        # output filename
        outname = self.params['outfile'].split('.')[0]
        outname+='_L'+''.join([str(i)+'_' for i in self.params['dim']])
        outname = outname[:-1]
        self.get_constraints()
        d = nupack.Design(self.tubes, soft_constraints=self.soft, hard_constraints=self.hard, options=self.options, model=self.model)
        if self.params['pickle']:
            for t in range(self.params['trials']):
                fname = outname+'_trial'+str(t)+'.spec'
                logging.info('pickling '+fname)
                d.save(fname)


