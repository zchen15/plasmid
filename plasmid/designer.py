#!/usr/bin/env python
# Class functions to design primers

# computing libraries
import numpy as np
import scipy as sp
import scipy.optimize
import pandas as pd

# reading and writing genbank files
import Bio
import Bio.SeqUtils
import Bio.SeqUtils.MeltingTemp
import parasail

# system libraries
import time

from .misc import revcomp
from .misc import json_pretty
from .misc import codon_table
from .misc import translate
from .misc import synonymous_codon_table

class Designer:
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
            out = json_pretty(self.params)
            print(out)
    
    def optimizer(self, cost_func, args, bounds, x0, discrete=False):
        cons = ()
        method = self.params['method']
        if discrete:
            cons = (sp.optimize.NonlinearConstraint(self.integer_deviation,-0.1,0.1))
        if method=='shgo':
            kwargs = {'method':'Powell'}
            return sp.optimize.shgo(cost_func, bounds=bounds, args=args, minimizer_kwargs=kwargs)
        elif method=='dual_annealing':
            kwargs = {'method':'Powell'}
            return sp.optimize.dual_annealing(cost_func, bounds=bounds, args=args, local_search_options=kwargs)
        elif method=='basinhopping':
            kwargs = {'args':args, 'method':'Powell', 'bounds':bounds}
            return sp.optimize.basinhopping(cost_func, x0, minimizer_kwargs=kwargs)
        elif method=='differential_evolution':
            return sp.optimize.differential_evolution(cost_func, bounds=bounds, args=args, constraints=cons)
    
    def integer_deviation(self, x):
        '''
        Calculate deviation from integer values
        '''
        return np.sum([abs(i-int(i)) for i in x])

    def get_Tm_NN(self, seq, conc, table=Bio.SeqUtils.MeltingTemp.DNA_NN4, Na=None):
        '''
        Use nearest neighbor model for Tm calculation
        seq = DNA sequence to get melting temperature
        Na = concentration of Na+, default = 50 or 50mM
        conc = primer concentration in nM, default = 500nM
        table = parameter table to use
                see help(Bio.SeqUtils.MeltingTemp)
        return float value of Tm in Celsius
        '''
        if Na==None:
            Na = self.params['xtPCR']['Na']
        return Bio.SeqUtils.MeltingTemp.Tm_NN(seq, dnac1=conc)

    def primer_alternate(self, df):
        '''
        Alternate direction of primers
        df = dataframe with columns [label, sequence, Tm, strand, annealed]
        returns dataframe with alternating fwd and reverse sequences
        '''
        # only work on sets with more than one primer
        if len(df) == 1:
            return df
        for i in range(len(df)-2,0,-2):
            seq = df.iloc[i]['sequence']
            df.loc[i,'sequence'] = revcomp(seq)
            df.loc[i,'strand']*=-1
            df.loc[i,'locus_tag']+='c'
        return df

    def xtPCR(self, fL, seq, fR=None, padding=[2,2], niter=3, w=[10,100,1,1,2], get_cost=False):
        '''
        Find primers which can seed and extend a PCR fragment
        fL = flanking sequence on 5' end
        seq = sequence on 3' end which gets amplified
        fR = flanking sequence on 3' end
        padding = number of extra primers to try
        w = weights for cost function
        method = optimization method
        returns list of primers
        '''
        # handle alternative datatypes safely
        fL = str(fL)
        seq = str(seq)
        
        if fR==None:
            if len(fL)==0:
                padding[0] = 1
            # find best primer set
            N_L = np.max([np.ceil(len(fL)/(self.params['xtPCR']['len'][1] - self.params['xtPCR']['len'][0])), 1]).astype(int)
            best_res = [N_L, None, np.inf] 
            start = time.time()                
            for i in range(padding[0]):
                N_left = N_L + i
                args = (fL, seq, w)
                bounds = [(self.params['xtPCR']['len'][0], self.params['xtPCR']['len'][1])]*N_left + [(0, self.params['xtPCR']['len'][1] - self.params['xtPCR']['len'][0])]*N_left
                x0 = [self.params['xtPCR']['len'][0]]*N_left + [0]*N_left
                res = self.optimizer(self.xtPCR_cost, args, bounds, x0) 
                # save the best result
                if best_res[2] > res.fun:
                    best_res = [N_left, res, res.fun]                    
                # status output
                if self.params['verbose']:
                    cost = self.xtPCR_cost(res.x, fL, seq, w, verbose=True)
                    print('N', N_left)
                    print('res.x:', res.x)
                    print('res.fun', res.fun)
                    print('cost:', cost)
                    print('elapse', time.time()-start) 
            # get strands
            N_left, res, cost = best_res
            oligos = self.xtPCR_get_oligo(res.x, fL, seq)
            col = ['sequence','annealed']
            oligos = pd.DataFrame(oligos, columns=col)
            annealed = oligos['annealed'].values
            conc = [self.params['xtPCR']['nM'][0]]* (len(annealed)-1) + [self.params['xtPCR']['nM'][1]]
            Tm = [self.get_Tm_NN(annealed[i], conc=conc[i]) for i in range(len(annealed))]
            oligos['Tm'] = Tm
            if get_cost:
                cost = self.xtPCR_cost(res.x, fL, seq, w, verbose=True)
                return [oligos, cost]
            return [oligos]
        else:
            seqR = revcomp(fL+seq)
            fR = revcomp(fR)
            print('running fwd')
            df_F = self.xtPCR(fL, seq+fR, padding=padding, w=w, get_cost=get_cost)
            fwd = df_F[0]
            print('running rev')
            df_R = self.xtPCR(fR, seqR, padding=padding[::-1], w=w, get_cost=get_cost)
            rev = df_R[0]
            fwd['locus_tag'] = [str(i)+'_F' for i in range(len(fwd)-1)] + ['fin_F']
            fwd['strand'] = 1
            rev['locus_tag'] = [str(i)+'_R' for i in range(len(rev)-1)] + ['fin_R']
            rev['strand'] = -1
            oligos = pd.concat([fwd, rev])
            col = ['locus_tag','Tm','sequence','annealed','strand']
            oligos = oligos[col]
            if get_cost:
                cost = df_F[1] + df_R[1]
                return [oligos, cost]
            return oligos
    
    def xtPCR_cost(self, pos, fL, seq, w, verbose=False):
        '''
        Compute oligo cost
        pos = oligo positions
        fL = left extension
        seq = seed sequence
        w = weights for each cost function
        '''
        cost = np.zeros(len(w))
        # Deviation from min and max primer length
        N = int(len(pos)/2)
        Lb = [self.params['xtPCR']['len'][0]]*N + [0]*N
        c1 = Lb - pos 
        p3 = pos[:N].astype(int) + pos[N:].astype(int)
        c2 = p3 - self.params['xtPCR']['len'][1]
        c1[c1 < 0] = 0
        c2[c2 < 0] = 0
        cost[0] = np.sum(c1) + np.sum(c2)
        # amount of sequence covered
        ucL = len(fL) - np.sum(pos[N:].astype(int))
        cost[1] =  np.max([ucL, 0])
        # total primer length
        cost[2] = np.sum(pos)
        # Tm deviation from target temperature
        oligos = self.xtPCR_get_oligo(pos, fL, seq)
        Tm = [self.get_Tm_NN(i[1], conc=self.params['xtPCR']['nM'][0]) for i in oligos[:-1]]
        Tm+= [self.get_Tm_NN(oligos[-1][1], conc=self.params['xtPCR']['nM'][1])]
        Tm = np.array(Tm)
        cost[3] = np.sum((self.params['xtPCR']['Tm'] - Tm)**2)
        # cross reactivity
        if w[4]!=0:
            cost[4] = self.xtPCR_crosscheck(pos, fL, seq)
        if verbose:
            return cost
        return np.sum(w*cost)

    def xtPCR_get_oligo(self, pos, fL, seq, verbose=False):
        '''
        Get oligo list
        pos = primer positions
        fL = left extension
        seq = seed sequence
        '''
        N = int(len(pos)/2)
        x1 = pos[:N].astype(int)
        x2 = pos[N:].astype(int)
        out = []
        p1 = seq
        p2 = fL
        for i in range(N):
            s1 = x1[i]
            s2 = x2[i]
            if s1 > len(p1):
                s1 = len(p1)
            if s2 > 0:
                annealed = p1[:s1]
                xt = p2[-s2:]
                out.append([xt+' '+annealed, annealed])
                p1 = xt + p1
                p2 = p2[:-s2]
            else:
                annealed = p1[:s1]
                out.append([annealed, annealed])
            # terminate when no sequences are left extend 
            if len(p2) == 0:
                break
        return out

    def xtPCR_crosscheck(self, pos, fL, seq):
        '''
        Check cross reactivity of xtPCR primers
        '''
        g1, g2 = self.params['parasail']['gap_penalty']
        seqR = revcomp(fL+seq)
        N = int(len(pos)/2)
        x1 = pos[:N].astype(int)
        x2 = pos[N:].astype(int)
        p1 = seq
        p2 = fL
        score = []
        for i in range(N):
            s1 = x1[i]
            s2 = x2[i]
            # cut off out of bound lengths
            if s1 > len(p1):
                s1 = len(p1)
            
            exclude = [p2, p1[s1:], seqR]
            annealed = p1[:s1]
            if s2 > 0:
                p1 = p2[-s2:] + p1
                p2 = p2[:-s2]
            # score cross reactivity
            res1 = parasail.sg_scan_16(annealed, annealed, g1, g2, parasail.nuc44)
            score.append(res1.score - self.get_reactivity(annealed, exclude, res1.score))
        
            # terminate all sequences are covered
            if len(p2)==0:
                break

        # exact matches are highly penalized
        score = min(score)
        if score==0:
            return np.inf
        else:
            return -score

    def iPCR(self, insert, seq, ggsite='', ggN=4, idx=None, w=[10,100,1,1,1,1]):
        '''
        Designs iPCR primers that ligate via blunt or with golden gate
        insert = sequence to insert
        ggsite = golden gate enzyme restriction site
        ggN = number of bases for golden gate overlap
        seq = linearized fragment where new sequences will be appended to ends
        w = weights to the cost function
        '''
        if type(insert).__name__ == 'Plasmid': insert = str(insert)
        if type(seq).__name__ == 'Plasmid': seq = str(seq)
        if len(insert)==0: idx=0
 
        if idx==None:
            # find best primers by iterating through all possibilities
            best_cost = np.inf
            for idx in range(len(insert)):
                y = self.iPCR_cost(idx, insert, seq, ggsite, ggN, w, get_cost=False)
                if y[1] < best_cost:
                    best_cost = y[1]
                    best_oligos = y[0]
                print('processing index at', idx)
                print('[best, current] = ', str([best_cost, y[1]]))
        else:
            y = self.iPCR_cost(idx, insert, seq, ggsite, ggN, w, get_cost=False)
            best_oligos = y[0]
        return best_oligos

    def iPCR_cost(self, idx, insert, seq, ggsite, ggN, w, get_cost=True):
        '''
        Compute cost of iPCR
        '''
        # get primers for deletion
        if len(insert) == 0:
            insert = ' '
        
        # add golden site
        if ggsite!='':
            if idx-ggN > -1:
                fF = ggsite + insert[idx-ggN:]
                fR = insert[:idx] + revcomp(ggsite)
            else:
                fF = ggsite + seq[-ggN:] + insert
                fR = revcomp(ggsite)
        
        # blunt ligation
        else:
            fF = insert[idx:]
            fR = insert[:idx]
        
        # run xtPCR
        y = self.xtPCR(fF, seq, fR, w=w[:-1], get_cost=True)
        
        # weight evenness of split
        cost = [i for i in y[1]]
        cost+= [abs(len(insert)/2 - idx)]
        cost = np.array(cost)
        if get_cost:
            return np.sum(w*cost)
        else:
            return [y[0], np.sum(w*cost)]

    def GoldenGate(self, seqlist, exclude=[], w=[0,1], circular=True):
        '''
        Design primers for goldengate assembly
        seqlist = list of sequences to assemble
        exclude = sites to exclude
        circular = assemble fragments into a circular construct
        returns list of primers
        '''
        # format input
        frags = []
        pLseq = ''
        for i in range(len(seqlist)):
            f = seqlist[i]
            # convert from plasmid type
            if type(f).__name__=='Plasmid':
                frags.append(['',str(f),''])
                pLseq = pLseq + f
            # is a list of [fL, seq, fR]
            elif type(f) == list:
                x = []
                for j in f:
                    x.append(str(j))
                    pLseq = pLseq + j
                frags.append(x)
            else:
                frags.append(['',f,''])

        # check if fully assembled construct has golden gate site
        pLseq = pLseq.annotate(name='ggsite', sequence=self.params['goldengate']['ggsite']) 
        ggloc = ['ggsite' in L for L in pLseq['locus_tag']]
        if np.sum(ggloc) > 0:
            print('Goldengate site found in assembled construct')
            print(pLseq[['locus_tag','location']][ggloc])
        
        # generate paired frags
        overlaps = []
        L = self.params['goldengate']['window']
        for i in range(len(frags)-1):
            f1 = ''.join(frags[i])
            f2 = ''.join(frags[i+1])
            overlaps.append(f1[-L:]+f2[:L])
        if circular:
            f1 = ''.join(frags[-1])
            f2 = ''.join(frags[0])
            overlaps.append(f1[-L:]+f2[:L])
        
        # optimize overlaps
        start = time.time()
        Lb = [0]*len(overlaps)
        Ub = [len(i) - self.params['goldengate']['ggN'] for i in overlaps]
        args = (overlaps, exclude, w, False)
        bounds = [(Lb[i], Ub[i]) for i in range(len(Lb))]
        x0 = Lb
        res = self.optimizer(self.goldengate_cost, args, bounds, x0)
        
        # results of overlap opt
        print('res.x', res.x)
        print('res.fun', res.fun)
        print('exclude:',exclude)
        olap = self.goldengate_get_oligo(res.x, overlaps)
        Tm_olap = [self.get_Tm_NN(olap[i], conc=self.params['goldengate']['nM']) for i in range(len(olap))]
        print('overlaps:', olap)
        print('Tm overlap:', Tm_olap)
        
        # generate seed primers
        site = self.params['goldengate']['padding'] + self.params['goldengate']['ggsite']
        siteR = revcomp(site)
        ggN = self.params['goldengate']['ggN']
        L = self.params['goldengate']['window']
        y = res.x.astype(int)
        for i in range(len(frags)-1):
            # trim left side extension or main sequence if needed
            offset = L - y[i] - ggN - len(frags[i][2])
            if offset >= 0:
                frags[i][1] = frags[i][1][:len(frags[i][1]) - offset]
                frags[i][2] = siteR
            elif y[i] + ggN < L:
                frags[i][2] = frags[i][2][:len(frags[i][2]) - L + y[i] + ggN] + siteR
            else:
                frags[i][2] = frags[i][2] + overlaps[i][L:y[i]+ggN] + siteR
            # trim right side extension or main sequence if needed 
            offset = y[i] - L - len(frags[i+1][0]) 
            if offset >= 0:
                frags[i+1][1] = frags[i+1][1][offset:]
                frags[i+1][0] = site
            elif y[i] > L:
                frags[i+1][0] = site + frags[i+1][0][y[i]-L:]
            else:
                frags[i+1][0] = site + overlaps[i][y[i]:L] + frags[i+1][0]
                
        # process first and last
        if circular:
            # trim left side extension or main sequence if needed
            offset = L - y[-1] - ggN - len(frags[-1][2])
            if offset >= 0:
                frags[-1][1] = frags[-1][1][:len(frags[-1][1]) - offset]
                frags[-1][2] = siteR
            elif y[-1] + ggN < L:
                frags[-1][2] = frags[-1][2][:len(frags[-1][2]) - L + y[-1] + ggN] + siteR
            else:
                frags[-1][2] = frags[-1][2] + overlaps[-1][L:y[-1]+ggN] + siteR
            # trim right side extension or main sequence if needed 
            offset = y[-1] - L - len(frags[0][0]) 
            if offset >= 0:
                frags[0][1] = frags[0][1][offset:]
                frags[0][0] = site
            elif y[-1] > L:
                frags[0][0] = site + frags[0][0][y[-1]-L:]
            else:
                frags[0][0] = site + overlaps[-1][y[-1]:L] + frags[0][0]
        
        # run optimization for xtPCR
        out = []
        for i in range(len(frags)):
            print('processing primers for frag',i)
            [fL, seq, fR] = frags[i]
            p = self.xtPCR(fL, seq, fR)
            p['locus_tag'] = ['frag'+str(i)+'_'+j for j in p['locus_tag']]
            out.append(p)
        out = pd.concat(out)
        # add full length frags
        x = []
        for i in range(len(frags)):
            seq = frags[i][0] + frags[i][1] + frags[i][2]
            seq = str(seq)
            x.append(['seq'+str(i), seq])
        x = pd.DataFrame(x, columns=['locus_tag','sequence'])
        out = pd.concat([out,x])
        return out.sort_values(by='locus_tag').reset_index(drop=True)

    def goldengate_cost(self, x, frags, exclude, w, verbose):
        '''
        Compute design costs for golden gate
        '''
        cost = np.zeros(len(w))
        oligos = self.goldengate_get_oligo(x, frags)
        # deviation from target melting temperature      
        Tm = np.array([self.get_Tm_NN(oligos[i], conc=self.params['goldengate']['nM']) for i in range(len(oligos))])
        cost[0] = np.sum((self.params['goldengate']['Tm'] - Tm)**2)
        # check cross reactivity of overlaps
        cost[1] = self.overlap_check(oligos, exclude=exclude)
        if verbose:
            return cost
        return np.sum(w*cost)

    def goldengate_get_oligo(self, x, frags):
        '''
        Translate position to oligos
        '''
        x = x.astype(int)
        oligos = []
        for i in range(len(frags)):
            oligos.append(frags[i][x[i]:x[i] + self.params['goldengate']['ggN']])
        return oligos

    def overlap_check(self, seq, exclude=[]):
        '''
        Check if overlap sequences bind to each other
        '''
        g1, g2 = self.params['parasail']['gap_penalty']
        seqR = [revcomp(i) for i in seq]
        score = np.inf
        for i in range(len(seq)-1):
            res1 = parasail.sg_scan_16(seq[i], seq[i], g1, g2, parasail.nuc44)
            test = [seq[i+1:], seqR[i:], exclude]
            for tseq in test:
                score = np.min([score, res1.score - self.get_reactivity(seq[i], tseq, res1.score)])
                # premature exit
                if score==0:
                    return np.inf
        return -score

    def get_reactivity(self, seq, seqlist, best=None):
        score = -np.inf
        g1, g2 = self.params['parasail']['gap_penalty']
        for j in range(len(seqlist)):
            res2 = parasail.sg_scan_16(seq, seqlist[j], g1, g2, parasail.nuc44)
            score = np.max([score, res2.score])
            if best!=None and best==res2.score:
                return best
        return score

    def Gibson(self, seqlist, w=[10,1], method='differential_evolution', circular=True):
        '''
        Design primers for gibson assembly
        seqlist = list of sequences to assemble via gibson in order 
        circular = assemble fragments into a circular construct
        returns list of primers
        '''
        # use goldengate assembler
        self.params['goldengate']['Tm'] = self.params['gibson']['Tm']
        self.params['goldengate']['ggN'] = self.params['gibson']['len']
        self.params['goldengate']['window'] = self.params['gibson']['window']
        self.params['goldengate']['ggsite'] = ''
        self.params['goldengate']['padding'] = ''
        return self.GoldenGate(seqlist, circular=circular)

    def LAMP(self, seq, w=[100,1,1], verbose=False):
        '''
        Find optimal annealing sites for lamp primers
        seq = sequence to amplify
        fwd_only = get primers only in fwd orientation
        verbose = show optimization results
        '''
        # format input sequence
        seq = str(seq)

        # find best primer set
        start = time.time()                
        Lb = [self.params['LAMP']['outer']['gap'][0],
              self.params['LAMP']['outer']['len'][0], 
              self.params['LAMP']['inner']['gap'][0],
              self.params['LAMP']['inner']['len'][0],
              self.params['LAMP']['stem']['gap'][0],
              self.params['LAMP']['stem']['len'][0]]*2
        Ub = [self.params['LAMP']['outer']['gap'][1],
              self.params['LAMP']['outer']['len'][1],
              self.params['LAMP']['inner']['gap'][1],
              self.params['LAMP']['inner']['len'][1],
              self.params['LAMP']['stem']['gap'][1],
              self.params['LAMP']['stem']['len'][1]]*2
        conc = [self.params['LAMP']['outer']['nM'],
                self.params['LAMP']['inner']['nM'],
                self.params['LAMP']['stem']['nM']]*2
        tTm = [self.params['LAMP']['outer']['Tm'],
               self.params['LAMP']['inner']['Tm'],
               self.params['LAMP']['stem']['Tm']]*2
        args = (Lb, Ub, conc, tTm, seq, revcomp(seq), w, False)
        bounds = [(Lb[i], Ub[i]) for i in range(len(Lb))]
        x0 = Lb
        res = self.optimizer(self.LAMP_cost, args, bounds, x0)
        # status output
        cost = self.LAMP_cost(res.x, Lb, Ub, conc, tTm, seq, revcomp(seq), w, verbose=True)
        if verbose:
            print('res.x:', res.x)
            print('res.fun', res.fun)
            print('cost:', cost)
            print('elapse', time.time()-start) 
        
        # get forward primers
        fwd = self.LAMP_results(res.x[:6], seq, conc)
        fwd['locus_tag'] = [i+'_F' for i in fwd['locus_tag']]
        # get reverse primers
        rev_seq = revcomp(seq)
        rev = self.LAMP_results(res.x[6:], rev_seq, conc)
        rev['locus_tag'] = [i+'_R' for i in rev['locus_tag']]
        out = pd.concat([fwd, rev])
        out = out.sort_values(by=['locus_tag'])
        return out
    
    def LAMP_results(self, pos, seq, conc):
        '''
        Format primer results
        '''
        # get strands
        oligos = self.LAMP_get_oligo(pos, seq)
        out = pd.DataFrame()
        out['locus_tag'] = ['outer','inner','hairpin']
        # get stem
        stem = revcomp(oligos[2])
        # get hairpin sequence
        out['sequence'] = [oligos[0], stem+' '+oligos[1], stem+' '+oligos[3]+' '+oligos[2]]
        out['annealed'] = oligos[:-1]
        # get Tm
        Tm = [self.get_Tm_NN(oligos[i], conc=conc[i]) for i in range(len(out))]
        out['Tm'] = Tm
        return out

    def LAMP_cost(self, pos, Lb, Ub, conc, tTm, seq, seqR, w, verbose=False):
        '''
        Find optimal annealing sites for lamp primers
        pos = primer positions
        seq = sequence to amplify
        w = weights on the cost function
        '''
        cost = np.zeros(len(w))
        # Deviation from min and max primer length bounds
        c1 = Lb - pos
        c2 = pos - Ub
        c1[c1 < 0] = 0
        c2[c2 < 0] = 0
        Lmax = self.params['LAMP']['inner']['len'][2]
        c3 = max(pos[3] + pos[5] - Lmax, 0) + max(pos[9] + pos[11] - Lmax, 0)
        cost[0] = np.sum(c1) + np.sum(c2) + c3
        # Tm deviation from target temperature
        oligosF = self.LAMP_get_oligo(pos[:6], seq)
        oligosR = self.LAMP_get_oligo(pos[6:], seqR)
        oligos = oligosF[:-1] + oligosR[:-1]
        Tm = np.array([self.get_Tm_NN(oligos[i], conc=conc[i]) for i in range(len(oligos))])
        cost[1] = np.sum((tTm-Tm)**2)
        # cross reactivity
        if w[2]!=0:
            oligos = [oligosF[0], oligosR[0], oligosF[3], oligosR[3], oligosF[2],  oligosR[2]]
            cost[2] = self.LAMP_crosscheck(oligos)
        if verbose:
            return cost
        return np.sum(w*cost)

    def LAMP_get_oligo(self, pos, seq):
        '''
        Returns oligos associated with positions
        '''
        x = pos.astype(int)
        x[x < 0] = 0
        outer = seq[x[0]:x[0]+x[1]]
        s = x[0]+x[1]+x[2]
        inner = seq[s:s+x[3]]
        loop = seq[s:s+x[3]+x[4]]
        s+= x[3]+x[4]
        stem = seq[s:s+x[5]]
        return [outer, inner, stem, loop]

    def LAMP_crosscheck(self, seq, exclude=[]):
        '''
        Check cross reactivity of LAMP primers
        '''
        g1, g2 = self.params['parasail']['gap_penalty']
        seqR = [revcomp(i) for i in seq]
        query = seq[:-2] + seqR[-2:]
        ref = seqR[:-2] + seq[-2:]
        score = np.inf
        for i in range(len(query)):
            res1 = parasail.sg_scan_16(query[i], query[i], g1, g2, parasail.nuc44)
            test = [ref[i:], exclude]
            for tseq in test:
                score = np.min([score, res1.score - self.get_reactivity(query[i], tseq, res1.score)])
                # premature exit
                if score==0:
                    return np.inf
        return -score

    def overlay_primers(self, primers, seq, mask=''):
        '''
        overlay primers onto plasmid
        primers = dataframe with [label, sequence]
        seq = plasmid object
        returns primers overlay on plasmid
        '''
        primers = [[L, s.replace(' ','').replace(mask,'')] for L, s in primers[['locus_tag','sequence']].values]
        for label, sequence in primers:
            seq = seq.annotate(name=label, sequence=sequence)
        return seq

    def domesticate_CDS(seq, enzymes, max_opt=100):
        '''
        Change codons to remove restriction enzyme sites from a 
        coding gene while preserving the amino acid sequence
        seq = plasmid dataframe
        enzymes = list of restriction enzyme sequences
        max_opt = number of times to try optimization
        return plasmid dataframe without restriction enzymes
        '''
        # check which orientation it is
        strand = 1
        z = seq[seq['type']=='CDS']['strand'].values
        if len(z) > 0 and z[0]==-1:
            strand = z[0]
            seq = seq.reverse_complement()
        
        while max_opt > 0:
            print('max_opt=',max_opt)
            df = seq.copy()
            # mark restriction sites
            for eseq in enzymes:
                df = df.annotate('enzyme', eseq)

            # iterate through enzyme sites
            out = seq.copy()
            df = df[df['locus_tag']=='enzyme']
            for s1, s2 in df[['start','end']].values:
                s1 = int(np.floor(s1/3)*3)
                s2 = int(np.ceil(s2/3)*3)
                x = df.__str__()[s1:s2]
                mut = Designer.mutate_codon(x)

                # make the mutation
                s1 = mut['start']+s1
                s2 = s1+3
                key = range(s1,s2)
                val = mut['sequence']
                out.replace_sequence(key, val, disrupt=False, inplace=True)
                name = 'domesticate_CDS'
                out.annotate_by_location(name, s1, s2, 1, inplace=True)

            # check if enzyme sites are still there
            for eseq in enzymes:
                out = out.annotate('enzyme', eseq)
            if sum(out['locus_tag']=='enzyme') == 0:
                # restore original orientation if necessary
                if strand==-1:
                    out = out.reverse_complement()
                return out
            
            # repeat if optimization not successful up to max_opt
            max_opt-=1
        raise ValueError('Failed to domesticate CDS sequence')
        
    def mutate_codon(seq):
        '''
        Change codons for given amino acid sequence
        seq = input nucleotide sequence
        return 
        '''
        # select a codon at random from the sequence
        idx = np.random.randint(int(len(seq)/3))
        s1 = idx*3
        x = seq[s1:s1+3]
        
        # look up codon in synonymous codon table
        y = synonymous_codon_table[x]
        idx = np.random.randint(len(y))
        mut = y[idx]
        out = {'start':s1,
               'stop':s1+3,
               'sequence':mut}
        return out