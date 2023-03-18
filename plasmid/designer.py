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

# system libraries
import time

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
            print(params['verbose'])
        
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
            df.loc[i,'sequence'] = base.revcomp(seq)
            df.loc[i,'strand']*=-1
            df.loc[i,'locus_tag']+='c'
        return df

    def xtPCR(self, fL, seq, fR=None, padding=[2,2], niter=3, w=[10,100,1,1,2], verbose=False, get_cost=False):
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
        if type(fL).__name__ == 'plasmid': fL = fL.seq()
        if type(seq).__name__ == 'plasmid': seq = seq.seq()
        if type(fR).__name__ == 'plasmid': fR = fR.seq()

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
                if verbose:
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
            seqR = base.revcomp(fL+seq)
            fR = base.revcomp(fR)
            print('running fwd')
            df_F = self.xtPCR(fL, seq+fR, padding=padding, w=w, verbose=verbose, get_cost=get_cost)
            fwd = df_F[0]
            print('running rev')
            df_R = self.xtPCR(fR, seqR, padding=padding[::-1], w=w, verbose=verbose, get_cost=get_cost)
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
        seqR = base.revcomp(fL+seq)
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

    def iPCR(self, insert, seq, ggsite='', ggN=4, idx=None, w=[10,100,1,1,1,1], verbose=False):
        '''
        Designs iPCR primers that ligate via blunt or with golden gate
        insert = sequence to insert
        ggsite = golden gate enzyme restriction site
        ggN = number of bases for golden gate overlap
        seq = linearized fragment where new sequences will be appended to ends
        w = weights to the cost function
        '''
        if type(insert).__name__ == 'plasmid': insert = insert.seq()
        if type(seq).__name__ == 'plasmid': seq = seq.seq()
        if len(insert)==0: idx=0
 
        if idx==None:
            # find best primers by iterating through all possibilities
            best_cost = np.inf
            for idx in range(len(insert)):
                y = self.iPCR_cost(idx, insert, seq, ggsite, ggN, w, get_cost=False, verbose=verbose)
                if y[1] < best_cost:
                    best_cost = y[1]
                    best_oligos = y[0]
                print('processing index at', idx)
                print('[best, current] = ', str([best_cost, y[1]]))
        else:
            y = self.iPCR_cost(idx, insert, seq, ggsite, ggN, w, get_cost=False, verbose=verbose)
            best_oligos = y[0]
        return best_oligos

    def iPCR_cost(self, idx, insert, seq, ggsite, ggN, w, get_cost=True, verbose=False):
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
                fR = insert[:idx] + base.revcomp(ggsite)
            else:
                fF = ggsite + seq[-ggN:] + insert
                fR = base.revcomp(ggsite)
        # blunt ligation
        else:
            fF = insert[idx:]
            fR = insert[:idx]
        # run xtPCR
        y = self.xtPCR(fF, seq, fR, w=w[:-1], get_cost=True, verbose=verbose)
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
        verbose = print out assembled construct and highlight plasmid locations
        returns list of primers
        '''
        # format input
        frags = []
        pLseq = ''
        for i in range(len(seqlist)):
            f = seqlist[i]
            # convert from plasmid type
            if type(f).__name__=='plasmid':
                frags.append(['',f.seq(),''])
                pLseq = pLseq + f
            # is a list of [fL, seq, fR]
            elif type(f) == list:
                x = []
                for j in f:
                    if type(j).__name__ == 'plasmid':
                        x.append(j.seq())
                    else:
                        x.append(j)
                    pLseq = pLseq + j
                frags.append(x)
            else:
                frags.append(['',f,''])

        # check if fully assembled construct has golden gate site
        pLseq = pLseq.annotate(label='ggsite', sequence=self.params['goldengate']['ggsite']) 
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
        siteR = base.revcomp(site)
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
        seqR = [base.revcomp(i) for i in seq]
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
        if type(seq).__name__ == 'plasmid':
            seq = seq.seq()

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
        args = (Lb, Ub, conc, tTm, seq, base.revcomp(seq), w, False)
        bounds = [(Lb[i], Ub[i]) for i in range(len(Lb))]
        x0 = Lb
        res = self.optimizer(self.LAMP_cost, args, bounds, x0)
        # status output
        cost = self.LAMP_cost(res.x, Lb, Ub, conc, tTm, seq, base.revcomp(seq), w, verbose=True)
        if verbose:
            print('res.x:', res.x)
            print('res.fun', res.fun)
            print('cost:', cost)
            print('elapse', time.time()-start) 
        
        # get forward primers
        fwd = self.LAMP_results(res.x[:6], seq, conc)
        fwd['locus_tag'] = [i+'_F' for i in fwd['locus_tag']]
        # get reverse primers
        rev_seq = base.revcomp(seq)
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
        stem = base.revcomp(oligos[2])
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
        seqR = [base.revcomp(i) for i in seq]
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
            seq = seq.annotate(label=label, sequence=sequence)
        return seq

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


