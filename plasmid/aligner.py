#!/usr/bin/env python
# Class function for aligners

# computing libraries
import numpy as np
import pandas as pd

# reading and writing genbank files
import Bio
import Bio.SeqIO
import Bio.SeqUtils
import spoa
import parasail
import mappy

# system and time
import re

# import custom libraries
from .misc import read_to_df

class Aligner:
    '''
    This class holds functions pertaining to sequence alignment and analysis
    '''
    params = {'verbose':False}
    params['parasail'] = {'g1':5,
                          'g2':1,
                          'mat':parasail.nuc44,
                          'method':'sg_trace_scan_16'}
    params['minimap'] = {'k':7,
                         'w':1}

    def __init__(self, args=None):
        if args!=None:
            self.params = load_args(args, self.params)
        if self.params['verbose']:
            print(self.params)

    def run(self):
        qry = self.params['query']
        ref = self.params['reference']
        df_q = read_to_df(qry)
        df_r = read_to_df(ref)
        func = getattr(self, self.params['method'])
        df = func(df_q, df_r)
        sort = self.params['sort']
        df = aligner.filter_idx(df, 'query_id', sort[0], sort[1])
        df.to_csv(self.params['ofile'], index=False)
        return df

    def search_DNA(qry, ref, fwd_only=False, exact=False, circular=False):
        '''
        Find all instances of the query sequence in the reference.
        qry = DNA subsequence to find
        ref = DNA sequence containing the query subsequence
        fwd_only = search in forward orientation 5'->3' only. Defaults to searching forward and reverse complement
        exact = if False, then allow matching to ambiguous bases
        circular = search circularly
        returns a dataframe with columns [start, end, strand]
        '''
        # make it case in sensitive
        qry = str(qry).upper()
        ref = str(ref).upper()
        qry = qry.replace('-','')
        qry = qry.replace(' ','')
        L = len(ref)

        if fwd_only:
            if circular:
                ref+=ref
            out = []
            if exact:
                # look for stuff in forward 5'->3'
                x = ref.find(qry)
                idx = 0
                # look for every occurrence
                while x > -1:
                    s1 = x+idx
                    s2 = s1 + len(qry)
                    out.append([s1, s2, 1])
                    # shift the search window
                    idx = idx + x + 1
                    x = ref[idx:].find(qry)
            else:
                # search with ambiguous bases
                fwd = Bio.SeqUtils.nt_search(ref, qry)
                for i in fwd[1:]:
                    out.append([i, i+len(qry), 1])
            if circular and len(out) > 0:
                out = np.array(out)
                out = out[out[:,0] < L]
                out[:,1] = out[:,1]%L
                out[out[:,1]==0,1] = L
        else:
            # fwd and reverse search
            fwd_out = aligner.search_DNA(qry, ref, fwd_only=True, exact=exact, circular=circular)
            rev_ref = base.revcomp(ref)
            rev_out = aligner.search_DNA(qry, rev_ref, fwd_only=True, exact=exact, circular=circular)
            # correct reverse locations
            for i in range(len(rev_out)):
                s1 = rev_out.iat[i,0]
                s2 = rev_out.iat[i,1]
                rev_out.iat[i,0] = len(ref) - s2 
                rev_out.iat[i,1] = len(ref) - s1
                rev_out.iat[i,2] = -1
            return pd.concat([fwd_out, rev_out])

        # write the output
        return pd.DataFrame(out, columns=['start','end','strand'])

    def search_protein(qry, ref, fwd_only=False, circular=False):
        '''
        Finds all instances of the query protein sequence in the reference DNA sequence. Forward and reverse complements are searched.
        qry = protein sequence
        ref = DNA sequence
        circular = assume reference sequence is circular and search around origin
        returns a dataframe with columns [start, end, strand]
        '''
        if type(qry) == Bio.SeqRecord.SeqRecord:
            qry = str(qry.seq)
        else:
            qry = str(qry)
        if type(ref) == Bio.SeqRecord.SeqRecord:
            ref = str(ref.seq)
        else:
            ref = str(ref)
        
        # look for stuff in forward 5'->3' at each frame shift
        if fwd_only:
            out = []
            for i in range(0,3):
                if circular:
                    db = Bio.Seq.Seq(ref+ref)[i:].translate()
                else:
                    db = Bio.Seq.Seq(ref)[i:].translate()
                df = aligner.search_DNA(qry, db, fwd_only=True, exact=True, circular=False)
                for [start, end, strand] in df[['start','end','strand']].values:
                    # convert back to base pair location
                    start = start*3 + i
                    end = start + len(qry)*3
                    out.append([start, end, 1])
            # correct for wraps around the origin
            if circular and len(out) > 0:
                L = len(ref)
                out = np.array(out)
                out = out[out[:,0] < L]
                out[:,1] = out[:,1]%L
                out[out[:,1]==0,1] = L
            # write the output
            return pd.DataFrame(out, columns=['start','end','strand'])
        else:
            # fwd and reverse search
            fwd_out = aligner.search_protein(qry, ref, fwd_only=True, circular=circular)
            rev_ref = base.revcomp(ref)
            rev_out = aligner.search_protein(qry, rev_ref, fwd_only=True, circular=circular)
            # correct reverse locations
            for i in range(len(rev_out)):
                s1 = rev_out.iat[i,0]
                s2 = rev_out.iat[i,1]
                rev_out.iat[i,0] = len(ref) - s2 
                rev_out.iat[i,1] = len(ref) - s1
                rev_out.iat[i,2] = -1
            return pd.concat([fwd_out, rev_out])

    def search_ORF(self, sequence):
        '''
        Search for open reading frames, todo
        '''
        return -1

    def blast(self, sequence):
        '''
        Search sequence on NIH blast
        '''
        return -1

    def spoa(seq, verbose=True, width=100):
        '''
        Runs multi-sequence alignment on provided sequences with spoa
        seq = list of sequences
        verbose = print out results
        width = characters per line
        returns consensus, msa
        '''
        consensus, msa = spoa.poa(seq)
        if verbose:
            x = np.arange(0,len(msa[0]),width)
            if x[-1] < len(msa[0]):
                x = np.append(x,len(msa[0]))
            output = ''
            for i in range(len(x)-1):
                for j in range(len(msa)):
                    output+= msa[j][x[i]:x[i+1]]+' '+str(x[i+1])+' seq'+str(j)+'\n'
                # find mismatches
                match = np.array([True]*(x[i+1]-x[i]))
                for j in range(len(msa)-1):
                    seq1 = np.array([k for k in msa[j][x[i]:x[i+1]]])
                    seq2 = np.array([k for k in msa[j+1][x[i]:x[i+1]]])
                    match*= (seq1==seq2)
                mchar = np.array([' ']*(x[i+1]-x[i]))
                mchar[~match] = '.'
                output+= ''.join(mchar)+'\n'
            print(output) 
        return consensus, msa

    def string_match(self, qry, ref):
        '''
        Exact string search against query and reference
        '''
        qry = aligner.str_to_df(qry, 'query')
        ref = aligner.str_to_df(ref, 'ref')
        out = []
        for rname, rseq in ref[['name','sequence']].values:
            for qname, qseq in qry[['name','sequence']].values:
                x = aligner.search_DNA(qseq, rseq)
                x['query_id'] = qname
                x['ref_id'] = rname
                out.append(x)
        return pd.concat(out)

    def minimap2(self, qry, ref):
        '''
        Run minimap aligner
        qry = dataframe of query sequences with columns [name, sequence]
        ref = dataframe of reference sequences with columns [name, sequence]
        k = kmer length
        w = window size   
        '''
        k = self.params['minimap']['k']
        w = self.params['minimap']['w']
        qry = aligner.str_to_df(qry, 'query')
        ref = aligner.str_to_df(ref, 'ref')

        out = []
        for rname, rseq in ref[['name','sequence']].values:
            a = mappy.Aligner(seq=rseq)
            for qname, qseq in qry[['name','sequence']].values:
                for h in a.map(qseq):
                    out.append([qname, len(qseq), h.q_st, h.q_en, rname, h.ctg_len, h.r_st, h.r_en, h.trans_strand, h.read_num,
                                h.blen, h.mlen, h.mlen/h.blen, h.NM, h.mapq, h.is_primary, h.cigar_str])
        col = ['query_id','q_len','q_start','q_end','ref_id','t_len','t_start','t_end','strand','read_num',
               'blen','mlen','match_score','NM','mapq','is_primary','cigar']
        return pd.DataFrame(out, columns=col)

    def parasail(self, qry, ref):
        '''
        Run parasail aligner
        qry = dataframe of query sequences with columns [name, sequence]
        ref = dataframe of reference sequences with columns [name, sequence]
        '''
        g1 = self.params['parasail']['g1']
        g2 = self.params['parasail']['g2']
        mat = self.params['parasail']['mat']
        method = self.params['parasail']['method']

        qry = aligner.str_to_df(qry, 'query')
        ref = aligner.str_to_df(ref, 'ref')

        out = []
        # select aligner function to use
        paligner = getattr(parasail, method)
        for rname, rseq in ref[['name','sequence']].values:
            for qname, qseq in qry[['name','sequence']].values:
                qseqR = base.revcomp(qseq)
                res = paligner(qseq, rseq, g1, g2, mat)
                resR = paligner(qseqR, rseq, g1, g2, mat)
                cigarF = ''
                cigarR = ''
                if 'trace' in method:
                    cigarF = res.cigar.decode.decode('UTF-8')
                    cigarR = resR.cigar.decode.decode('UTF-8')
                sF = aligner.cigar_to_stats(cigarF)
                sR = aligner.cigar_to_stats(cigarR)
                out.append([qname, len(qseq), rname, len(rseq), 1, res.score, sF[1]+sF[2], sF[0], sF[1], sF[2], cigarF])
                out.append([qname, len(qseq), rname, len(rseq), -1, resR.score, sR[1]+sR[2], sR[0], sR[1], sR[2], cigarR])
        col = ['query_id','q_len','ref_id','t_len','strand','AS','NM','match','mismatch','indel','cigar']
        return pd.DataFrame(out, columns=col)

    def filter_idx(df, col, value='score', opt='idxmax'):
        dfi = df.reset_index(drop=True)
        idx = df.groupby(by=col).agg({value:opt}).values[:,0]
        return df.iloc[idx]

    def str_to_df(x, header='seq'):
        # format query and ref to dataframe if not already
        if type(x) == str:
            x = [x]
        if type(x) == list or type(x) == np.ndarray:
            x = pd.DataFrame(x, columns=['sequence'])
            x['name'] = [header+'_'+str(i) for i in range(len(x))]
        return x

    def cigar_to_list(cigar):
        pattern = re.compile('([0-9]*)([DIM=SX])')
        return [[int(n), c] for n,c in pattern.findall(cigar)]

    def cigar_to_stats(cigar):
        '''
        Return total [match, mismatch, indel] in cigar string
        '''
        m = 0
        x = 0
        indel = 0
        for n,c in aligner.cigar_to_list(cigar):
            if c in ['M','=']:
                m+=n
            if c in ['X']:
                x+=n
            if c in ['D','I']:
                indel+=n
        return [m, x, indel]

    def cigar_to_alignment(query, ref, cigar):
        '''
        Converts cigar to alignment text
        query = query sequence string
        ref = reference sequence string
        cigar = cigar string
        '''
        qa = ''
        ra = ''
        align = ''
        i1 = 0
        i2 = 0
        for n, c in aligner.cigar_to_list(cigar):
            span = int(n)
            if c == 'M' or c=='=':
                qa+=query[i1:i1+span]
                ra+=ref[i2:i2+span]
                align+=''.join(['|']*span)
                i1+=span
                i2+=span
            elif c=='X':
                qa+=query[i1:i1+span]
                ra+=ref[i2:i2+span]
                align+=''.join([' ']*span)
                i1+=span
                i2+=span
            elif c=='S' or c=='H':
                qa+=query[i1:i1+span]
                ra+=ref[i2:i2+span]
                align+=''.join([' ']*span)
                i1+=span
                i2+=span
            elif c=='I':
                qa+=query[i1:i1+span] # contains insertion
                ra+=''.join(['-']*span) # adds gap
                align+=''.join([' ']*span)
                i1+=span
            elif c=='D':
                qa+=''.join(['-']*span) # contains deletion --> add a gap
                ra+=ref[i2:i2+span]
                align+=''.join([' ']*span)
                i2+=span
        return [qa, ra, align]
      
