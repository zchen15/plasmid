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
from .misc import revcomp

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
    params['bin'] = {'minimap':'minimap2',
                     'bwa':'bwa',
                     'bowtie2':'bowtie2'}

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
        df = Aligner.filter_idx(df, 'query_id', sort[0], sort[1])
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
            fwd_out = Aligner.search_DNA(qry, ref, fwd_only=True, exact=exact, circular=circular)
            rev_ref = revcomp(ref)
            rev_out = Aligner.search_DNA(qry, rev_ref, fwd_only=True, exact=exact, circular=circular)
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
                df = Aligner.search_DNA(qry, db, fwd_only=True, exact=True, circular=False)
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
            fwd_out = Aligner.search_protein(qry, ref, fwd_only=True, circular=circular)
            rev_ref = revcomp(ref)
            rev_out = Aligner.search_protein(qry, rev_ref, fwd_only=True, circular=circular)
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
        qry = Aligner.str_to_df(qry, 'query')
        ref = Aligner.str_to_df(ref, 'ref')
        out = []
        for rname, rseq in ref[['name','sequence']].values:
            for qname, qseq in qry[['name','sequence']].values:
                x = Aligner.search_DNA(qseq, rseq)
                x['query_id'] = qname
                x['ref_id'] = rname
                out.append(x)
        return pd.concat(out)
    
    def fastqc(self, qry):
        '''
        Wrapper for fastqc
        '''
        return -1

    def prodigal(self, qry):
        '''
        Wrapper for prodigal
        '''
        return -1
    
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
        qry = Aligner.str_to_df(qry, 'query')
        ref = Aligner.str_to_df(ref, 'ref')

        out = []
        for rname, rseq in ref[['name','sequence']].values:
            a = mappy.aligner(seq=rseq)
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

        qry = Aligner.str_to_df(qry, 'query')
        ref = Aligner.str_to_df(ref, 'ref')

        out = []
        # select aligner function to use
        paligner = getattr(parasail, method)
        for rname, rseq in ref[['name','sequence']].values:
            for qname, qseq in qry[['name','sequence']].values:
                qseqR = revcomp(qseq)
                res = paligner(qseq, rseq, g1, g2, mat)
                resR = paligner(qseqR, rseq, g1, g2, mat)
                cigarF = ''
                cigarR = ''
                if 'trace' in method:
                    cigarF = res.cigar.decode.decode('UTF-8')
                    cigarR = resR.cigar.decode.decode('UTF-8')
                sF = Aligner.cigar_to_stats(cigarF)
                sR = Aligner.cigar_to_stats(cigarR)
                out.append([qname, len(qseq), rname, len(rseq), 1, res.score, sF[1]+sF[2], sF[0], sF[1], sF[2], cigarF])
                out.append([qname, len(qseq), rname, len(rseq), -1, resR.score, sR[1]+sR[2], sR[0], sR[1], sR[2], cigarR])
        col = ['query_id','q_len','ref_id','t_len','strand','AS','NM','match','mismatch','indel','cigar']
        return pd.DataFrame(out, columns=col)

    def filter_idx(df, col, value='score', opt='idxmax'):
        dfi = df.reset_index(drop=True)
        idx = df.groupby(by=col).agg({value:opt}).values[:,0]
        return df.iloc[idx]

    def str_to_df(x, header='seq'):
        '''
        Format query and ref to dataframe if not already        
        x = string or list of strings
        returns dataframe
        '''
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
        for n,c in Aligner.cigar_to_list(cigar):
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
        for n, c in Aligner.cigar_to_list(cigar):
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

    def run_minimap2(query, database, workspace='./minimap2/', config='-x map-ont', cigar=True, build_index=False, use_index=False, cleanup=False):
        '''
        This is a wrapper for the minimap2 aligner. This aligner is better suited for long reads (>100bp)

        query = list of sequences to search for in database
        query must be in pandas dataframe format with columns:
        [id, sequence] or [id, sequence, quality]
        [fwd_id, fwd_seq, rev_id, rev_seq] or
        [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]

        database = database of sequences being search through
        database must be in pandas dataframe format with columns:
        [id, sequence] or [id, sequence, quality]

        build_index = determines if a burrow-wheeler transform index should be build for the database
        use_index = determines if index will be used or just read from database fasta file
        workspace = defines where input and output files for bowtie2 resides
        config = defines the config options to pass to bowtie2
        '''
        workspace = check_dir(workspace)
        # load sequence incase low mem option is used 
        query = add_seq(query) 
        database = add_seq(database)
        # Define locations of input and output files
        qfile = workspace+'read1.fq'
        mfile = workspace+'read2.fq'
        dbfile = workspace+'database.fq'
        outfile = workspace+'results.paf'
        btfile = workspace+'index.mmi'
        # remove spaces from id and save the mapping
        if 'fwd_id' in query.columns:
            query, fmap = aligner_fix_name(query, col='fwd_id')
            query, rmap = aligner_fix_name(query, col='rev_id')
        else:
            query, qmap = aligner_fix_name(query, col='id')
        database, dbmap = aligner_fix_name(database, col='id')

        # write aligner input 
        paired = aligner_write_input(qfile, mfile, dbfile, query, database)
        if cigar: # adds flag to generate cigar string
            config = '-c '+config
        # Only build index when we need to. Index building takes a lot of time
        if build_index:
            cmd = _minimap2+' '+config+' -d '+btfile+' '+dbfile
            subprocess.run(cmd.split(' '))
        # Switch to index if needed
        if use_index:
            dbfile = btfile
        # Do paired end search if needed
        if paired:
            cmd = _minimap2+' '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
        else:
            cmd = _minimap2+' '+config+' '+dbfile+' '+qfile+' -o '+outfile
        # call the aligner
        subprocess.run(cmd.split(' '))
        # extract the PAF file information
        data = parse_PAF(outfile)
        # add spaces back to the sequences
        if 'fwd_id' in query.columns:
            data = aligner_restore_name(data, fmap, col='query_id')
            data = aligner_restore_name(data, rmap, col='query_id')
        else:
            data = aligner_restore_name(data, qmap, col='query_id')
        data = aligner_restore_name(data, dbmap, col='database_id')

        # Format the fields to integer
        x = ['q_len','q_start','q_end','t_len','t_start','t_end','match','tot','mapq','cm','s1','s2','NM','AS','ms','nn','rl']
        for i in x:
            data[i] = data[i].astype(int)
        data['match_score'] = data['match']/data['tot']*1.0 # compute blast like score
        # compute similarity
        data = aligner_get_similarity(data)
        # drop cigar if it is not present
        if cigar == False:
            data.drop(columns = ['CIGAR'], inplace = True)
        # remove work folder
        if cleanup and ~use_index:
            subprocess.run(['rm','-r',workspace])
        return data

    def run_bowtie2(query, database, workspace='./bowtie2/', config='-a --very-sensitive-local --threads 1 --quiet', build_index=True, cleanup=False):
        '''
        This is a wrapper for the bowtie2 aligner.

        query = list of sequences to search for in database
        query must be in pandas dataframe format with columns:
        [id, sequence] or [id, sequence, quality]

        database = database of sequences being search through
        database must be in pandas dataframe format with columns:
        [id, sequence] or [id, sequence, quality]
        [fwd_id, fwd_seq, rev_id, rev_seq] or
        [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]

        build_index = determines if a burrow-wheeler transform index should be build for the database
        workspace = defines where input and output files for bowtie2 resides
        config = defines the config options to pass to bowtie2
        '''
        workspace = check_dir(workspace)
        # load sequence incase low mem option is used 
        query = add_seq(query) 
        database = add_seq(database)
        # Define locations of input and output files for bowtie
        qfile = workspace+'read1.fa'
        mfile = workspace+'read2.fa'
        dbfile = workspace+'database.fa'
        outfile = workspace+'results.sam'
        btfile = workspace+'index'
        # write aligner input
        paired = aligner_write_input(qfile, mfile, dbfile, query, database)
        # Only build index when we need to. Index building takes a lot of time
        if build_index == True:
            value = check_seqdf(database)
            if value == 2:
                cmd = _bowtie2+'-build --quiet -q '+dbfile+' '+btfile+' --threads 2'
            elif value == 1:
                cmd = _bowtie2+'-build --quiet -f '+dbfile+' '+btfile+' --threads 2'
            subprocess.run(cmd.split(' '))
        # call bowtie2
        subprocess.run(cmd.split(' '))
        data = parse_SAM(outfile)
        # add q_len and t_len info
        q_len = np.transpose([query['id'].values, [len(i) for i in query['sequence']]])
        q_len = pd.DataFrame(q_len, columns = ['query_id','q_len'])
        data = data.merge(q_len, on='query_id', how='left')
        t_len = np.transpose([database['id'].values, [len(i) for i in database['sequence']]])
        t_len = pd.DataFrame(t_len, columns = ['database_id','t_len'])
        # debug
        data = data.merge(t_len, on='database_id', how='left')
        data = data.dropna() # drop shit that doesn't have sequence length
        # Format the fields to integer
        x = ['q_len','t_len','t_start','t_end','mapq','AS','XS','XN','XM','XO','XG','NM']
        for i in x:
            data[i] = data[i].astype(int)
        # compute similarity
        data = aligner_get_similarity(data)
        # remove work folder
        if cleanup:
            subprocess.run(['rm', '-r', workspace])
        return data

    def run_bwa(query, database, workspace='./bwa/', config=' mem -a ', build_index=True, cleanup=False):
        '''
        This is a wrapper for the bwa aligner.

        query = list of sequences to search for in database
        query must be in pandas dataframe format with columns:
        [id, sequence] or [id, sequence, quality]
        [fwd_id, fwd_seq, rev_id, rev_seq] or
        [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]

        database = database of sequences being search through
        database must be in pandas dataframe format with columns:
        [id, sequence] or [id, sequence, quality]

        build_index = determines if a burrow-wheeler transform index should be build for the database
        workspace = defines where input and output files for bowtie2 resides
        config = defines the config options to pass to bowtie2
        '''
        workspace = check_dir(workspace)
        # load sequence incase low mem option is used 
        query = add_seq(query) 
        database = add_seq(database)
        # Define locations of input and output files for bowtie
        qfile = workspace+'read1.fq'
        mfile = workspace+'read2.fq'
        dbfile = workspace+'database.fa'
        outfile = workspace+'results.sam'
        # write aligner input
        paired = aligner_write_input(qfile, mfile, dbfile, query, database)
        # Only build index when we need to. Index building takes a lot of time
        if build_index:
            cmd = _bwa+' index '+dbfile
            subprocess.run(cmd.split(' '))
        # do paired end search if needed
        if paired:
            cmd = _bwa+' '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
        else:
            cmd = _bwa+' '+config+' '+dbfile+' '+qfile+' -o '+outfile
        # call the aligner
        subprocess.run(cmd.split(' '))
        # extract the SAM file information
        data = parse_SAM(outfile)
        # add q_len and t_len info
        q_len = np.transpose([query['id'].values, [len(i) for i in query['sequence']]])
        q_len = pd.DataFrame(q_len, columns = ['query_id','q_len'])
        data = data.merge(q_len, on='query_id', how='left')
        t_len = np.transpose([database['id'].values, [len(i) for i in database['sequence']]])
        t_len = pd.DataFrame(t_len, columns = ['database_id','t_len'])
        data = data.merge(t_len, on='database_id', how='left')
        data = data.dropna() # drop shit that doesn't have sequence length
        # Format the fields to integer
        x = ['q_len','t_len','t_start','t_end','mapq','AS','XS','XN','XM','XO','XG','NM']
        for i in x:
            data[i] = data[i].astype(int)
        # compute similarity
        data = aligner_get_similarity(data)
        # remove work folder
        if cleanup:
            subprocess.run(['rm','-r',workspace])
        return data
