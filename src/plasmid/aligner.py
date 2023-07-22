#!/usr/bin/env python
# Class function for aligners

# computing libraries
import numpy as np
import pandas as pd

# reading and writing genbank files
import Bio
import Bio.SeqUtils
import spoa
import parasail
import mappy
import sklearn

# system and time
import re
import tempfile
import subprocess

# import custom libraries
from .misc import read_to_df
from .misc import revcomp
from .misc import translate
from .misc import substring_search
from .misc import isDNA
from .misc import load_args
from .fileIO import fileIO
from .graphic import Graphic
from .graphic import colormaps

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
                         'w':1,
                         'config':'-x map-ont -P',
                         'cigar':True}

    params['bwa'] = {'config':'mem -a'}
    params['bowtie2'] = {'config':'-a --very-sensitive-local --threads 1 --quiet'}
    
    params['spoa'] = {'algorithm':0,
                      'genmsa':True,
                      'm':5,
                      'n':-4,
                      'g':-8,
                      'e':-6,
                      'q':-10,
                      'c':-4}
    
    params['work_dir'] = None
    params['binaries'] = {'minimap':'minimap2',
                          'bwa':'bwa',
                          'bowtie2':'bowtie2',
                          'mmseqs2':'mmseqs',
                          'ngmerge':'./ngmerge'}

    def __init__(self, args=None):
        if args!=None:
            self.params = load_args(args, self.params)
        if self.params['verbose']:
            print(self.params)
        self.create_tempdir()

    def __del__(self):
        '''
        Destructor class to clean up temporary files
        '''
        folder = self.params['work_dir']
        folder.cleanup()
    
    def __enter__(self):
        '''
        Entry function for use in with statements
        '''
        return self
    
    def __exit__(self):
        '''
        Exit function for use in with statements
        '''
        del(self)
    
    def create_tempdir(self):
        '''
        Create a temporary working directory for aligners
        '''
        self.params['work_dir'] = tempfile.TemporaryDirectory(prefix='aligner_')
        
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

    def string_match(self, qry:list, ref:list) -> pandas.core.frame.DataFrame:
        '''
        Performs exact string search of query sequences against reference sequences
        Args:
            qry (list) : list of query sequences 
            ref (list) : list of reference sequences

        Returns:
            pandas.core.frame.DataFrame: A pandas DataFrame with columns [query_id, ref_id, start, end, strand]
        '''
        out = []
        for i in range(len(ref)):
            for j in range(len(qry)):
                x = Aligner.search_DNA(qseq[j], rseq[i])
                x['query_id'] = 'query_'+str(i)
                x['ref_id'] = 'ref_'+str(j)
                out.append(x)
        return pd.concat(out)
 
    def search_DNA(qry:str, ref:str, fwd_only:bool=False, circular:bool=False, exact:bool=False) -> pandas.core.frame.DataFrame:
        '''
        Find all instances of the query sequence in the reference
        Args:
            qry (str): DNA subsequence to find
            ref (str): DNA sequence containing the query subsequence
            fwd_only (bool): search in forward orientation 5'->3' only. Defaults to searching forward and reverse complement
            circular (bool): search string as though its a circular plasmid sequence
            exact (bool): use exact substring search

        Returns:
            pandas.core.frame.DataFrame: pandas dataframe with columns [start, end, strand]
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
            if isDNA(qry) and exact==False:
                # search with ambiguous bases
                fwd = Bio.SeqUtils.nt_search(ref, qry)
                for i in fwd[1:]:
                    out.append([i, i+len(qry), 1])
            else:
                # use exact string search
                fwd = substring_search(ref, qry)
                for i in fwd:
                    out.append([i, i+len(qry), 1])

            if circular and len(out) > 0:
                out = np.array(out)
                out = out[out[:,0] < L]
                out[:,1] = out[:,1]%L
                out[out[:,1]==0,1] = L
        else:
            # fwd and reverse search
            fwd_out = Aligner.search_DNA(qry, ref, fwd_only=True, circular=circular, exact=exact)
            rev_ref = revcomp(ref)
            rev_out = Aligner.search_DNA(qry, rev_ref, fwd_only=True, circular=circular, exact=exact)
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
                df = Aligner.search_DNA(qry, db, fwd_only=True, circular=False)
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

    def search_ORF(seq:str, table:str='Standard', fwd_only:bool=False) -> pandas.core.frame.DataFrame:
        '''
        Generates a dataframe of open reading frames for a given sequence
        
        Args:
            seq (str): input sequence
            table (str): codon table to use
            fwd_only (bool): search only in forward orientation

        Returns:
            pandas.core.frame.DataFrame: pandas dataframe with columns [name, start, stop, orientation, amino acid sequence]
        '''
        # search in forward direction
        out = []
        for i in range(3):
            # translate and get ORF
            x = translate(seq, frame=i)
            # convert to nt positions
            df = Aligner.get_ORF(x)
            df['start'] = df['start']*3+i 
            df['stop'] = (df['stop'])*3+i
            out.append(df)
        df = pd.concat(out)
        df['orientation'] = 1
        if fwd_only:
            return df
        
        # search in reverse direction
        rev_seq = revcomp(seq)
        df2 = Aligner.search_ORF(rev_seq, table=table, fwd_only=True)
        df2['orientation'] = -1
        stop = df2.stop.values
        df2['stop'] = len(seq) - df2['start']
        df2['start'] = len(seq) - stop
        return pd.concat([df,df2])
    
    def get_ORF(seq, start='M', stop='*'):
        '''
        Get list of ORF based on start and stop codons
        return list of sequences
        '''
        idx = np.arange(len(seq))
        seq_arr = np.array([i for i in seq])
        # get all positions of start and stop codons
        start = idx[seq_arr==start]
        stop = idx[seq_arr==stop]
        
        # generate ORF frags from start and stop positions
        out = []
        for i in range(len(start)):
            s1 = start[i]
            s2 = stop[stop > s1]
            # get first stop codon
            if len(s2) > 0:
                s2 = s2[0]
                out.append([s1, s2+1, seq[s1:s2+1]])
            # exit loop if no more stop codons are present
            else:
                break
        col = ['start','stop','sequence']
        out = pd.DataFrame(out, columns=col)
        return out
    
    def blast(self, sequence):
        '''
        Search sequence on NIH blast, to do
        '''
        return -1

   
    def spoa(self, seq):
        '''
        Runs multi-sequence alignment on provided sequences with spoa
        seq = list of sequences
        returns dictionary containing {'consensus':sequence, 'msa':<list of sequences>}
        '''
        algo = self.params['spoa']['algorithm']
        genmsa = self.params['spoa']['genmsa']
        m = self.params['spoa']['m']
        n = self.params['spoa']['n']
        g = self.params['spoa']['g']
        e = self.params['spoa']['e']
        q = self.params['spoa']['q']
        c = self.params['spoa']['c']
        
        consensus, msa = spoa.poa(sequences=seq, genmsa=genmsa, algorithm=algo, m=m, n=n, g=g, e=e, q=q, c=c)
        out = {'consensus':consensus,
               'msa':msa,
               'params':{}}
        for k in self.params['spoa'].keys():
            out['params'][k] = self.params['spoa'][k]
        return out
    
    def format_msa(msa, width=100):
        '''
        Returns a multi-sequence alignment formatted 
        to a certain width
        msa = list of sequences
        width = characters per line
        returns a string output
        '''
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
        return output
    
    def colorize_msa(msa, width=100, cmap=None):
        '''
        Returns a multi-sequence alignment formatted 
        to a certain width
        msa = list of sequences
        width = characters per line
        cmap = colormap from characters to colors
        returns a string output
        '''
        # use default nucleic acid colors
        if cmap==None:
            cmap = colormaps['nucleic']
        
        # initialize coloring engine
        gr = Graphic()
        
        # get character frequencies
        stats = Aligner.get_msa_stats(msa)
        
        # colorize characters in the minor population
        seq = np.array([[i for i in m] for m in msa])
        dim = seq.shape
        carr = np.array([['']*dim[1]]*dim[0])
        
        for i in range(len(msa[0])):
            chars = stats[i][0]
            freq = stats[i][1]
            b = 1
            # catch cases where at least two characters are frequent
            if len(chars) > 1 and freq[0] == freq[1]:
                b = 0
            for j in range(b,len(chars)):
                c = seq[:,i]==chars[j]
                carr[c,i] = chars[j]

        # format lines
        x = np.arange(0,len(msa[0]),width)
        if x[-1] < len(msa[0]):
            x = np.append(x,len(msa[0]))
            
        # format to colors
        output = ''
        for i in range(len(x)-1):
            for j in range(len(msa)):
                start = x[i]
                stop = x[i+1]
                text = msa[j][start:stop]
                fgarr = carr[j, start:stop]
                text = gr.color_by_array(text, cmap, bgarr=fgarr)
                output+= text+' '+str(x[i+1])+' seq'+str(j)+'\n'
        return output
    
    def get_msa_stats(msa):
        '''
        Compute counts for each character in the multi-alignment
        msa = list of sequences
        return list of [characters, frequency]
        '''
        seq = np.array([[i for i in m] for m in msa])
        seq = seq.T
        hist = [np.unique(s, return_counts=True) for s in seq]
        # sort the characters by frequency
        out = []
        for i in range(len(hist)):
            freq = hist[i][1]
            chars = hist[i][0]
            s = np.argsort(freq)
            s = s[::-1]
            out.append([chars[s], freq[s]])
        return out
    
    def mappy(self, qry, ref):
        '''
        Run minimap aligner
        qry = dataframe of query sequences with columns [name, sequence]
        ref = dataframe of reference sequences with columns [name, sequence]
        k = kmer length
        w = window size
        '''
        k = self.params['minimap']['k']
        w = self.params['minimap']['w']
        
        out = []
        for rname, rseq in ref[['name','sequence']].values:
            a = mappy.Aligner(seq=rseq, k=k, w=w)
            for qname, qseq in qry[['name','sequence']].values:
                for h in a.map(qseq):
                    out.append([qname, len(qseq), h.q_st, h.q_en, rname, h.ctg_len, h.r_st, h.r_en, 
                                h.trans_strand, h.read_num, h.blen, h.mlen, h.mlen/h.blen,
                                h.NM, h.mapq, h.is_primary, h.cigar_str])
        col = ['query_id','q_len','q_start','q_end','ref_id','t_len','t_start','t_end','strand','read_num',
               'blen','mlen','match_score','NM','mapq','is_primary','CIGAR']
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
                sF = Aligner.cigar_to_stats(cigarF, qseq, rseq)
                sF = [qname, rname, 1, res.score]+sF+[cigarF]
                sR = Aligner.cigar_to_stats(cigarR, qseq, rseq)
                sR = [qname, rname, -1, res.score]+sR+[cigarR]
                out.append(sF)
                out.append(sR)
        col = ['query_id','ref_id','strand','AS',
               'q_len','q_start','q_end',
               't_len','t_start','t_end',
               'match','mismatch','indel','CIGAR']
        df = pd.DataFrame(out, columns=col)
        df['NM'] = df['mismatch']+df['indel']
        return df
        
    def cigar_to_stats(cig, qseq, rseq):
        '''
        Get alignment statistics associated with parasail output
        '''
        cig = Aligner.cigar_decode(cig)
        
        # compute match boundaries
        q_start = 0
        q_end = len(qseq)
        q_len = len(qseq)
        t_start = 0
        t_end = len(rseq)
        t_len = len(rseq)
        if cig[0][1]=='I':
            t_start = cig[0][0]
        elif cig[0][1]=='D':
            q_start = cig[0][0]
        if cig[-1][1]=='I':
            t_end-= cig[-1][0]
        elif cig[-1][1]=='D':
            q_end-= cig[-1][0]
        
        x = np.array([cig[i][0] for i in range(len(cig))])
        y = np.array([cig[i][1] for i in range(len(cig))])
        
        # get number of matches
        M = np.sum(x[y=='='])
        M+= np.sum(x[y=='M'])
        
        # get number of mismatches
        Xm = np.sum(x[y=='X'])

        # get indel count
        indel = np.sum(x[y=='I'])
        indel+= np.sum(x[y=='D'])
        return [q_len, q_start, q_end, t_len, t_start, t_end, M, Xm, indel]
    
    def cigar_decode(cigar: str, key: str='([0-9]*)([DIM=SX])') -> list:
        '''
        Decode cigar strings into counts of indels and mismatches

        Args:
            cigar = cigar string
            key = key used to decrypt the cigar string
            
        Returns:
            list containing [[key value, counts],...]
        '''
        pattern = re.compile(key)
        return [[int(i),j] for i,j in pattern.findall(cigar)]

    def cigar_to_alignment(query: str, ref: str, cigar: str) -> list:
        '''
        Function to align query and reference strings based on CIGAR string

        Args:
            query = query sequence --> this must be relative to how it was used in aligner
            ref = reference sequence
            cigar = CIGAR string from aligner
        
        Returns:
            return [aligned query, aligned reference, alignment ticks]
        '''
        pattern = re.compile('([0-9]*)([DIMSX])')
        query = query.upper()
        ref = ref.upper()
        q_aligned = ''
        ref_aligned = ''
        align = ''
        i1 = 0
        i2 = 0
        for n, c in pattern.findall(cigar):
            span = int(n)
            if c == 'M':
                q_aligned+=query[i1:i1+span]
                ref_aligned+=ref[i2:i2+span]
                align+=''.join(['|']*span)
                i1+=span
                i2+=span
            elif c=='X':
                q_aligned+=query[i1:i1+span]
                ref_aligned+=ref[i2:i2+span]
                align+=''.join([' ']*span)
                i1+=span
                i2+=span
            elif c=='S' or c=='H':
                q_aligned+=query[i1:i1+span].lower() # query contains soft clipped sequence
                ref_aligned+=ref[i2:i2+span]
                align+=''.join([' ']*span)
                i1+=span
                i2+=span
            elif c=='I':
                q_aligned+=query[i1:i1+span] # query contains an insert
                ref_aligned+=''.join(['-']*span) # adds gap to reference
                align+=''.join([' ']*span)
                i1+=span
            elif c=='D':
                q_aligned+=''.join(['-']*span) # query contains a deletion --> add a gap
                ref_aligned+=ref[i2:i2+span]
                align+=''.join([' ']*span)
                i2+=span
        return [q_aligned, ref_aligned, align]
    
    def filter_idx(df: pd.DataFrame, col: list, value: str='score', opt: str='idxmax') -> pandas.core.frame.DataFrame:
        '''
        Get max, mean, median, or min values organized by a certain column in a pandas dataframe

        Args:
            col (str): column to group by
            value (str) : column to sort
            opt (str) : idxmax or idxmin

        Returns:
            pandas.core.frame.DataFrame: A pandas DataFrame of filtered values
        '''
        df=df.reset_index(drop=True)
        idx = df.groupby(by=col).agg({value:opt}).reset_index()
        return df.iloc[idx[value].values].reset_index()
    
   
    def minimap2_get_index(self, database, filename='index.mmi'):
        '''
        This is a wrapper for the minimap2 to build the fn_index
        fn_index file is built in the temporary file directory
        
        database = database of sequences being search through
        database must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        or a string defining the fasta or fastq file to index
        
        filename = defines the filename of the minimap index
        
        return filename of fn_index
        '''
        work_dir = self.params['work_dir']
        if type(database) == str:
            db_file = database
        else:
            db_file = work_dir.name + '/db'
            db_file = fileIO.write_fastx(db_file, database)

        # set the parameters
        w = self.params['minimap']['w']
        k = self.params['minimap']['k']
        mm2bin = self.params['binaries']['minimap']
        fn_index = work_dir.name + '/' + filename
        
        # execute the command
        cmd = mm2bin +' -k '+str(k)+' -w '+str(w)+' -d '+fn_index+' '+db_file
        subprocess.run(cmd.split(' '))
        self.params['minimap']['index'] = fn_index
        return fn_index
    
    def minimap2(self, query, query2=None, database=None):
        '''
        This is a wrapper for the minimap2 aligner.
        This aligner is better suited for long reads (>100bp) from nanopore sequencing

        query = list of sequences to search for in database
        query2 = reverse pair of paired end read

        query and query2 must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        
        database = database of sequences being search through

        database must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        or a string containing the filename of the fn_index

        returns a pandas dataframe of alignment results from the paf file
        '''
        # Define locations of input and output files
        work_dir = self.params['work_dir']
        if type(query) == str:
            query = fileIO.read_fastx(query)
        query, qmap = Aligner.aligner_fix_name(query, col='name')
        qfile = work_dir.name + '/' + 'read1'
        qfile = fileIO.write_fastx(qfile, query)
        
        paired = False
        if type(query2).__name__!='NoneType':
            paired = True
            if type(query2) == str:
                query2 = fileIO.read_fastx(query2)

            query2, qmap2 = Aligner.aligner_fix_name(query2, col='name')
            mfile = work_dir.name + '/' + 'read2'
            mfile = fileIO.write_fastx(mfile, query2)
        
        if type(database).__name__=='NoneType':
            dbfile = self.params['minimap']['index']
        elif type(database) == str:
            dbfile = database
        else:
            database, dbmap = Aligner.aligner_fix_name(database, col='name')
            dbfile = work_dir.name + '/' + 'database'
            dbfile = fileIO.write_fastx(dbfile, database)
        
        outfile = work_dir.name+'/results.paf'
        
        # set the parameters
        w = self.params['minimap']['w']
        k = self.params['minimap']['k']
        mm2bin = self.params['binaries']['minimap']
        
        config = self.params['minimap']['config']
        if self.params['minimap']['cigar']: # adds flag to generate cigar string
            config+= ' -c'
        config+= ' -k '+str(k)
        config+= ' -w '+str(w)

        # Do paired end search if needed
        if paired:
            cmd = mm2bin+' '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
        else:
            cmd = mm2bin+' '+config+' '+dbfile+' '+qfile+' -o '+outfile
        # call the aligner
        subprocess.run(cmd.split(' '))
        # extract the PAF file information
        data = fileIO.parse_PAF(outfile)
        # add spaces back to the sequences
        data = Aligner.aligner_restore_name(data, qmap, col='query_id')
        if query2!=None:
            data = Aligner.aligner_restore_name(data, qmap2, col='query_id')
        if type(database).__name__=='DataFrame':
            data = Aligner.aligner_restore_name(data, dbmap, col='database_id')

        # Format the fields to integer
        cols = ['q_len','q_start','q_end','t_len','t_start','t_end',
                'match','tot','mapq','cm','s1','s2','NM','AS','ms','nn','rl']
        for i in cols:
            data[i] = data[i].astype(int)
        data['match_score'] = data['match']/data['tot']*1.0 # compute blast like score

        # compute similarity
        data = Aligner.get_similarity(data)
        
        return data

    def aligner_fix_name(df, col):
        '''
        Removes spacing in the names that causes minimap2 to truncate read ids
        df = dataframe column containing the old id
        col = column with in the old id
        returns fixed dataframe and id mapping between old and new names
        '''
        out = df.copy()
        revmap = {str(i)+'_'+out.iloc[i][col].replace(' ','_'):out.iloc[i][col] for i in range(0,len(df))}
        out[col] = list(revmap.keys())
        return out, revmap

    def aligner_restore_name(df, revmap, col):
        '''
        Restore spaces in the names
        df = dataframe with column containing the new id
        revmap = dictionary containing mapping of new id to old id
        col = column name containing the id restore
        '''
        df[col] = [revmap[df.iloc[i][col]] for i in range(0,len(df))]
        return df

    def get_similarity(data):
        '''
        Add similarity metric to output dataframe from aligner
        data = dataframe with columns [NM, q_len, t_len]
        return dataframe with similarity added to columns
        '''
        v1 = np.min([data['q_len'].values, data['t_len'].values], axis = 0)
        data['similarity'] = 1 - data['NM']/v1
        data.loc[data['similarity'] < 0, 'similarity'] = 0
        return data

    def bwa_get_index(self, database):
        '''
        This is a wrapper for bwa to build the fn_index
        
        database = database of sequences being searched through
        
        database must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        or a string defining the fasta or fastq file to index
        
        return filename of bwa index
        '''
        bwa = self.params['binaries']['bwa']
        work_dir = self.params['work_dir']

        # write input files
        if type(database) == str:
            dbfile = database
        else:
            dbfile = work_dir.name + '/db'
            dbfile = fileIO.write_fastx(dbfile, database)
        cmd = bwa+' index '+dbfile
        subprocess.run(cmd.split(' '))
        self.params['bwa']['index'] = dbfile
        return dbfile

    def bwa(self, query, query2=None, database=None):
        '''
        This is a wrapper for the bwa aligner.

        query = list of sequences to search for in database
        query2 = reverse pair of paired end read

        query and query2 must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
                
        database = database of sequences being search through
        
        database must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        or a string containing the filename of the fn_index
        
        returns a pandas dataframe of alignment results from the sam file
        '''
        # Define locations of input and output files
        work_dir = self.params['work_dir']
        if type(query) == str:
            qfile = query
        else:
            qfile = work_dir.name + '/' + 'read1'
            qfile = fileIO.write_fastx(qfile, query)
        
        paired = False
        if type(query2).__name__!='NoneType':
            paired = True
            if type(query2) == str:
                mfile = query2
            else:
                mfile = work_dir.name + '/' + 'read2'
                mfile = fileIO.write_fastx(mfile, query2)
        
        if type(database).__name__!='NoneType':
            self.bwa_get_index(database)
        dbfile = self.params['bwa']['index']
        outfile = work_dir.name+'/results.sam'
        config = self.params['bwa']['config']
        bwa = self.params['binaries']['bwa']
        
        # do paired end search if needed
        if paired:
            cmd = bwa+' '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
        else:
            cmd = bwa+' '+config+' '+dbfile+' '+qfile+' -o '+outfile
        subprocess.run(cmd.split(' '))
        
        # extract the SAM file information
        data = fileIO.parse_SAM(outfile)
        data = Aligner.bwa_format_SAM(query, query2, database, data)
        return data

    def bwa_format_SAM(query, query2, database, data):
        '''
        Add q_len and t_len information to sam output
        query = str or pandas dataframe
        query2 = str or pandas dataframe
        database = str or pandas dataframe
        data = pandas dataframe from parse_SAM

        returns dataframe with q_len and t_len info
        '''
        paired = False
        if type(query2).__name__!='NoneType':
            paired = True
        
        if type(query)==str:
            query = fileIO.read_fastx(query)
        if type(query2)==str:
            query2 = fileIO.read_fastx(query2)
        if type(database)==str:
            database = fileIO.read_fastx(database)
        
        # add q_len and t_len info
        q_len = np.transpose([query['name'].values, [len(i) for i in query['sequence']]])
        q_len = pd.DataFrame(q_len, columns = ['query_id','q_len'])
        if paired:
            q_len2 = np.transpose([query2['name'].values, [len(i) for i in query2['sequence']]])
            q_len2 = pd.DataFrame(q_len, columns = ['query_id','q_len'])
            q_len = pd.concat([q_len, q_len2])
        data = data.merge(q_len, on='query_id', how='left')
        
        t_len = np.transpose([database['name'].values, [len(i) for i in database['sequence']]])
        t_len = pd.DataFrame(t_len, columns = ['database_id','t_len'])
        data = data.merge(t_len, on='database_id', how='left')
        data = data.dropna() # drop shit that doesn't have sequence length
        
        # Format the fields to integer
        x = ['q_len','t_len','t_start','t_end','mapq','AS','XS','XN','XM','XO','XG','NM']
        for i in x:
            data[i] = data[i].astype(int)
        
        # compute similarity
        data = Aligner.get_similarity(data)
        return data
    
    def bowtie2_get_index(self, database, filename='bowtie2_index'):
        '''
        This is a wrapper for bowtie2 to build the fn_index
        database = database of sequences being search through
        database must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        or a string defining the fasta or fastq file to index
        
        filename = defines the filename of the bwa index
        
        return filename of bwa index
        '''
        bowtie2 = self.params['binaries']['bowtie2']
        work_dir = self.params['work_dir']

        # write input files
        if type(database) == str:
            if '.fq' in database or '.fastq' in database:
                database = fileIO.read_fastq(database)
                dbfile = work_dir.name + '/db'
                fileIO.write_fasta(dbfile, database)
            else:
                dbfile = database
        else:
            dbfile = work_dir.name + '/db'
            fileIO.write_fasta(dbfile, database)

        btfile = work_dir.name+'/'+filename
        # bowtie2 only builds indexes against fasta files
        cmd = bowtie2+'-build --quiet -f '+dbfile+' '+btfile+' --threads 2'
        subprocess.run(cmd.split(' '))
        
        self.params['bowtie2']['index'] = btfile
        return btfile
    
    def bowtie2(self, query, query2=None, database=None):
        '''
        This is a wrapper for the bowtie2 aligner.

        query = list of sequences to search for in database
        query2 = reverse pair of paired end read
        
        query and query2 must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        
        database = database of sequences being search through
        database must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        or a string containing the filename of the fn_index
        
        returns a pandas dataframe of alignment results from the sam file
        '''
        # Define locations of input and output files
        work_dir = self.params['work_dir']
        if type(query) == str:
            qfile = query
        else:
            qfile = work_dir.name + '/' + 'read1'
            qfile = fileIO.write_fastx(qfile, query)
        
        paired = False
        if type(query2).__name__!='NoneType':
            paired = True
            if type(query2) == str:
                mfile = query2
            else:
                mfile = work_dir.name + '/' + 'read2'
                mfile = fileIO.write_fastx(mfile, query2)
        
        if type(database).__name__!='NoneType':
            self.bowtie2_get_index(database)
        btfile = self.params['bowtie2']['index']
        outfile = work_dir.name+'/results.sam'
        config = self.params['bowtie2']['config']
        bowtie2 = self.params['binaries']['bowtie2']
        
        if paired:
            cmd = bowtie2+' '+config+' -x '+btfile+' -1 '+qfile+' -2 '+mfile+' -S '+outfile            
        else:
            if '.fq' in qfile:
                cmd = bowtie2+' '+config+' -x '+btfile+' -q '+qfile+' -S '+outfile
            elif '.fa' in qfile:
                cmd = bowtie2+' '+config+' -x '+btfile+' -f '+qfile+' -S '+outfile
        # call bowtie2
        subprocess.run(cmd.split(' '))
        
        # extract the SAM file information
        data = fileIO.parse_SAM(outfile)
        data = Aligner.bwa_format_SAM(query, query2, database, data)
        return data

    def merge_paired_end(self, read1, read2, outfile=None):
        '''
        Merge paired end reads using NGmerge
        read1  = fastq reads fwd pair
        read2 = fastq reads rev pair
        
        return filename of output fastq
        '''
        # do paired alignment
        ngmerge = self.params['binaries']['ngmerge']
        work_dir = self.params['work_dir']
        if type(outfile)==type(None):
            outfile = work_dir.name+'/results.fq.gz'
        
        # run the read merge tool
        cmd = ngmerge+' -1 '+read1+' -2 '+read2+' -o '+outfile
        subprocess.run(cmd.split(' '))
        return outfile

    def get_fastq_statistics(query):
        '''
        Compute quality statistics about each read
        query and query2 must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        '''
        out = []
        for v in query['quality'].values:
            v = Aligner.phred_to_int(v)
            q = [0, 0.25, 0.5, 0.75, 1]
            y = np.quantile(v, q)
            out.append(y)

        col = ['Q_min','Q_lower','Q_median','Q_upper','Q_max']
        df = pd.DataFrame(out, columns=col)
        for c in col:
            query[c] = df[c]
        return query

    def phred_to_int(x):
        '''
        Translate phred ASCII to quality integers
        '''
        return [ord(i)-33 for i in x]

    def mmseq2_get_database(self, database, outdir):
        '''
        Wrapper for mmseq2 fetching a database
        database = database name options include
        
        UniRef100, Amino Acid db
        UniProtKB, Amino Acid db
        NR, Amino Acid Non-Redundant
        NT, Nucleotide, Nucleotide Non-Redundant
        GTDB, Genome Taxonomy Database
        Pfam-A.full, Profile
        SILVA, ribosomal RNA
        
        outdir = output directory
        '''
        mmseq = self.params['binaries']['mmseq2']
        work_dir = self.params['work_dir']
        # run clustering
        cmd = mmseq +' database '+database+' '+outdir+' '+work_dir.name
        subprocess.run(cmd.split(' '))

    def mmseq2_cluster(self, query, config='easy-cluster'):
        '''
        Wrapper for mmseq2 clustering
        '''
        mmseq = self.params['binaries']['mmseq2']
        work_dir = self.params['work_dir']
        
        if type(query) == str:
            qfile = query
        else:
            qfile = work_dir.name + '/' + 'query'
            qfile = fileIO.write_fastx(qfile, query)
        
        # run clustering
        cmd = mmseq +' '+config+' '+qfile+' '+ofile+' '+work_dir.name
        subprocess.run(cmd.split(' '))

    def mmseq2_search(self, query, config='easy-linsearch'):
        '''
        Wrapper for mmseq2 search
        '''
        mmseq = self.params['binaries']['mmseq2']
        
        cmd = mmseq +' '+config
        #mmseqs easy-cluster examples/DB.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1

    def mmseq2_taxon(self, query, config='easy-taxonomy'):
        '''
        Wrapper for mmseq2 search
        '''
        mmseq = self.params['binaries']['mmseq2']
        
        cmd = mmseq +' '+config
        #mmseqs easy-cluster examples/DB.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1

    def prodigal(self, query):
        '''
        Wrapper for prodigal
        '''
        return -1


    
    
