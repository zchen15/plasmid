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

# system and time
import re
import tempfile
import subprocess

# import custom libraries
from .misc import read_to_df
from .misc import revcomp
from .misc import load_args
from .fileIO import fileIO

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
                          'bowtie2':'bowtie2'}

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

    def search_DNA(qry, ref, fwd_only=False, circular=False):
        '''
        Find all instances of the query sequence in the reference.
        qry = DNA subsequence to find
        ref = DNA sequence containing the query subsequence
        fwd_only = search in forward orientation 5'->3' only. Defaults to searching forward and reverse complement
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
            fwd_out = Aligner.search_DNA(qry, ref, fwd_only=True, circular=circular)
            rev_ref = revcomp(ref)
            rev_out = Aligner.search_DNA(qry, rev_ref, fwd_only=True, circular=circular)
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

    def search_ORF(self, sequence):
        '''
        Search for open reading frames, to do
        '''
        return -1

    def blast(self, sequence):
        '''
        Search sequence on NIH blast, to do
        '''
        return -1

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
        # to do add query start and end on ref
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
            qfile = query
        else:
            query, qmap = Aligner.aligner_fix_name(query, col='name')
            qfile = work_dir.name + '/' + 'read1'
            qfile = fileIO.write_fastx(qfile, query)
        
        if type(query2).__name__!='NoneType':
            paired = True
            if type(query2) == str:
                mfile = query2
            else:
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
        data['similarity'] = 1-data['NM']/v1
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
            cmd = bwa+' '+config+' '+dbfile+' '+qfile+' -o '+outfile
        else:
            cmd = bwa+' '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
        subprocess.run(cmd.split(' '))
        
        # extract the SAM file information
        data = fileIO.parse_SAM(outfile)
        
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

    def bowtie2_get_index(self, database, filename='bt_index'):
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
            dbfile = database
        else:
            dbfile = work_dir.name + '/db'
            dbfile = fileIO.write_fastx(dbfile, database)

        btfile = work_dir.name+'/'+filename    
        if '.fq' in dbfile:
            cmd = bowtie2+' -build --quiet -q '+dbfile+' '+btfile+' --threads 2'
        elif '.fa' in dbfile:
            cmd = bowtie2+' -build --quiet -f '+dbfile+' '+btfile+' --threads 2'
        subprocess.run(cmd.split(' '))
        self.params['bowtie2']['index'] = btfile
        return btfile
    
    def bowtie2(query, query2=None, database=None):
        '''
        This is a wrapper for the bowtie2 aligner.

        query = list of sequences to search for in database
        query must be in pandas dataframe format with columns:
        [name, sequence] or [name, sequence, quality]
        
        query2 = reverse pair of paired end read
        query must be in pandas dataframe format with columns:
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
        
        if query2!=None:
            if type(query2) == str:
                mfile = query2
            else:
                mfile = work_dir.name + '/' + 'read2'
                mfile = fileIO.write_fastx(mfile, query2)
        
        if database!=None:
            self.bowtie2_get_index(database)
        btfile = self.params['bowtie2']['index']
        outfile = work_dir.name+'/results.sam'
        config = self.params['bowtie2']['config']
        bowtie2 = self.params['binaries']['bowtie2']
        
        if query2==None and '.fa' in qfile:
            cmd = bowtie2+' '+config+' -x '+bt_file+' -q '+qfile+' -S '+outfile
        elif query2==None and '.fq' in qfile:
            cmd = bowtie2+' '+config+' -x '+bt_file+' -f '+qfile+' -S '+outfile
        else:
            cmd = bowtie2+' '+config+' -x '+bt_file+' -1 '+qfile+' -2 '+mfile+' -S '+outfile
        # call bowtie2
        subprocess.run(cmd.split(' '))
        data = fileIO.parse_SAM(outfile)
        
        # add q_len and t_len info
        q_len = np.transpose([query['name'].values, [len(i) for i in query['sequence']]])
        q_len = pd.DataFrame(q_len, columns = ['query_id','q_len'])
        data = data.merge(q_len, on='query_id', how='left')
        t_len = np.transpose([database['name'].values, [len(i) for i in database['sequence']]])
        t_len = pd.DataFrame(t_len, columns = ['database_id','t_len'])
        
        # debug to do
        data = data.merge(t_len, on='database_id', how='left')
        data = data.dropna() # drop shit that doesn't have sequence length

        # Format the fields to integer
        x = ['q_len','t_len','t_start','t_end','mapq','AS','XS','XN','XM','XO','XG','NM']
        for i in x:
            data[i] = data[i].astype(int)
        
        # compute similarity
        data = Aligner.get_similarity(data)

        return data

