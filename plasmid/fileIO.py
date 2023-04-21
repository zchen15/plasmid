#!/usr/bin/env python
# Class function for handling sequence input output

# computing libraries
import numpy as np
import pandas as pd

# compression libraries
import bz2
import gzip
import tempfile

# system and time
import re
import os
import logging

# custom lib
from .misc import load_args

class fileIO:
    '''
    This class holds functions pertaining to sequence alignment and analysis
    '''
    params = {'verbose':False}
    
    def __init__(self, args=None):
        if args!=None:
            self.params = load_args(args, self.params)
        if self.params['verbose']:
            print(self.params)
    
    def add_seq(df):
        '''
        Adds sequence and quality information back to dataframe
        df = dataframe containing at least the columns:
            [filename, id, seek1, rlen] for fasta files
            [filename, id, seek1, seek2, rlen] for fastq files
        return dataframe with columns [id, sequence, quality] or [id, sequence]
        to do testing
        '''
        # if data is already there, dont do anything
        if 'sequence' in df.keys():
            return df
        # check if relevant columns are there
        cols = ['filename','name','seek1','rlen']
        for c in cols:
            if (c in df.keys())==False:
                logging.warning('columns '+c+' not found in dataframe')
                return df
        # sort and iterate
        df = df.sort_values(by=['filename','seek1']).copy()
        data = []
        if 'seek2' in df.keys():
            # work on fastq
            x = df[['filename','name','seek1','seek2','rlen']].values
            files = np.unique(x[:,0])
            for fname in files:
                y = x[x[:,0]==fname]
                # decompress file if needed
                with fileIO.decompress(fname) as f:
                    for i in range(0,len(y)):
                        f.seek(y[i,2])
                        seq = f.read(y[i,4])
                        f.seek(y[i,3])
                        q = f.read(y[i,4])
                        data.append([seq, q])
            df['sequence'] = [i[0] for i in data]
            df['quality'] = [i[1] for i in data]
            cols = ['seek1','seek2','rlen','filename']
            df = df.drop(columns=cols)
        else:
            # work on fasta
            x = df[['filename','name','seek1','rlen']].values
            files = np.unique(x[:,0])
            for fname in files:
                y = x[x[:,0]==fname]
                # decompress the file
                with fileIO.decompress(fname) as f:
                    for i in range(0,len(y)):
                        f.seek(y[i,2])
                        seq = f.read(y[i,3])
                        seq = seq.replace('\n','')
                        data.append(seq)
            df['sequence'] = data
            cols = ['seek1','rlen','filename']
            df = df.drop(columns=cols)
        return df

    def decompress(fname):
        '''
        Opens a compressed file
        fname = name of compressed file with .bz2 or .gz
        returns filehandle passed through the decompression 
        '''
        ftype = fname.split('.')
        if ftype[-1]=='bz2':
            return bz2.open(fname, mode='rt')
        elif ftype[-1]=='gz':
            return gzip.open(fname, mode='rt')
        else:
            return open(fname,'r')

    def compress(fname):
        '''
        Checks if a filename will be bz2 or gz compressed.
        fname = filename with .bz2 or .gz or no extension
        return a file handle
        '''
        ftype = fname.split('.')
        if ftype[-1]=='bz2':
            return bz2.open(fname, mode='wt')
        elif ftype[-1]=='gz':
            return gzip.open(fname, mode='wt')
        else:
            return open(fname,'w')

    def read_fastq(fname, low_mem=False):
        '''
        Reads list of sequences from a fastq file format.

        fname = fname to read the information
        low_mem = saves memory by not storing sequence and quality strings.
                Instead this function returns columns [filename, id, seek1, seek2, rlen]
                Use add_seq function to add back sequence and quality to the dataframe
        returns list of names, sequence, quality
        '''
        data = []
        with tempfile.TemporaryFile(mode='w+', encoding='utf-8') as f:
            # decompress and dump contents into a temporary file
            with fileIO.decompress(fname) as f2:
                out = f2.read()
                fsize = len(out)
                f.write(out)
            # go to beginning of file
            f.seek(0)
            while fsize > f.tell():
                rid = f.readline()[:-1]
                s1 = f.tell()
                seq = f.readline()[:-1]
                f.readline()
                s2 = f.tell()
                q = f.readline()[:-1]
                L = len(seq)
                if low_mem:
                    data.append([rid[1:], s1, s2, L])
                else:
                    data.append([seq, q, rid[1:]])
        # format the dataframe
        col = ['sequence','quality','name','seek1', 'seek2', 'rlen']
        if low_mem:
            data = pd.DataFrame(data, columns=col[2:])
            data['filename'] = fname
        else:
            data = pd.DataFrame(data, columns=col[:3])
        return data

    def write_fastq(fname, data):
        '''
        Write list of sequences into a fastq file format
        data = dataframe or numpy array with columns [id, sequence, quality]
        fname = file to write the information.
                File extension can be .bz2 or .gz if compression is desired
        '''
        # add some error handling in case dataframe is wrong
        if type(data) == np.ndarray and data.shape[1]!=3:
            logging.error('write_fastq: data for '+fname+' column count is wrong')
            sys.exit(1)
        elif type(data) == pd.core.frame.DataFrame:
            data = data[['name','sequence','quality']].values
        with fileIO.compress(fname) as f:
            for i in data:
                text = '@' + str(i[0]) + '\n'
                text+= str(i[1]) + '\n'
                text+= '+\n'
                text+= str(i[2]) + '\n'
                f.write(text)

    def read_fasta(fname, low_mem=False):
        '''
        Reads list of sequences from a fasta file format
        fname = file to read the information
        returns list of sequence and names
        '''
        data = []
        
        with tempfile.TemporaryFile(mode='w+', encoding='utf-8') as f:
            # decompress and dump contents into a temporary file
            with fileIO.decompress(fname) as f2:
                out = f2.read()
                fsize = len(out)
                f.write(out)
            # go to beginning of file
            f.seek(0)

            seq = ''
            rid = ''
            while fsize > f.tell():
                line = f.readline()
                if line[0]=='>':
                    # record it once new line is hit
                    if seq!='' and rid!='':
                        L = len(seq)
                        if low_mem:
                            data.append([rid, s1, L])
                        else:
                            seq = seq.replace('\n','')
                            data.append([seq, rid])
                    # start the new entry 
                    rid = line[1:-1]
                    s1 = f.tell()
                    seq = ''
                else:
                    seq+= line
        
        # end of file reached, run last loop
        L = len(seq)
        seq = seq.replace('\n','')

        # export the data
        col = ['sequence','name','seek1','rlen']
        if low_mem:
            data.append([rid, s1, L])
            data = pd.DataFrame(data, columns=col[1:])
            data['filename'] = fname
        else:
            data.append([seq, rid])
            data = pd.DataFrame(data, columns=col[:2])
        return data

    def write_fasta(fname, data):
        '''
        Write list of sequences into a fasta file format.
        data =  dataframe or numpy array with columns [id, sequence]
        fname = file to write the information
                .bz2 and .gz extension can be added to filename if compression is desired
        '''
        if type(data) == np.ndarray and data.shape[1]!=2:
            logging.error('write_fasta: data for '+fname+' column count is wrong')
            sys.exit(1)
        elif type(data) == pd.core.frame.DataFrame:
            data = data[['name','sequence']].values
        with fileIO.compress(fname) as f:
            for i in data:
                text = '>' + str(i[0]) + '\n'
                text+= str(i[1]) + '\n'
                f.write(text)
    
    def read_fastx(fname):
        '''
        Reads data from fasta or fastq
        returns filename
        '''
        # write data to fastq file
        try:
            return fileIO.read_fastq(fname)
        except:
            logging.info('failed to read '+fname+' as fastq')
        
        try:
            return fileIO.read_fasta(fname)
        except:
            logging.info('failed to read '+fname+' as fasta')
        
        # throws an error to do
        sys.exit(1)
        
    def write_fastx(prefix, data):
        '''
        Writes data to fasta or fastq
        returns filename
        '''
        # write data to fastq file
        try:
            fname = prefix+'.fq'
            fileIO.write_fastq(fname, data)
            return fname
        except:
            logging.info('failed to write fastq')
        
        # write data to fasta file
        try:
            fname = prefix+'.fa'
            fileIO.write_fasta(fname, data)
            return fname
        except:
            logging.info('failed to write fasta')
        
        # throws an error to do
        sys.exit(1)

    def parse_SAM(fname):
        '''
        This function parses the information from a SAM file and returns a dataframe
        fname = filename of .sam file
        returns dataframe with columns
        ['query_id','database_id','orientation','flag','t_start','t_end',
        'mapq','CIGAR','AS','XS','XN','XM','XO','XG','NM']
        '''
        fsize = os.path.getsize(fname)
        record = False
        col1 = ['query_id','database_id','orientation','flag','t_start','t_end','mapq','CIGAR']
        col2 = ['AS','XS','XN','XM','XO','XG','NM']
        key = {col2[i]:i for i in range(0,len(col2))}
        data = []
        with open(fname,'r') as f:
            while fsize > f.tell():
                text = f.readline()
                # dont record anything until we read past the header
                if ~record and text[:3]=='@PG':
                    record = True
                    text = f.readline()
                # start recording data
                if record:
                    out, info = fileIO.get_SAM_info(text, key)
                    if out[3][-2] == '1':
                        text = f.readline()
                        out2, info2 = fileIO.get_SAM_info(text, key)
                        out[0]+= ' '+out2[0]
                        out[7]+= ' '+out2[7]
                        info = info+info2
                    data.append(out+[j for j in info])
        return pd.DataFrame(data, columns=col1+col2)
    
    def get_SAM_info(read, key):
        '''
        Function to parse useful info from a line of the SAM file
        read = each line of sam file
        key = columns of data to look for
        '''
        read = read.split('\t')
        qname = read[0]
        flag = '{0:10b}'.format(int(read[1])) # decode the sam flag via binary to see if reverse complement is aligned
        rname = read[2]
        if flag[-5] == '1':
            orientation = '-'
            pos2 = int(read[3])
            pos1 = pos2 + int(read[8])
        else:
            orientation = '+'
            pos1 = int(read[3])
            pos2 = pos1 + int(read[8])
        mapq = int(read[4])
        cigar = read[5]
        # get another sam info
        info = np.array([0]*len(key))
        if rname != '*':
            for i in range(11,len(read)):
                x = read[i].split(':') # parse it backwards
                if x[0] in key.keys():
                    info[key[x[0]]] = x[-1]
        data = [qname, rname, orientation, flag, pos1, pos2, mapq, cigar]
        return data, info

    def parse_PAF(fname):
        '''
        This function parses the information from a PAF file and returns a dataframe
        fname = filename of paf file output from minimap2
        returns a dataframe
        '''
        with open(fname,'r') as f:
            text = f.read().split('\n')

        if len(text[-1])==0:
            text = text[:-1] # removes the last line that is empty    
        data = []

        # information we are parsing
        col1 = ['query_id','q_len','q_start','q_end','orientation',
                'database_id','t_len','t_start','t_end','match','tot','mapq']
        col2 = ['tp','cm','s1','s2','NM','AS','ms','nn','rl','cg']
        key = {col2[i]:i for i in range(0,len(col2))}
        for line in text:
            out = line.split('\t')
            # parses for other info
            info = [0]*len(key)
            for i in range(len(col1),len(out)):
                x = out[i].split(':')
                if x[0] in key:
                    info[key[x[0]]] = x[-1]
            # save the data
            data.append(out[:len(col1)]+info)
        data = pd.DataFrame(data, columns=col1+col2)
        data = data.rename(columns={'cg':'CIGAR'})
        return data
