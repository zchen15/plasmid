#!/usr/bin/env python
# This is a python library for editing and plotting genbank files

# computing libraries
import numpy as np
import pandas as pd

# biopython imports
import Bio
import Bio.Seq
import Bio.SeqIO

# system imports
import datetime
import time
import logging
import os

def get_codons():
    # generate codon table
    codons = {'I':['ATT', 'ATC', 'ATA'],
              'L':['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
              'V':['GTT', 'GTC', 'GTA', 'GTG'],
              'F':['TTT', 'TTC'],
              'M':['ATG'],
              'C':['TGT', 'TGC'],
              'A':['GCT', 'GCC', 'GCA', 'GCG'],
              'G':['GGT', 'GGC', 'GGA', 'GGG'],
              'P':['CCT', 'CCC', 'CCA', 'CCG'],
              'T':['ACT', 'ACC', 'ACA', 'ACG'],
              'S':['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
              'Y':['TAT', 'TAC'],
              'W':['TGG'],
              'Q':['CAA', 'CAG'],
              'N':['AAT', 'AAC'],
              'H':['CAT', 'CAC'],
              'E':['GAA', 'GAG'],
              'D':['GAT', 'GAC'],
              'K':['AAA', 'AAG'],
              'R':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
              '*':['TAA', 'TAG', 'TGA']}
    xcodon = []
    for k in codons.keys():
        if k!='*':
            xcodon+= codons[k]
    codons['X'] = xcodon
    return codons

codons = get_codons()

def get_AAlib(AAseq):
    '''
    Define a library of codons for each amino acid
    '''
    catalog = []
    codons = get_codons()
    for j in AAseq.upper():
        catalog.append(codons[j])
    return catalog

def load_args(args, params=None):
    '''
    Writes values from args into params dictionary
    argnew = new parameters
    params = old parameter set
    return updated parameter set as a dictionary
    '''
    if params==None:
        params = {}
    if args!=None:
        argvars = args.__dict__
        for k in argvars.keys():
            if argvars[k]!=None:
                params[k] = argvars[k]
    return params

def read_to_df(fname):
    '''
    Function used to read csv, genbank, fasta, or text file of sequences
    fname = filename to read
    return pandas dataframe of sequences in the file
    '''
    out = []
    if type(fname) == str:
        try:
            print('Try reading',fname,'as genbank')
            rec = Bio.SeqIO.read(fname, 'genbank')
            out = [[fname, fname, str(rec.seq)]]
        except:
            print('Failed to read',fname,'as genbank')
            
            try:
                print('Try reading',fname,'as csv')
                out = read_csv(fname)
            except:
                print('Failed to read',fname,'as csv')

                try:
                    print('Try reading',fname,'as fasta')
                    out = read_fasta(fname)
                except:
                    print('Failed to read',fname,'as fasta')

    # format into a pandas dataframe
    out = pd.DataFrame(out, columns=['filename','name','sequence'])
    return out

def read_csv(fname):
    '''
    Read sequences from a csv file
    fname = file to open
    returns a list with [[filename, name, sequence]]
    '''
    df = pd.read_csv(fname)
    # figure out which has name
    name = [fname]*len(df)
    col = ['name','label','Strand','locus_tag']
    for c in col:
        if c in df.columns:
            name = df[c].values
            break
    # figure out which has sequence
    col = ['Sequence','sequence','seq']
    for c in col:
        if c in df.columns:
            seq = df[c].values
            break
    # assign filename 
    if ('filename' in df.columns) == False:
        df['filename'] = fname
    filename = df.filename
    return [[filename[i], name[i], seq[i]] for i in range(len(df))]

def read_fasta(fname):
    '''
    Read sequences from a fasta file
    fname = file to open
    returns a list with [[filename, name, sequence]]
    '''
    with open(fname, 'r') as f: 
        rec = list(Bio.SeqIO.parse(f, 'fasta'))
    df = [[fname, i.name, str(i.seq)] for i in rec]
    if len(df) > 0:
        return df
    else:
        sys.error()

def revcomp(seq):
    '''
    Returns reverse complement of a nucleic acid string
    '''
    return str(Bio.Seq.Seq(seq).reverse_complement())

def complement(seq):
    '''
    Return complement of a nucleic acid string
    '''
    return str(Bio.Seq.Seq(seq).complement())

def format_to_dna(sequence):
    '''
    Convert a string to have only IUPAC dna codes. Other erroneous letters are removed.
    sequence = string to change to IUPAC only codes
    return the cleaned string
    '''
    allowed = 'ATGCRYSWKMBDHVN'
    out = str(sequence).upper()
    # convert U to T
    out = out.replace('U','T')
    # convert gap characters to N
    out = out.replace('.','N')
    out = out.replace('-','N')
    # remove spaces
    out = out.replace(' ','')
    # remove erroneous characters
    remove = set(out) - set(allowed)
    for i in remove:
        out = out.replace(i,'')
    return out

def isDNA(seq):
    return set(seq.upper()).issubset(set('ATGCU'))

def isAA(seq):
    return set(seq).issubset(set(codons))

def search_key(data, keyword, rev=False):
    '''
    Searches dictionary for first instance of keyword
    data = dictionary to search for the keyword
    keyword = keyword value to find in keys
    rev = search key in keyword instead of keyword in key
    return value of the keyword if found. Otherwise return None
    '''
    for key in data.keys():
        if rev:
            value = (key.lower() in keyword.lower())
        else:
            value = (keyword.lower() in key.lower())
        if value:
            return data[key]
    return None

def get_rand_letters(N, method=0):
    '''
    Generate nucleotide domain of length N
    N = length of nucleotides to generate
    method = method to use
    returns string of nucleotides
    '''
    # no repeating bases
    m0 = {'A':'TGC', 'T':'GCA', 'G':'ATC', 'C':'ATG'}
    # opposite day
    m1 = {'A':'GC', 'T':'GC', 'G':'AT', 'C':'AT'}
    # loop back on predecessor bases
    m2 = {'AT':'CA', 'AG':'CA', 'AC':'GA',
          'TA':'GT', 'TG':'CT', 'TC':'GT',
          'GA':'TG', 'GT':'AG', 'GC':'AG',
          'CA':'TC', 'CT':'AC', 'CG':'TC'}
    # loop back on predecessor bases
    m3 = {'AT':'CA', 'AG':'CA', 'AC':'GA',
          'TA':'GT', 'TG':'CT', 'TC':'GT',
          'GA':'TG', 'GT':'AG', 'GC':'TG',
          'CA':'TC', 'CT':'AC', 'CG':'AC'}
    # loop back on predecessor bases
    m4 = {'AT':'GA', 'AG':'CA', 'AC':'GA',
          'TA':'CT', 'TG':'CT', 'TC':'GT',
          'GA':'TG', 'GT':'AG', 'GC':'AG',
          'CA':'TC', 'CT':'AC', 'CG':'TC'}
    # loop back on predecessor bases
    m5 = {'AT':'GA', 'AG':'CA', 'AC':'GA',
          'TA':'CT', 'TG':'CT', 'TC':'GT',
          'GA':'TG', 'GT':'AG', 'GC':'TG',
          'CA':'TC', 'CT':'AC', 'CG':'AC'}
    letters = [m0, m1, m2, m3, m4, m5]
    letters = letters[method]
    out = list(letters.keys())
    out = out[np.random.randint(0,len(out))]
    N_in = len(out)
    N_rand = len(letters[out])
    for j in np.random.randint(0, N_rand, N-len(out)):
        key = out[-N_in:]
        out+= letters[key][j]
    return out

def get_timestamp():
    return '._'+datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')

def set_log_level(filename='run.log', level='INFO'):
    fmt = 'pid['+str(os.getpid())+'] %(asctime)s.%(msecs)03d %(levelname)s: %(message)s'
    dfmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(fmt, datefmt=dfmt)
    # print to stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(level)
    stdout_handler.setFormatter(formatter)
    logging.getLogger().addHandler(stdout_handler)
    # print to log file
    handle = logging.FileHandler(filename=filename)
    handle.setLevel(level)
    handle.setFormatter(formatter)
    logging.getLogger().addHandler(handle)
    logging.getLogger().setLevel(level) 

def get_memory_usage():
    pid = os.getpid()
    process = psutil.Process(pid)
    return process.memory_full_info().rss/1024/1024

def generate_biobricks(database, keyword, verbose=True):
    '''
    Searches an annotation database for sequences matching the keyword
    database = list or andas dataframe with columns [label, feature_type, color, strand, sequence]
    keyword = keyword matching the label label in the annotation database
    verbose = print list of parts that were generated
    return a list of plasmid objects which can be inserted into another plasmid or replace existing features
    '''
    output = []
    data = database[['locus_tag','feature_type','color','strand','sequence']].values
    for [label, feature_type, color, strand, sequence] in data:
        if (keyword.lower() in label.lower()):
            loc = Bio.SeqFeature.FeatureLocation(0, len(sequence), int(strand))
            qualifiers = {'locus_tag':[label], 'color':color}
            gene = Bio.SeqFeature.SeqFeature(type=feature_type, location=loc, qualifiers=qualifiers)
            SeqRecord = Bio.SeqRecord.SeqRecord(seq=sequence, features=[gene])
            output.append(plasmid(SeqRecord))
    # print results
    if verbose and len(output)>0:
        df = pd.concat([i.to_dataframe() for i in output])
        df = df.drop(columns=['sequence','location','start','end']).reset_index(drop=True)
        print(df)
    return output

def unpack_FeatureLocation(loc):
    '''
    Unpack locations nested in locations object
    '''
    if 'parts' in loc.__dict__.keys():
        out = []
        for p in loc.parts:
            out+= Plasmid.unpack_FeatureLocation(p)
        return out
    else:
        return [loc]

