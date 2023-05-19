#!/usr/bin/env python
# This is a python library for editing and plotting genbank files

# computing libraries
import numpy as np
import pandas as pd

# biopython imports
import Bio
import Bio.Seq
import Bio.SeqIO
import json

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

codon_table = get_codons()

def get_codon_synonymous():
    # generate table of synonymous codons
    codons = get_codons()
    out = {}
    for k in codons.keys():
        if k!='X':
            x = codons[k]
            if len(x) > 1:
                x = np.array(x)
                for i in x:
                    out[i] = [str(j) for j in x[x!=i]]
            else:
                out[x[0]] = [x[0]]
    return out

synonymous_codon_table = get_codon_synonymous()

def get_IUPAC_nt():
    return 'ATGCRYSWKMBDHVN'

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

def translate(seq, frame=0, table='Standard'):
    '''
    Translates a nucleotide sequence into amino acids
    seq = input nucleotide sequence
    frame = reading frame
    table = codon table to use
    return amino acid sequence
    '''
    seq = Bio.Seq.Seq(seq)
    return str(seq[frame:].translate(table=table))

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
    allowed = get_IUPAC_nt()
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
    dna = get_IUPAC_nt()
    return set(seq.upper()).issubset(set(dna))

def isAA(seq):
    return set(seq).issubset(set(codon_table))

def isBool(value):
    '''
    Check if list of a values are only booleans
    '''
    values = np.unique(value)
    # splice based on boolean array
    if len(values) < 3 and type(values[0])==np.bool_ and type(values[-1])==np.bool_:
        return True
    else:
        return False

def isStringArray(value):
    '''
    Check if list contains only strings
    '''
    x = [type(i)==str for i in value]
    return sum(x) == len(value)

def isIntArray(value):
    '''
    Check if list contains only integers
    '''
    x = [type(i)==int for i in value]
    return sum(x) == len(value)

def substring_search(ref, qry):
    """
    Searches for all occurrences of the given substrings in a string
    ref = string to search
    qry = substring to locate
    return a list with start locations
    """
    out = []
    x = ref
    offset = 0
    while len(x) > 0:
        idx = x.find(qry)
        if idx == -1:
            break
        else:
            out.append(idx+offset)
            x = x[idx+1:]
            offset+= idx+1
    return out

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

def json_pretty(data):
    return json.dumps(data, indent=4)

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

def load_annotation_database(filename, lib_format='csv'):
    '''
    Loads an annotation into a pandas dataframe, so it can be used to annotate sequences
    lib_format = type of formatting the annotation database is stored as.
                 valid options are csv or ApE
    returns a pandas dataframe with columns [label, feature, color, start, end, strand, sequence]
    '''
    if lib_format=='csv':
        df = pd.read_csv(filename)
        if ('strand' in df.columns)==False:
            print('Sequence orientation is not defined. Please add strand column to csv file')
            return None
        df['strand'] = df['strand'].astype(int)
    elif lib_format=='ApE':
        df = pd.read_csv('Default_Features.txt', delimiter='\t', header=None)
        df = df.drop(columns=df.columns[4:])
        df.columns = ['locus_tag','sequence','feature_type','color']
        df['strand'] = 1
    # do some safety checks on the dna sequences
    data = []
    for i in range(0,len(df)):
        seq = format_to_dna(df.iloc[i]['sequence'])
        # correct strand orientation to 5'->3'
        if df.iloc[i]['strand'] < 0:
            seq = revcomp(seq)
        data.append(seq)
    df['sequence'] = data    
    df['length'] = [len(seq) for seq in df['sequence']]
    df['start'] = None
    df['end'] = None
    df['strand'] = 1
    return df

def unpack_FeatureLocation(loc):
    '''
    Unpack locations nested in locations object
    '''
    if 'parts' in loc.__dict__.keys():
        out = []
        for p in loc.parts:
            out+= unpack_FeatureLocation(p)
        return out
    else:
        return [loc]

def unwrap_FeatureLocation(loc, length):
    '''
    Get proper start and end locations of a genetic feature that wraps around the origin
    loc = Bio.SeqFeature location
    length = length of circular Plasmid
    returns [start, end, strand]
    '''
    strand = loc.strand
    start = loc.start%length
    end = loc.end%length
    if end==0:
        end=length
    return [start, end, strand]

def wrap_FeatureLocation(start, end, strand, length):
    '''
    Create a CompoundLocation or FeatureLocation object that wrap the annotation around the origin of a circular Plasmid
    start = start base position of the gene
    end = end base position of the gene
    strand = reading frame of the gene 5'->3' is +1
    length = base pair length of the Plasmid
    returns Bio.SeqFeature.CompoundLocation or Bio.SeqFeature.FeatureLocation object
    '''
    start = start%length
    end = end%length
    if end==0:
        end = length

    if end < start:
        loc1 = Bio.SeqFeature.FeatureLocation(int(start), int(length), int(strand))
        loc2 = Bio.SeqFeature.FeatureLocation(0, int(end), int(strand))
        loc = Bio.SeqFeature.CompoundLocation([loc1,loc2])
    else:
        loc = Bio.SeqFeature.FeatureLocation(int(start), int(end), int(strand))
    return loc

def shift_CompoundLocation(loc, N, L, circular=True):
    '''
    Applies nucleotide shift to compound locs and merges adjacent FeatureLocations that arise after the shift
    loc = Bio.SeqFeature.FeatureLocation or Bio.SeqFeature.CompoundLocation object
    N = number of nucleotides to shift
    L = length of the Plasmid sequence
    circular = apply relevant corrections if the genome is circular
    returns a Bio.SeqFeature.CompoundLocation or Bio.SeqFeature.FeatureLocation
    '''
    loc+=int(N)
    if circular:
        if 'parts' in loc.__dict__:
            # correct the shift in each sub loc
            out = []
            for p in loc.__dict__['parts']:
                x = shift_CompoundLocation(p, 0, L, circular=True)
                # unwrap joint locations
                if 'parts' in x.__dict__:
                    out+= x.__dict__['parts']
                else:
                    out.append(x)
            # merge adjacent locations
            keep = []
            start = -1
            for i in range(0, len(out)):
                # find the start
                if start == -1:
                    start = out[i].start
                # record the new position when a discontinuity is found
                if i+1 == len(out) or (out[i].end == out[i+1].start and out[i].strand == out[i+1].strand) == False:
                    keep.append(Bio.SeqFeature.FeatureLocation(start, out[i].end, out[i].strand))
                    start = -1
            if len(keep) > 1:
                loc = Bio.SeqFeature.CompoundLocation(keep)
            elif len(keep) == 1:
                loc = keep[0]
        else:
            loc = wrap_FeatureLocation(loc.start, loc.end, loc.strand, L)
    return loc

