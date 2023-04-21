#!/usr/bin/env python
# Class functions to examine genbank files and sequences

# computing libraries
import numpy as np
import pandas as pd
import copy

# reading and writing genbank files
import Bio
import Bio.SeqIO

# local libaries
from .misc import search_key
from .misc import format_to_dna
from .misc import revcomp

from .misc import isDNA
from .misc import isAA
from .misc import isBool
from .misc import isStringArray
from .misc import isIntArray

from .misc import read_to_df
from .misc import get_timestamp
from .misc import unwrap_FeatureLocation
from .misc import wrap_FeatureLocation
from .misc import shift_CompoundLocation

from .graphic import Graphic
from .aligner import Aligner

def read_genbank(fname):
    '''
    Reads a genbank file into plasmid class object
    fname = genbank file to read
    returns Plasmid object with data loaded
    '''
    # load data with SeqIO
    print('reading ',fname,' as genbank file')
    SeqRecord = Bio.SeqIO.read(fname, 'genbank')

    # add some info to the genbank files
    for gene in SeqRecord.features:
        # find the label ordering is the label that takes precedance
        name = ['ApEinfo_label','label','locus_tag','gene']
        label = 'unknown'
        for k in name:
            if k in gene.qualifiers:
                if type(gene.qualifiers[k]) == list:
                    label = gene.qualifiers[k][0]
                elif type(gene.qualifiers[k]) == str:
                    label = gene.qualifiers[k]
        gene.qualifiers['locus_tag'] = [label]

    # correct genbank read issue in parsing features with CompoundLocations which wrap around the origin
    for j in range(0, len(SeqRecord.features)):
        d = SeqRecord.features[j].location.__dict__ 
        if 'parts' in d and len(d['parts'])==2:
            # if ordering is wrong, swap it
            if d['parts'][0].start==0 and d['parts'][1].end==len(SeqRecord.seq) and d['parts'][0].strand == d['parts'][1].strand:
                SeqRecord.features[j].location.__dict__['parts'] = [d['parts'][1], d['parts'][0]]

    # create the Plasmid dataframe
    return Plasmid(SeqRecord)

def read_fasta(fname):
    '''
    Read sequences from csv or fasta files
    fname = file to read
    returns list of Plasmid objects with data loaded
    '''
    df = read_to_df(fname)
    out = []
    # generate list of plasmids
    for name, seq in df[['name','sequence']].values: # todo test debug
        x = Plasmid(seq)
        x = x.annotate(name=name, sequence=seq)
        out.append(x)
    return out

class Plasmid:
    '''
    This class holds data and functions for manipulating sequence records read from a genbank file
    '''
    # default color palette
    params = {}
    # default colors for fluorescent proteins
    default_annotations = {'molecule_type':'DNA',
                           'topology':'circular'}

    def __init__(self, seq):
        '''
        Initialize some data into the Plasmid dataframe
        seq = Bio.SeqRecord or string
        '''
        # initialize the Bio.SeqRecord object
        try:
            if type(seq) == str:
                self.SeqRecord = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(seq))
            else:
                self.SeqRecord = seq
            self.features = self.SeqRecord.features
        except:
            print('error object is not a Bio.SeqRecord or String')
        
        # set default cmap
        features = ['rbs','cds','promoter','misc_binding','ncRNA','terminator','primer_bind','unknown','rep_origin','protein_bind']
        palette = 'Category20_20'
        self.params['colors'] = {}
        self.params['colors']['features'] = Graphic.get_colormap(features, palette)
        # set default color for known or unknown feature or label
        self.params['colors']['locus_tag'] = {'bfp':'blue','cfp':'dodgerblue','gfp':'green','yfp':'gold','rfp':'red','mirfp':'purple'}

        # update the feature list ordering
        self.reset_colors(ow_type=False, ow_locus=False, inplace=True)
        self.update(inplace=True)

    def reset_colors(self, ow_type=True, ow_locus=True, inplace=False):
        '''
        Apply default or custom color map to features or labels
        ow_type = overwrite existing colors for features
        ow_locus = overwrite existing colors for locus_tags
        inplace = apply inplace
        return Plasmid with new colors
        '''
        label_cmap = self.params['colors']['locus_tag']
        ftype_cmap = self.params['colors']['features']
        default_color = self.params['colors']['features']['unknown']
        out = self.inplace(inplace=inplace)
        # update colors in SeqRecord
        for gene in out.SeqRecord.features:
            # check for label color
            label_color = search_key(label_cmap, gene.qualifiers['locus_tag'][0], rev=True)
            # check for feature color
            feature_color = search_key(ftype_cmap, gene.type)
            # check for saved colors
            fwd_color = search_key(gene.qualifiers, 'ApEinfo_fwdcolor')
            rev_color = search_key(gene.qualifiers, 'ApEinfo_revcolor')
            # use label_color if available
            if (ow_locus or (fwd_color==None and rev_color==None)) and label_color!=None:
                Plasmid.set_color(gene, label_color)
            # use feature color if available
            elif (ow_type or (fwd_color==None and rev_color==None)) and feature_color!=None:
                Plasmid.set_color(gene, feature_color)
            # apply default colors if nothing found
            elif fwd_color==None or rev_color==None:
                Plasmid.set_color(gene, default_color)
        return out

    def set_color(gene, color):
        # sets the color of the gene
        if type(color) == list:
            gene.qualifiers['ApEinfo_fwdcolor'] = [color[0]]
            gene.qualifiers['ApEinfo_revcolor'] = [color[-1]]
        elif type(color) == str:
            gene.qualifiers['ApEinfo_fwdcolor'] = [color]
            gene.qualifiers['ApEinfo_revcolor'] = [color]

    def set_column(self, key, col, value):
        '''
        Sets columns of Plasmid dataframe
        '''
        key = self.unpack_key(key)
        if type(value)!= list:
            value = [value]*len(key)
        for i in range(len(key)):
            k = key[i]
            val = value[i]
            # set color values of each feature
            if col in ['color','fwd_color','rev_color']:
                if col=='color':
                    strand = self.SeqRecord.features[k].location.strand
                elif col=='fwd_color':
                    strand = 1
                else:
                    strand = -1
                if strand==1:
                    self.SeqRecord.features[k].qualifiers['ApEinfo_fwdcolor'] = [val]
                else:
                    self.SeqRecord.features[k].qualifiers['ApEinfo_revcolor'] = [val]
            elif col=='location':
                self.SeqRecord.features[k].location = val
            elif col=='locus_tag':
                self.SeqRecord.features[k].qualifiers['locus_tag'] = [val]
            elif col=='type':
                self.SeqRecord.features[k].type = val
            else:
                sys.error()
        self.update(inplace=True)

    def unpack_key(self, key):
        '''
        Formats key into a list of values
        key = slice, range, bool array, int array, or list
        Returns a list of values
        '''
        if type(key) == int:
            return [key]
        elif type(key) == slice or type(key)==range:
            start = key.start
            stop = key.stop
            if start==None:
                start = 0
            if stop==None:
                stop = len(self.SeqRecord.features) 
            return list(range(start, stop))
        elif isIntArray(key):
            return [int(k) for k in key]
        elif isBool(key):
            key = np.arange(len(key))[key]
            return [int(k) for k in key]
        elif type(key) == list:
            out = []
            for k in key:
                out+= self.unpack_key(k)
            return out

    def update(self, inplace=False):
        '''
        Synchronize entries from SeqRecord to Plasmid parts list
        Run sort_features and sort gene features by starting base index
        '''
        if inplace:
            out = self
        else:
            out = self.copy()
        out.sort_features(inplace=True)
        # sync gene labels for backward compatibility with ApE
        for gene in out.SeqRecord.features:
            label = gene.qualifiers['locus_tag'][0]
            gene.qualifiers['ApEinfo_label'] = [label]
            # delete useless tags
            tags = ['color','label','ApEinfo_graphicformat']
            for t in tags:
                if t in gene.qualifiers:
                    del(gene.qualifiers[t])
        out.check_annotations()
        return out
    
    def sort_features(self, inplace=False):
        '''
        Sorts the features in the SeqRecord object
        returns a modified copy of the object
        '''
        if inplace:
            out = self
        else:
            out = self.copy()
        # do sorting with pandas
        df = out.to_dataframe()
        idx = df.sort_values(by=['start','end','locus_tag','type','color']).index
        out.SeqRecord.features = [out.SeqRecord.features[i] for i in idx]
        return out

    def check_annotations(self):
        '''
        Check features for consistency against the sequence length
        returns false if there are issues with the gene annotation
        '''
        # check seq length consistency
        seq_length = len(self.SeqRecord.seq)
        good = True
        for gene in self.SeqRecord.features:
            start = gene.location.start
            end = gene.location.end
            if start > seq_length:
                print('Warning: annotation issues with '+gene.qualifiers['locus_tag'][0]+' start < 0')
                good = False
            if end > seq_length:
                print('Warning: annotation issues with '+gene.qualifiers['locus_tag'][0]+' end > seq_length')
                good = False
        
        # check header consistency
        for k in self.default_annotations.keys():
            if ~(k in self.SeqRecord.annotations.keys()):
                self.SeqRecord.annotations[k] = self.default_annotations[k]
        
        if good==False:
            sys.error()
        return good

    def to_genbank(self, filename, timestamp=False):
        '''
        Writes the Plasmid SeqRecord to a filename
        filename = filename of output file
        timestamp = adds a timestamp to the output file
        '''
        if timestamp:
            filename+= get_timestamp()
        # sync info and update it
        self.update(inplace=True)
        with open(filename,'w') as f:
            Bio.SeqIO.write(self.SeqRecord, f, 'genbank')
    
    def copy(self):
        '''
        Returns a copy of self
        '''
        return copy.deepcopy(self)

    def inplace(self, inplace):
        '''
        Decides if a copy of self should be returned
        '''
        if inplace:
            return self
        else:
            return self.copy()

    def get_best_origin(self):
        '''
        Find best index to set the origin of the Plasmid
        return index to set the origin
        '''
        genes = self.SeqRecord.features
        for i in range(0, len(genes)):
            [s1, s2, strand] = unwrap_FeatureLocation(genes[i].location, len(self.SeqRecord.seq))
            for k in [s1, s2]:
                good = True
                for j in range(i, len(genes)):
                    [start, end, strand] = unwrap_FeatureLocation(genes[j].location, len(self.SeqRecord.seq))
                    # break out of loop if origin is not in good location
                    if (end > start and k > start and k < end) or (start > end and (k > start or k < end)):
                        good = False
                        break
                if good:
                    return k

    def set_origin(self, index=None, inplace=False):
        '''
        Sets a new position for the origin coordinates of a circular Plasmid
        index = base index to set new origin. If index==None, then origin is set at the first contiguous genetic feature
        inplace = performs operation inplace or returns a copy with the modified data
        returns the Plasmid dataframe
        '''
        out = self.inplace(inplace)
        # find best index to set the origin
        if index==None:
            index = self.get_best_origin()
        # do nothing if better index is not found 
        if index==None:
            return out
        elif index > len(self.SeqRecord.seq):
            index = index - len(self.SeqRecord.seq)
        elif index < 0:
            index = len(self.SeqRecord.seq) + index
        out.SeqRecord.seq = Bio.Seq.Seq(str(out.SeqRecord.seq[index:]) + str(out.SeqRecord.seq[0:index]))
        out.shift_features(-index, circular=True, inplace=True)
        out.update(inplace=True)
        return out

    def to_graphic(self):
        rec = Graphic(data=self)
        return rec
 
    def reverse_complement(self, inplace=False):
        '''
        Returns reverse complement of the Plasmid
        '''
        if inplace:
            out = self
        else:
            out = self.copy()
        out.SeqRecord = out.SeqRecord.reverse_complement()
        return out 

    def insert_sequence(self, index, seq, disrupt=True, inplace=False):
        '''
        Insert sequences into the Plasmid
        index = position to insert the sequence or Plasmid construct
        seq = sequence or linearized Plasmid construct to insert into the Plasmid
        disrupt = if sequence is inserted into an existing feature, break the feature into two subsequences
        inplace = modifications are done inplace
        return the modified Plasmid
        '''
        # format to integer
        index = int(index)
        if inplace:
            out = self
        else:
            out = self.copy()
        # wrap location 
        if index < 0:
            index = len(out.SeqRecord.seq) - index
        elif index > len(out.SeqRecord.seq):
            index = index - len(out.SeqRecord.seq)
        # format the values
        if type(seq)==str:
            sequence = seq
            insert_features = []
        elif type(seq).__name__=='Plasmid':
            sequence = str(seq.SeqRecord.seq)
            insert_features = seq.shift_features(index, inplace=False).copy().SeqRecord.features
        else:
            print('insert_sequence: no sequence provided')
            return out
        # modify the sequence
        out.SeqRecord.seq = Bio.Seq.Seq(str(out.SeqRecord.seq[:index]) + sequence + str(out.SeqRecord.seq[index:]))
        # shift features around
        for j in range(0,len(out.SeqRecord.features)):
            # handle compound locations
            if 'parts' in out.SeqRecord.features[j].location.__dict__:
                f = out.SeqRecord.features[j].location.__dict__['parts']
                for i in range(0,len(f)):
                    start = f[i].start
                    end = f[i].end
                    strand = f[i].strand
                    if start > index:
                        start+=len(sequence)
                    if end > index:
                        end+=len(sequence)
                    # special condition for genes which wrap around the origin
                    c2 = i+1 < len(f) and (end == len(out.SeqRecord.seq)-len(sequence) and f[i+1].start == 0)
                    if c2:
                        end+=len(sequence)
                    out.SeqRecord.features[j].location.__dict__['parts'][i] = Bio.SeqFeature.FeatureLocation(start, end, strand)
            # handle single locations
            else:
                start = out.SeqRecord.features[j].location.start
                end = out.SeqRecord.features[j].location.end
                strand = out.SeqRecord.features[j].location.strand
                if start >= index:
                    start+=len(sequence)
                if end > index:
                    end+=len(sequence)
                out.SeqRecord.features[j].location = Bio.SeqFeature.FeatureLocation(start, end, strand)

        # record disruptions to existing features
        keep = []
        if disrupt:
            for j in range(0,len(out.SeqRecord.features)):
                # handle compound locations
                if 'parts' in out.SeqRecord.features[j].location.__dict__:
                    f = out.SeqRecord.features[j].location.__dict__['parts']
                    # partition sublocation to before or after the insert
                    isDisrupted = False
                    new = []
                    for i in range(0,len(f)):
                        start = f[i].start
                        end = f[i].end
                        strand = f[i].strand
                        # split the subsequences, detect that genes wrapping around origin were disrupted
                        c1 = (start <= index) and (end > index + len(sequence))
                        c2 = (start < index) and (end >= index + len(sequence))
                        if c1 or c2:
                            isDisrupted = True
                            if index > start: # special logic 
                                new.append(Bio.SeqFeature.FeatureLocation(start, index, strand))
                            if end > index + len(sequence):
                                new.append(Bio.SeqFeature.FeatureLocation(index + len(sequence), end, strand))
                        else:
                            new.append(f[i])
                    # flag that the gene was split
                    out.SeqRecord.features[j].location.__dict__['parts'] = new
                    if isDisrupted: out.SeqRecord.features[j].qualifiers['locus_tag'][0]+=' disrupted'
                # handle single locations
                else:
                    start = out.SeqRecord.features[j].location.start
                    end = out.SeqRecord.features[j].location.end
                    strand = out.SeqRecord.features[j].location.strand
                    # split the subsequences
                    if start < index and end > index + len(sequence):
                        loc1 = Bio.SeqFeature.FeatureLocation(start, index, strand)
                        loc2 = Bio.SeqFeature.FeatureLocation(index + len(sequence), end, strand)
                        out.SeqRecord.features[j].location = Bio.SeqFeature.CompoundLocation([loc1, loc2])
                        out.SeqRecord.features[j].qualifiers['locus_tag'][0]+=' disrupted'
        # record the disruption
        out.SeqRecord.features+=insert_features
        out.update(inplace=True)
        return out

    def delete_sequence(self, key, disrupt=True, inplace=False):
        '''
        Delete sequences on the Plasmid
        key = range(start, end) base pair positions to delete.
        disrupt = if a sequence is deleted from a feature, flag the feature which had the deletion
        return copy of the modified Plasmid
        '''
        out = self.inplace(inplace)
        d1 = key.start
        d2 = key.stop
        dlen = d2-d1
        if dlen < 0:
            out.delete_sequence(range(d1, len(out.SeqRecord.seq)), disrupt=disrupt, inplace=True)
            out.delete_sequence(range(0, d2), disrupt=disrupt, inplace=True)
        elif dlen > 0:
            out.SeqRecord.seq = Bio.Seq.Seq(str(out.SeqRecord.seq[:d1]) + str(out.SeqRecord.seq[d2:]))
            # shift locations around
            for j in range(0,len(out.SeqRecord.features)):
                isDisrupted = False
                # handle compound location
                if 'parts' in out.SeqRecord.features[j].location.__dict__:
                    f = out.SeqRecord.features[j].location.__dict__['parts']
                    remove = []
                    for i in range(0,len(f)):
                        s1 = f[i].start
                        s2 = f[i].end
                        strand = f[i].strand
                        # genes above the area need to shift location
                        if (s1 >= d2 and s2 > d2)  :
                            f[i] = Bio.SeqFeature.FeatureLocation(s1 - dlen, s2 - dlen, strand)
                        # genes caught in the rear half are truncated
                        elif (s1 < d1 and s2 > d1  and s2 < d2):
                            f[i] = Bio.SeqFeature.FeatureLocation(s1, d1, strand)
                            isDisrupted = True
                        # genes caught in front half are truncated
                        elif (s1 < d2 and s1 > d1  and s2 > d2):
                            f[i] = Bio.SeqFeature.FeatureLocation(d1, s2 - dlen, strand)
                            isDisrupted = True
                        # gene overlapping the deletion area truncated in middle
                        elif (s1 < d1 and s2 > d2):
                            isDisrupted = True
                            f[i] = Bio.SeqFeature.FeatureLocation(s1, s2 - dlen, strand)
                        # genes caught in the deletion area are deleted
                        elif s1 >= d1 and s2 <= d2:
                            isDisrupted = True
                            remove.append(i)
                    keep = [f[i] for i in set(range(0,len(f))) - set(remove)]
                    if len(keep) > 1:
                        out.SeqRecord.features[j].location = Bio.SeqFeature.CompoundLocation(keep)
                    elif len(keep) == 1:
                        out.SeqRecord.features[j].location = keep[0]
                    else:
                        out.SeqRecord.features[j].location = None
                # handle single locations
                else:
                    s1 = out.SeqRecord.features[j].location.start
                    s2 = out.SeqRecord.features[j].location.end
                    strand = out.SeqRecord.features[j].location.strand
                    # genes above the area need to shift location
                    if (s1 >= d2 and s2 > d2):
                        out.SeqRecord.features[j].location = Bio.SeqFeature.FeatureLocation(s1 - dlen, s2 - dlen, strand)
                    # genes caught in the rear half are truncated
                    elif (s1 < d1 and s2 > d1  and s2 <= d2):
                        out.SeqRecord.features[j].location = Bio.SeqFeature.FeatureLocation(s1, d1, strand)
                        isDisrupted = True
                    # genes caught in front half are truncated
                    elif (s1 < d2 and s1 >= d1  and s2 > d2):
                        out.SeqRecord.features[j].location = Bio.SeqFeature.FeatureLocation(d1, s2 - dlen, strand)
                        isDisrupted = True
                    # gene overlapping the deletion area truncated in middle
                    elif (s1 < d1 and s2 > d2):
                        isDisrupted = True
                        out.SeqRecord.features[j].location = Bio.SeqFeature.FeatureLocation(s1, s2 - dlen, strand)
                    # genes caught in the deletion area are deleted
                    elif s1 >= d1 and s2 <= d2:
                        isDisrupted = True
                        out.SeqRecord.features[j].location = None
                if isDisrupted and disrupt:
                    out.SeqRecord.features[j].qualifiers['locus_tag'][0]+=' subsequence'
            # purge genes without valid locations
            j=0
            while j < len(out.SeqRecord.features):
                if out.SeqRecord.features[j].location==None:
                    out.SeqRecord.features.pop(j)
                else:
                    j+=1
            out.update(inplace=True)
        return out

    def replace_sequence(self, key, value, disrupt=True, inplace=False):
        '''
        Replace a sequence at a given location with a new sequence
        key = range(start,end) base pair positions on the Plasmid to replace the sequence
        value = a Plasmid dataframe, which contains one sequence and features associated with this sequence
        disrupt = disrupt existing genetic features. Meaning if a gene is caught in the deletion zone, it will be truncated.
        inplace = perform modifications inplace
        '''
        out = self.inplace(inplace)        
        # delete the sequence
        out.delete_sequence(key, disrupt=disrupt, inplace=True)
        # insert the new sequence
        out.insert_sequence(key.start, value, disrupt=disrupt, inplace=True)
        return out

    def splice_sequence(self, key, inplace=False):
        '''
        splice a sequence from the Plasmid
        key = range(start, end) or slice(start, end) base pair position to cut Plasmid
        '''
        out = self.inplace(inplace)
        start = key.start
        stop = key.stop
        length = len(out.SeqRecord.seq)
        # trim the sequence
        if stop < start:
            out.set_origin(stop, inplace=True)
            out.delete_sequence(range(0, start-stop), disrupt=True, inplace=True)
        else:
            out.delete_sequence(range(stop, length), disrupt=True, inplace=True)
            out.delete_sequence(range(0, start), disrupt=True, inplace=True)
        return out

    def splice(self, inplace=False):
        '''
        Get a slice sequence from start of all features to end of all features
        '''
        gene = self.SeqRecord.features[0]
        [start, end, strand] = unwrap_FeatureLocation(gene.location, len(self.SeqRecord.seq))
        # find starting index
        for i in range(1, len(self.SeqRecord.features)):
            gene = self.SeqRecord.features[i]
            [s1, s2, strand] = unwrap_FeatureLocation(gene.location, len(self.SeqRecord.seq))
            c1 = s2 > s1
            c2 = end > start
            c3 = s2 > end
            c4 = s2 > start
            # exit if full Plasmid is covered
            if ~c1 & ~c2 & c4:
                return self.inplace(inplace)
            elif (c1 & c2 & c3) or (~c1 & c2) or (~c1 & ~c2 & c3):
                end = s2
        return self.splice_sequence(range(start, end), inplace)
    
    def __getitem__(self, key):
        '''
        Splice Plasmid sequence with range or filter through features on the Plasmid
        key = if key is a range or slice, then slice using base pair index from start to end
              if key is a str, then return the column values of pandas dataframe
              if key is an int, then return feature at the feature index
              if key is a list of integers, then return features only at the indexes in the list 
              if key is a boolean array, then return features only at indexes which evaluate True 
        '''
        out = self.copy()
        # splice based on slice or range
        if type(key)==range:
            out.splice_sequence(key, inplace=True)
        elif type(key)==slice:
            start = key.start
            stop = key.stop
            if start==None:
                start = 0
            if stop==None:
                stop = len(out.SeqRecord.seq)
            out.splice_sequence(range(start, stop), inplace=True)
        # splice based on feature index
        elif type(key)==int:
            out.SeqRecord.features = [out.SeqRecord.features[key]]
        # returns column entries from pandas dataframe
        elif type(key) == str or isStringArray(key):
            return self.to_dataframe()[key]
        elif isIntArray(key):
            out.SeqRecord.features = [out.SeqRecord.features[i] for i in key]
        # slice based on boolean array
        elif isBool(key):
            keep = np.arange(len(key))[key]
            out.SeqRecord.features = [out.SeqRecord.features[i] for i in keep]
        # work on output from pandas dataframe
        elif type(key) == pd.core.series.Series:
            return out[key.values]
        # return key and column
        elif type(key)==tuple:
            for k in key:
                out = out[k]
            return out
        out.update(inplace=True) 
        return out

    def __setitem__(self, key, value):
        '''
        Sets a feature to something
        key = 2D index of [<index>, <column>] or [<index>, <column>] or 1D [<index>] of features or sequences to replace
        value = slice of a Plasmid dataframe, which contains one sequence and features associated with this position
                the sequence at the location of the feature will be replaced with the provided sequence
                the SeqFeature at the index will be replaced with the new feature
                all modifications are done in place
        '''       
        if type(key)==int: # debug if features may get overriden
            [start, stop, strand] = unwrap_FeatureLocation(self.SeqRecord.features[key].location, len(self.SeqRecord.seq))
            self.replace_sequence(range(start, stop), value, disrupt=True, inplace=True)
        elif type(key)==range:
            self.replace_sequence(key, value, disrupt=True, inplace=True)
        elif type(key)==slice:
            start = key.start
            stop = key.stop
            if start==None:
                start = 0
            if stop==None:
                stop = len(out.SeqRecord.seq)
            self.replace_sequence(range(start, stop), value, disrupt=True, inplace=True)
        elif type(key)==tuple:
            if type(key[1])==str:
                self.set_column(key[0], key[1], value)
            elif type(key[0])==str:
                self.set_column(key[1], key[0], value)
        else:
            sys.error()

    def __delitem__(self, key):
        '''
        Delete feature, slice of features, or range of sequences
        key = if key is a slice, then features from start,end,stride are delete
              if key is an int, then the feature at the index is deleted
              if key is a range, then sequences from bp index start to end is deleted
              if key is a boolean array, then features at index with True are deleted
        '''
        if type(key)==int:
            del(self.SeqRecord.features[key])
        elif type(key)==range:
            self.delete_sequence(key, inplace=True)
        elif type(key)==slice:
            start = key.start
            stop = key.stop
            if start==None:
                start = 0
            if stop==None:
                stop = len(out.SeqRecord.seq)
            self.delete_sequence(range(start, stop), inplace=True)
        elif isBool(key):
            keep = np.arange(len(key))[key==False]
            self.SeqRecord.features = [self.SeqRecord.features[i] for i in keep]
        else:
            sys.error()

    def translate(self, strand=None, frame=0, table='Standard'):
        '''
        Translate the DNA sequence to amino acids
        frame = codon reading frame such as 0, 1, or 2
        table = codon table to use. Standard codon table is used by default
        return a string
        '''
        if len(self.SeqRecord.features) == 0:
            strand = 1 
        elif strand==None:
            strand = self.SeqRecord.features[0].location.strand
        
        if strand==1:
            AA_seq = self.SeqRecord.seq[frame:].translate(table=table)
        elif strand==-1:
            AA_seq = self.SeqRecord.seq.reverse_complement()[frame:].translate(table=table)
        return str(AA_seq)

    def to_dataframe(self, str_only=True):
        '''
        Obtain a dataframe of genetic features on the Plasmid
        str_only = export objects formatted into strings 
        '''
        data = []
        for gene in self.SeqRecord.features:
            label = gene.qualifiers['locus_tag'][0]
            feature = gene.type
            # figure out strand location
            [start, end, strand] = unwrap_FeatureLocation(gene.location, len(self.SeqRecord.seq))
            if strand > 0:
                color = gene.qualifiers['ApEinfo_fwdcolor'][0]
            else:
                color = gene.qualifiers['ApEinfo_revcolor'][0]
            ''' debug: will break auto annotations on Plasmid merge
            # handle sequences which wrap around the origin
            if start > end:
                seq = str(self.SeqRecord.seq[start:]) + str(self.SeqRecord.seq[:end])
            else: 
                seq = str(self.SeqRecord.seq[start:end])
            # reorient the sequence if needed
            if fwd_only and strand==-1:
                seq = reverse_complement(seq)
                strand = 1
            '''
            seq = gene.extract(self.SeqRecord.seq)
            length = len(seq)
            loc = str(gene.location)
            if str_only:
                data.append([label, feature, start, end, length, strand, str(gene.location), color])
            else:
                data.append([label, feature, start, end, length, strand, gene.location, color])
        col = ['locus_tag','type','start','end','length','strand','location','color']
        return pd.DataFrame(data, columns=col)

    def __repr__(self):
        '''
        Return string showing memory address and table contents
        '''
        # get memory address
        loc = str(type(self))+' at '+hex(id(self))
        return str(loc)+'\n' + self.get_table()

    def __str__(self):
        '''
        Returns nucleotide sequence of Plasmid dataframe
        '''
        return str(self.SeqRecord.seq)

    def __len__(self):
        '''
        Return the number of genetic features
        '''
        return len(self.SeqRecord.features)
    
    def get_colored(self, width=0):
        '''
        Print out colored representation of Plasmid
        width = number of characters per line
        '''
        rec = Graphic(data=self)
        out = rec.colorize_plasmid(width=width)
        return out
 
    def get_table(self):
        # add genbank header info
        table = ''
        for k in self.SeqRecord.annotations.keys():
            table+=str(k+':'+str(self.SeqRecord.annotations[k]))+'\n' 
        df = self.to_dataframe()
        table+= str(df[['locus_tag','type','location','length','color']])+'\n'
        # colorize the colors
        for c in np.unique(df['color']):
            val = Graphic().get_colored_text(c,fg=c)
            table = table.replace(c,val)
        table+= 'total length:'+str(len(self.__str__()))+'\n'
        return table

    def annotate(self, name, sequence, feature='unknown', color=None, circular=True, inplace=False):
        '''
        Adds annotations to a plasmid using a parts library
        name = name of the gene
        sequence = DNA or amino acid sequence of the feature
        feature = type of genetic feature such as cds, mRNA, primer_bind
        color = [fwd_color, rev_color] to use
        inplace = performs modifications inplace
        returns a modified plasmid dataframe
        '''
        # remove gaps
        if type(sequence) == str:
            sequence = sequence.replace(' ','')
        if sequence=='':
            return self

        # apply default colors if it is not provided
        if color==None:
            color = search_key(self.params['colors']['locus_tag'], name, rev=True)
        if color==None:
            color = search_key(self.params['colors']['features'], feature, rev=True)

        # apply default colors
        if color==None:
            color = self.params['colors']['features']['unknown']
        
        # format to list
        if type(color) == str:
            color = [color, color]
        color = [color[0], color[-1]]
        
        out = self.inplace(inplace)
        # search for DNA sequence
        if isDNA(sequence):
            # convert U to T for DNA sequence only
            sequence = sequence.lower().replace('u','t')
            pos = Aligner.search_DNA(sequence, str(out.SeqRecord.seq), circular=circular)
        # search for protein sequence
        else:
            pos = Aligner.search_protein(sequence, str(out.SeqRecord.seq), circular=circular)
        
        new_features=[]
        # add the annotation
        for [start, end, strand] in pos[['start','end','strand']].values:
            # create compound location which can wrap around the origin
            loc = wrap_FeatureLocation(start, end, strand, len(out.SeqRecord.seq))
            qualifiers = {'locus_tag':[name], 'ApEinfo_fwdcolor':[color[0]], 'ApEinfo_revcolor':[color[1]]}
            gene = Bio.SeqFeature.SeqFeature(type=feature, location=loc, qualifiers=qualifiers)
            new_features.append(gene)

        # append to Plasmid and update
        out.SeqRecord.features+= new_features
        out.update(inplace=True)
        return out

    def annotate_by_reference(self, file, inplace=False):
        '''
        Annotate using a reference library of sequences
        '''
        out = self.inplace(inplace)
        ref = pd.read_csv(file)
        # check if color in the data column
        if ('color' in ref.columns)==False:
            ref['color'] = None
            
        cols = ['label','sequence','feature_type','color']
        if (cols in ref.columns)==False:
            cols = ['name','sequence','feature_type','color']
        for label, seq, ftype, color in ref[cols].values:
            out.annotate(name=label, sequence=seq, feature=ftype, color=color, inplace=True)
        return out

    def merge(self, value, how='left', inplace=False):
        '''
        Merge Plasmid dataframes together
        value = Plasmid or dataframe of features to merge into self
        how = if 'left' then features from value are appended to self
              if 'right' then features from self merged into value
              if 'and' then only intersection of features are returned
              if 'nand' then only non-intersecting features of self are returned
        returns Plasmid dataframe
        '''
        out = self.inplace(inplace)   
        if how=='left':
            if type(value).__name__=='Plasmid':
                value = value.to_dataframe()
            for v in value[['locus_tag','type','sequence','color']].values:
                out = out.annotate(name=v[0], feature=v[1], sequence=v[2], color=v[3])
            return out
        elif how=='right':
            return value.merge(self, how='left', inplace=inplace)
        elif how=='and': # debug todo refinement
            out.SeqRecord.features = []
            if type(value).__name__=='Plasmid':
                value = value.to_dataframe()
            for v in value[['locus_tag','type','sequence','color']].values:
                out = out.annotate(name=v[0], feature=v[1], sequence=v[2], color=v[3])
            return out 
        elif how=='nand': # needs refinement debug
            keep = [True]*len(out.SeqRecord.features)
            df = out.to_dataframe()
            for q in value['sequence'].values:
                for i in range(len(df)):
                    if keep[i] and q.lower()==df['sequence'].values[i].lower():
                        keep[i] = False
            x = np.arange(len(keep))
            out.SeqRecord.features = [out.SeqRecord.features[i] for i in x[keep]]
            return out
        else:
            sys.error()

    def drop_duplicates(self, inplace=False): 
        '''
        Shift through feature list and remove duplicated annotations
        '''
        out = self.inplace(inplace)
        return out[out.duplicates()==False]

    def duplicates(self): 
        '''
        Sort through feature list and return True/False index of entries which maybe duplicates
        '''
        df = self.update().to_dataframe()
        dupe = [False]*len(df)
        for i in range(len(df)-1):
            if dupe[i]==False:
                c1 = df.iloc[i][['start','end']].values
                for j in range(i+1,len(df)):
                    c2 = df.iloc[j][['start','end']].values
                    if c1[0]==c2[0] and c1[1]==c2[1]:
                        dupe[j] = True
                    if c2[0] > c1[0]:
                        break
        return np.array(dupe)
 
    def __radd__(self, other):
        '''
        Performs Plasmid operation other + self
        Other can be a str
        Self is a Plasmid dataframe
        returns a modified Plasmid
        '''
        # make Plasmids linear by splitting at the origin
        seq1 = self.linearize()
        if type(other) == str:
            seq2 = Plasmid(other)
        else:
            seq2 = other.linearize()
        
        # concatenate seq records
        sr = seq2.SeqRecord + seq1.SeqRecord
        seq1.SeqRecord = seq2.SeqRecord + seq1.SeqRecord
        seq1.update(inplace=True)
        return seq1
 
    def __add__(self, other):
        '''
        Performs Plasmid operation self + other
        Other can be a str
        Self is a Plasmid dataframe
        return a modified Plasmid
        '''
        # make Plasmids linear by splitting at the origin
        seq1 = self.linearize()
        if type(other) == str:
            seq2 = Plasmid(other)
        else:
            seq2 = other.linearize()
        
        # concatenate seq records
        seq1.SeqRecord = seq1.SeqRecord + seq2.SeqRecord
        seq1.update(inplace=True)
        return seq1
    
    def linearize(self, inplace=False):
        '''
        Edit SeqRecord so Plasmid is now linear
        returns a linearized Plasmid
        '''
        if inplace:
            out = self
        else:
            out = self.copy()
        N = len(out.SeqRecord.seq)
        new = []
        keep = []
        for gene in out.SeqRecord.features:
            [start, end, strand] = unwrap_FeatureLocation(gene.location, len(out.SeqRecord.seq))
            # check if gene overlaps origin
            if start > end:
                # partition the component locations
                loc1 = []
                loc2 = []
                for L in gene.location.__dict__['parts']:
                    if L.end <= end:
                        loc1.append(L)
                    else:
                        loc2.append(L)
                gene1 = copy.deepcopy(gene)
                gene2 = copy.deepcopy(gene)
                gene1.location = Bio.SeqFeature.CompoundLocation(loc1)
                gene2.location = Bio.SeqFeature.CompoundLocation(loc2)
                # add some information to label tag to show it was split
                gene1.qualifiers['locus_tag'][0]+=' subsequence'
                gene2.qualifiers['locus_tag'][0]+=' subsequence'
                new.append(gene1)
                new.append(gene2)
            else:
                keep.append(gene)
        out.SeqRecord.features = new + keep
        return out

    def shift_features(self, N, circular=False, inplace=False):
        '''
        Shifts base index of sequences by N base pairs
        record = SeqRecord to modify
        N = number of base pairs to shift
        circular = If SeqRecord is a Plasmid, which is circular, shift the index circularly around origin
        inplace = If false, a copy of modified self is returned and original data is not modified
        returns a Plasmid dataframe
        '''
        out = self.inplace(inplace)
        seq_length = len(out.SeqRecord.seq)
        for gene in out.SeqRecord.features:
            gene.location = shift_CompoundLocation(gene.location, N, seq_length, circular=circular)
        return out

