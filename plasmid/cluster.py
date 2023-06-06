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
from .aligner import Aligner

class Clust(Aligner):
    '''
    This class holds functions pertaining to sequence clustering
    '''
    params = {'verbose':False}
    params['work_dir'] = None

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
    
    def create_tempdir(self):
        '''
        Create a temporary working directory for aligners
        '''
        self.params['work_dir'] = tempfile.TemporaryDirectory(prefix='clust_')

    def get_symmetric_matrix(df, sym_larger=False):
        '''
        Makes a square distance matrix symmetric
        df = pandas dataframe with columns [id, d_i, ..., d_n]
             id columns should equal d_i, ..., d_n columns
        sym_larger = for symmetric matrices, favor the larger value
        '''
        cols = df.columns[df.columns!='id']
        dist = df[cols].values
        for i in range(0, len(dist)):
            v1 = dist[i,i:]
            v2 = dist[i:,i]
            c1 = (v1 == 0)
            c2 = (v2 == 0)
            c3 = (v2 > v1)
            if sym_larger:
                v1[c3] = v2[c3] # larger value is kept
                dist[i,i:] = v1
                dist[i:,i] = v1
            else: # favor smaller distance, 0 = invalid
                s1 = c1 & ~c2
                s2 = ~c2 & ~c3
                v1[s1] = v2[s1]
                v1[s2] = v2[s2]
                dist[i,i:] = v1
                dist[i:,i] = v1

        # make dist matrix into dataframe
        dist = pd.DataFrame(dist, columns = cols)
        dist['id'] = df['id'].values
        return dist
        
    def pca(x, n_comp=2, recenter=False):
        '''
        Runs principle component analysis to compression dimensions 
        df = pandas dataframe with columns [id, d_i, ..., d_n]
        n_comps = components to reduce to
        '''

        # initialize PCA settings
        logging.info('running pca with n_comp = '+str(n_comp))

        # get the vector
        col = df.columns[df.columns!='id']
        v = df[col].values

        # center the feature vectors
        if recenter:
            vT = v.T
            for i in range(0,len(vT)):
                vT[i] = vT[i] - np.mean(vT[i])
            v = vT.T

        # do PCA
        pca = sklearn.decomposition.SparsePCA(n_components = n_comp)
        v = pca.fit_transform(v)

        col = ['f_'+str(i) for i in range(0, n_comp)]
        data = pd.DataFrame(v, columns = col)
        data['id'] = df['id'].values
        return data

    def tsne(df, n_comp=2, metric='euclidean', perplexity=30, lr=200, n_iter=1000, verbose=0):
        '''
        Run tsne on input data
        df = pandas dataframe with columns [id, v_0, ..., v_i]
        '''
        logging.info('running sklearn tsne with n_comp = '+str(n_comp))
        tsne = sklearn.manifold.TSNE(n_components=n_comp, perplexity=perplexity, metric=metric,
                             early_exaggeration=12.0, learning_rate=lr,
                             n_iter=n_iter, n_iter_without_progress=300, verbose=verbose)
        # get the vector
        col = df.columns[df.columns!='id']
        v = tsne.fit_transform(df[col].values)
        cols = ['f_'+str(i) for i in range(0, n_comp)]
        data = pd.DataFrame(v, columns = cols)
        data['id'] = df['id'].values
        return data

    def optics(self, qry, csize=20, metric='similarity', outliers=False, pw_config ='-k15 -w10 -D', msa_config='-l 0 -r 0', workspace='./clust_run/'):
    '''
    Find cluster centers for a set of sequences using OPTICS
    qry = dataframe of query sequences to search for cluster centers
           Must have at least the following columns:
           [name, sequence]
    csize = number of sequences to take from cluster center for multi-seq alignment
    workspace = workspace folder for the aligner
    '''
    # run the aligner
    logging.info('cluster_compute: computing pairwise distance matrix')
    config = pw_config+' --dual=no --for-only'
    df = self.minimap2(query=qry, database=qry)

    # catch error in case pairwise alignment fails
    df = bpy.get_best(df, ['query_id','database_id'], 'AS')
    df = bpy.get_feature_vector(df[['query_id','database_id',metric]], symmetric=True)
    df = bpy.get_symmetric_matrix(df, sym_larger=False)
    
    # convert to distance
    logging.info('preparing precomputed data')
    for i in range(0,len(df)):
        df.iloc[:,i] = 1-df.iloc[:,i]
        df.iat[i,i] = 0
    
    # do optics on raw data
    logging.info('cluster_compute: running optics')
    df_clst = bpy.cluster_OPTICS(df, metric='precomputed', alt_label=False)
    df_clst = df_clst.merge(df_q, on='id', how='left')
    df_o = df_clst[df_clst['cluster_id'] == -1]
    df_o['outlier'] = True
    df_i = df_clst[df_clst['cluster_id'] > -1]
    
    # do merge via msa if needed
    if len(df_i) > 0:
        din = []
        for cid in np.unique(df_i['cluster_id']):
            df = df_i[df_i['cluster_id'] == cid]
            din.append(df.iloc[:csize])
        df_i = bpy.cluster_spoa_merge(pd.concat(din), config=msa_config, workspace=workspace, cleanup=True)
        df_i['outlier'] = False
    
    # return results
    if outliers:
        col = ['id','sequence','outlier']
        if len(df_i) > 0:
            return pd.concat([df_i,df_o])[col]
        else:
            return df_o[col]
    else:
        col = ['id','sequence']
        if len(df_i) == 0: 
            logging.warning('cluster_compute: no clusters found')
        return  df_i[col]
    