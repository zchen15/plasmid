#!/usr/bin/env python
# Class function for aligners

# computing libraries
import numpy as np
import pandas as pd
import scipy as sp

# reading and writing genbank files
import sklearn
import sklearn.cluster

# system and time
import re
import tempfile
import subprocess

# import custom libraries
from .misc import read_to_df
from .misc import revcomp
from .misc import load_args
from .fileIO import fileIO
from .aligner import Aligner

class Clust(Aligner):
    '''
    This class holds functions pertaining to sequence clustering
    '''
    def __init__(self, args=None):
        if args!=None:
            self.params = load_args(args, self.params)
        if self.params['verbose']:
            print(self.params)
        self.create_tempdir()
        
        self.params['tsne'] = {'metric':'euclidean',
                          'perplexity':30,
                          'lr':200,
                          'n_iter':1000}

        self.params['optics'] = {'min_samples':None,
                            'min_cluster_size':None,
                            'xi':0.05,
                            'max_eps':None,
                            'cluster_method':'dbscan',
                            'n_jobs':None,
                            'alt_label':False}
    
    def create_tempdir(self):
        '''
        Create a temporary working directory for aligners
        '''
        self.params['work_dir'] = tempfile.TemporaryDirectory(prefix='clust_')

    def pairs_to_adjacency(x, sym_larger=False):
        '''
        Converts a paired list into an adjacency matrix
        x = numpy array with columns [query, database, values]
        sym_larger = for symmetric matrices, favor the larger value
        return dictionary with entries 'matrix' and 'columns'
        '''
        cols = [i for i in x[:,0]] + [i for i in x[:,1]]
        cols = np.unique(cols)
        
        key = {cols[i]:i for i in range(len(cols))}
        dist = np.zeros([len(cols), len(cols)])
        
        for qry, db, val in x:
            i = key[qry]
            j = key[db]
            dist[i,j] = val
        
        for i in range(len(cols)):
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
        out = {'matrix':dist,
               'columns':cols}
        return out
        
    def get_closest_points(self, df, N):
        '''
        Get the N set of closest points
        '''
        
        
    def interpolate_distance(self, df, N):
        '''
        interpolate distances using nearest neighboring points
        '''
        
        
    def get_distance_matrix(self, query, database):
        '''
        Compute the distance matrix for a set of sequences
        query = dataframe with columns [name, sequence]
        database = dataframe with columns [name, sequence]
        return dictionary with 
                'matrix' as adjacency matrix of sequence similarity
                'columns' names of sequences representing each column
        '''
        df = self.minimap2(query=query, database=database)
        # filter to best alignments
        col = ['query_id','database_id']
        df = Aligner.filter_idx(df, col, 'match_score', 'idxmax')
        
        print('Converting paired list to adjacency matrix')
        metric = 'match_score'
        adj = Clust.pairs_to_adjacency(df[['query_id', 'database_id', metric]].values)
        mat = 1 - adj['matrix']
        # fill the diagonal values
        L = len(adj['columns'])
        for i in range(L):
            mat[i,i] = 0
        adj['matrix'] = mat
        return adj
        
    def pca(x, n_comp=2, recenter=False):
        '''
        Runs principle component analysis to compression dimensions 
        x = adjacency matrix
        n_comps = components to reduce to
        '''

        # initialize PCA settings
        print('running pca with n_comp = '+str(n_comp))

        # center the feature vectors
        if recenter:
            vT = x.T
            for i in range(0,len(vT)):
                vT[i] = vT[i] - np.mean(vT[i])
            x = vT.T

        # do PCA
        pca = sklearn.decomposition.SparsePCA(n_components = n_comp)
        v = pca.fit_transform(x)
        return v

    def tsne(df, n_comp=2):
        '''
        Run tsne on input data
        x = adjacency matrix
        '''
        print('running sklearn tsne with n_comp = '+str(n_comp))
        metric = self.params['tsne']['metric']
        perplexity = self.params['tsne']['perplexity']
        lr = self.params['tsne']['lr']
        n_iter = self.params['tsne']['n_iter']
        verbose = self.params['verbose']
        
        tsne = sklearn.manifold.TSNE(n_components=n_comp, perplexity=perplexity, metric=metric,
                             early_exaggeration=12.0, learning_rate=lr,
                             n_iter=n_iter, n_iter_without_progress=300, verbose=verbose)
        # get the vector
        v = tsne.fit_transform(x)
        return v

    def optics_sequence(self, qry):
        '''
        Find cluster centers for a set of sequences using OPTICS
        qry = dataframe of query sequences to search for cluster centers
               Must have at least the following columns:
               [name, sequence]
        csize = number of sequences to take from cluster center for multi-seq alignment
        '''
        # run the aligner
        print('cluster_compute: computing pairwise distance matrix')
        config = '-D --dual=no --for-only'
        self.params['minimap']['config'] = config
        x = self.get_distance_matrix(qry, qry)
        mat = x['matrix']

        # do optics on raw data
        clst = self.optics(mat, metric='precomputed')
        df = pd.DataFrame(x['columns'], columns=['name'])
        cols = ['labels','ordering','reachability']
        for c in cols:
            df[c] = clst[c]
        df = df.merge(qry, on='name', how='left')
        return df
    
    def optics(self, x, metric='precomputed'):
        '''
        Function to perform optics clustering on sequences
        x = adjacency matrix
        '''
        min_samples = self.params['optics']['min_samples']
        xi = self.params['optics']['xi']
        max_eps = self.params['optics']['max_eps']
        cluster_method = self.params['optics']['cluster_method']
        min_cluster_size = self.params['optics']['min_cluster_size']
        n_jobs = self.params['optics']['n_jobs']
        alt_label = self.params['optics']['alt_label']
        
        # running the clustering
        print('Running clust.OPTICS')

        # Safety check on number of samples given
        if len(x) < 3:
            ValueError('Sample size too small to use with OPTICS')

        # compute our own eps
        if max_eps == None:
            z = np.abs(x)
            max_eps = (np.max(z)+np.min(z))/2
        print('max_eps = '+str(max_eps))

        # do optimization on min_samples to get clustering with the least outliers
        if min_samples == None or min_samples > len(x):
            min_samples = np.ceil(len(x)/2)
            max_iter = 100
        else: # if min_samples is provided, dont optimize
            max_iter = 1
        
        # initialize settings
        prev_nout = len(x)
        prev_nclust = 0
        prev_min = min_samples*2
        best_clust = None
        best_nclust = 0
        best_min = prev_min
        best_nout = prev_nout
        delta = min_samples/2
        for i in range(0, max_iter):
            min_samples = int(min_samples)
            print('clust.OPTICS: iter='+str(i)+' using min_samples='+str(min_samples))
            clust = sklearn.cluster.OPTICS(min_samples=min_samples, min_cluster_size=min_cluster_size, xi=xi,
                                max_eps=max_eps, cluster_method=cluster_method, metric=metric, n_jobs=n_jobs)
            clust.fit(x)
            # check labels for cluster number and coverage
            nout = np.sum(clust.labels_ == -1)
            nclust = np.max(clust.labels_)+1
            print('clust.OPTICS: clusters='+str(nclust)+' outliers='+str(nout)+' delta='+str(delta))
            # always keep the change if we get more clusters
            if nclust > prev_nclust or (nclust == prev_nclust and nout <= prev_nout):
                tmp = min_samples
                delta = int((prev_min - min_samples)/2)
                min_samples-= delta
            # if cluster number goes down as we increase min_samples
            elif nclust < prev_nclust:
                tmp = min_samples
                delta = int((prev_min + min_samples)/2) - min_samples
                min_samples = prev_min + delta
            # record the best clustering
            if nclust > best_nclust or (nclust == best_nclust and nout <= best_nout):
                if np.min(clust.reachability_) >= 0:
                    best_clust = clust
                    best_min = min_samples
                    best_nclust = nclust
                    best_nout = nout
            # if shift in min_samples is zero, exit loop
            if delta==0 or min_samples <= 2 or nout==0:
                break
            prev_nclust = nclust
            prev_nout = nout
            prev_min = tmp

        # get the rest of the values
        s = best_clust.ordering_
        reach = best_clust.reachability_[s]
        labels = best_clust.labels_[s]
        r_id = np.arange(len(s))[s]

        # cap the reachability so fft can be computed
        y = np.copy(reach)
        y[reach >= np.inf] = -1
        reach[y == -1] = np.max(y)

        if np.min(reach) >= 0:
            # using custom relabeling on the reachability plot
            if alt_label and len(reach) > 5:
                print('using alt labeling')
                labels = Clust.reachability_alt_label(reach)
            print('n_clusters='+str(np.max(labels)+1)+' n_unclustered='+str(np.sum(labels==-1))+' N='+str(len(x)))
        else:
            print('clust.OPTICS: reachability < 0, please change the min_samples parameter you are using')

        # format and return the data        
        s = np.argsort(r_id)
        labels = labels[s]
        ordering = np.arange(len(s))[s]
        reach = reach[s]
        
        out = {'ordering':ordering,
               'reachability':reach,
               'labels':labels}
        return out

    def reachability_alt_label(reach):
        '''
        Custom implementation to mark cluster boundaries for optics
        '''
        # get the slope change
        d1 = sp.signal.fftconvolve(reach, [1,-1], mode='same')
        d1[:-1] = d1[1:]
        labels = np.zeros(len(reach))
        # flag steep drops as outliers
        km = sklearn.cluster.AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
        km.fit(np.transpose([d1]))
        c = []
        for j in np.unique(km.labels_):
            x = d1[km.labels_==j]
            c.append(np.mean(x))
        s = np.argsort(c)
        labels[(km.labels_==s[0])] = -1
        # flag high reachability as outliers
        y = np.copy(reach)
        y[(km.labels_==s[0])] = np.max(y)
        km = sklearn.cluster.AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
        km.fit(np.transpose([y]))
        c = []
        for j in np.unique(km.labels_):
            x = y[km.labels_==j]
            c.append(np.mean(x))
        s = np.argsort(c)
        labels[(km.labels_==s[-1])] = -1
        # assign new labels
        th = labels!=-1
        c = th[1:]!=th[:-1]
        cuts = np.arange(1,len(c)+1)[c]
        k = 0
        for i in range(1,len(cuts),2):
            s1 = cuts[i-1]
            s2 = cuts[i]
            labels[s1:s2] = k
            k+=1
        return labels

