#!/usr/bin/env python

# computing libraries
import pandas as pd

# system libraries
import os
import subprocess
import tempfile

from . import io


# columns requested from mmseqs convertalis
FORMAT_OUTPUT = ('query,qstart,qend,qlen,target,tstart,tend,tlen,'
                 'qframe,raw,evalue,fident,alnlen,mismatch')

RESULT_COLUMNS = ['q_id', 'q_start', 'q_end', 'q_len',
                  't_id', 't_start', 't_end', 't_len',
                  'qframe', 'AS', 'evalue', 'match_score',
                  'tot', 'NM']


def read_MMSEQ_results(fname: str) -> pd.DataFrame:
    '''
    Parse an mmseqs2 convertalis results file into a pandas dataframe.

    Args:
        fname: filename of the mmseqs2 tabular results file.

    Returns:
        A pandas dataframe of alignment results.
    '''
    columns = RESULT_COLUMNS + ['match', 'orientation']
    if not os.path.exists(fname) or os.path.getsize(fname) == 0:
        df = pd.DataFrame()
    else:
        df = pd.read_csv(fname, delimiter='\t', header=None)

    # return an empty, typed dataframe for empty results
    if len(df) == 0:
        data = {}
        int_cols = ['q_start', 'q_end', 'q_len', 't_start', 't_end', 't_len',
                    'qframe', 'AS', 'tot', 'NM', 'match']
        float_cols = ['evalue', 'match_score']
        for c in int_cols:
            data[c] = pd.Series(dtype='int64')
        for c in float_cols:
            data[c] = pd.Series(dtype='float64')
        for c in ['q_id', 't_id', 'orientation']:
            data[c] = pd.Series(dtype='object')
        return pd.DataFrame(data)[columns]

    df.columns = RESULT_COLUMNS

    # format numeric columns
    int_cols = ['q_start', 'q_end', 'q_len', 't_start', 't_end', 't_len',
                'qframe', 'AS', 'tot', 'NM']
    for c in int_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0).astype(int)
    for c in ['evalue', 'match_score']:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0.0).astype(float)

    # compute the number of matching residues
    df['match'] = df['tot'] - df['NM']

    # set the orientation and put query coordinates on the forward strand
    c = df['q_start'] > df['q_end']
    df['orientation'] = '+'
    df.loc[c, 'orientation'] = '-'
    q_start = df.loc[c, 'q_start'].copy()
    df.loc[c, 'q_start'] = df.loc[c, 'q_end']
    df.loc[c, 'q_end'] = q_start

    # make coordinates zero indexed
    for col in ['q_start', 'q_end', 't_start', 't_end']:
        df[col] = df[col] - 1

    df = df.sort_values(by=['q_id', 'q_start']).reset_index(drop=True)
    return df


class Aligner:
    '''
    This class holds functions pertaining to sequence clustering and
    database search with mmseqs2.

    The class takes fasta and fastq filenames as input or pandas dataframes
    with columns [name, sequence]. Clustering and search results are parsed
    into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        binary: mmseqs2 binary.
        config_createdb: mmseqs2 subcommand for database creation.
        config_cluster: mmseqs2 subcommand for clustering.
        config_blastx: mmseqs2 subcommand for translated search.
        config_linsearch: mmseqs2 subcommand for linear search.
        config_createindex: mmseqs2 subcommand for building a linear index.
        config_convertalis: mmseqs2 subcommand for converting alignments.
        min_seq_id: minimum sequence identity passed to search workflows.
        sensitivity: search sensitivity passed to search workflows.
        coverage: alignment coverage passed to search workflows.
        cov_mode: coverage mode passed to search workflows.
        evalue: maximum e-value passed to search workflows.
    '''
    verbose = False
    binary = 'mmseqs'
    config_createdb = 'createdb'
    config_cluster = 'easy-cluster'
    config_blastx = 'search'
    config_linsearch = 'easy-linsearch'
    config_createindex = 'createlinindex'
    config_convertalis = 'convertalis'
    min_seq_id = None
    sensitivity = None
    coverage = None
    cov_mode = None
    evalue = None

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='mmseqs2_')

    # ------------------------------------------------------------------
    # mmseqs2 database creation, clustering, and search
    # ------------------------------------------------------------------
    def createdb(self, dbfile: str, db_folder: str=None) -> str:
        '''
        Create an mmseqs2 database.

        Args:
            dbfile: fasta file of sequences to build a database from.
            db_folder: prefix of the output database. Defaults to a database
                in the working directory.

        Returns:
            A string containing the database prefix.
        '''
        if db_folder is None:
            db_folder = os.path.join(self.work_dir.name, 'database')
        cmd = [self.binary, self.config_createdb, str(dbfile), str(db_folder)]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 createdb failed. '+str(ret))
        return db_folder

    def createindex(self, db_folder: str, tmp_dir: str=None) -> str:
        '''
        Build the linear index used by repeated linsearch runs.

        Args:
            db_folder: fasta file or mmseqs2 database prefix of the target.
            tmp_dir: temporary directory for mmseqs2.

        Returns:
            A string containing the indexed database prefix.

        Raises:
            ValueError: if no database is given.
        '''
        if db_folder is None:
            raise ValueError('target database is required')
        if not os.path.exists(str(db_folder)+'.dbtype'):
            db_folder = self.createdb(db_folder)
        if tmp_dir is None:
            tmp_dir = os.path.join(self.work_dir.name, 'tmp')
        cmd = [self.binary, self.config_createindex, str(db_folder), tmp_dir]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 createindex failed. '+str(ret))
        return db_folder

    def cluster(self, qfile) -> pd.core.frame.DataFrame:
        '''
        Wrapper for mmseqs2 clustering.

        Args:
            qfile: fasta or fastq file of sequences to cluster or a dataframe with columns [name, sequence].

        Returns:
            A pandas dataframe of sequence clusters with columns [cluster_id, q_id, count].
        '''
        # write query to temporary file if a dataframe is provided
        if not isinstance(qfile, str):
            fname = os.path.join(self.work_dir.name, 'query.fa')
            io.write_fasta(fname, qfile)
            qfile = fname
        out_dir = os.path.join(self.work_dir.name, 'results')

        options = []
        if self.min_seq_id is not None:
            options += ['--min-seq-id', str(self.min_seq_id)]
        if self.sensitivity is not None:
            options += ['-s', str(self.sensitivity)]
        if self.coverage is not None:
            options += ['-c', str(self.coverage)]
        if self.cov_mode is not None:
            options += ['--cov-mode', str(self.cov_mode)]
        if self.evalue is not None:
            options += ['-e', str(self.evalue)]

        cmd = [self.binary, self.config_cluster, qfile, out_dir,
               self.work_dir.name] + options
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 cluster failed. '+str(ret))

        # load results
        ofile = out_dir + '_cluster.tsv'
        if not os.path.exists(ofile) or os.path.getsize(ofile) == 0:
            return pd.DataFrame(columns=['cluster_id', 'q_id', 'count'])
        df = pd.read_csv(ofile, header=None, delimiter='\t')
        df.columns = ['cluster_id', 'q_id']
        df = df.drop_duplicates()

        # get cluster statistics
        x = df.groupby(['cluster_id']).agg({'q_id': 'count'})
        x.columns = ['count']
        x = x.reset_index()

        # merge info
        df = df.merge(x, on='cluster_id', how='left')
        df = df.sort_values(by=['count', 'cluster_id']).reset_index(drop=True)
        return df

    def blastx(self, qfile: str, db_folder: str, out_file: str=None) -> pd.DataFrame:
        '''
        Run a translated nucleotide-to-protein search with mmseqs2.

        Nucleotide queries are translated in all six frames and searched
        against a protein database. If db_folder points to a fasta file, an
        mmseqs2 database is created.

        Args:
            qfile: fasta or fastq file of nucleotide reads.
            db_folder: fasta file or mmseqs2 database prefix of the protein
                reference.
            out_file: filename of the mmseqs2 results file.

        Returns:
            A pandas dataframe of read alignments.

        Raises:
            ValueError: if no database is given.
        '''
        work_dir = self.work_dir.name
        if db_folder is None:
            raise ValueError('database is required')
        if not os.path.exists(str(db_folder)+'.dbtype'):
            db_folder = self.createdb(db_folder)
        ref_db = str(db_folder)
        if out_file is None:
            out_file = os.path.join(work_dir, 'results.m8')

        # create the query database
        qry_db = os.path.join(work_dir, 'query_db')
        self.createdb(qfile, qry_db)

        options = []
        if self.min_seq_id is not None:
            options += ['--min-seq-id', str(self.min_seq_id)]
        if self.sensitivity is not None:
            options += ['-s', str(self.sensitivity)]
        if self.coverage is not None:
            options += ['-c', str(self.coverage)]
        if self.cov_mode is not None:
            options += ['--cov-mode', str(self.cov_mode)]
        if self.evalue is not None:
            options += ['-e', str(self.evalue)]

        # run the translated search
        aln_db = os.path.join(work_dir, 'aln_db')
        tmp_dir = os.path.join(work_dir, 'tmp')
        cmd = [self.binary, self.config_blastx, qry_db, ref_db,
               aln_db, tmp_dir, '-a'] + options
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 blastx failed. '+str(ret))

        # convert aln_db to a human readable format
        cmd = [self.binary, self.config_convertalis, '--format-output',
               FORMAT_OUTPUT, qry_db, ref_db, aln_db, out_file]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 convertalis failed. '+str(ret))

        return read_MMSEQ_results(out_file)

    def linsearch(self, qfile: str, db_folder: str, out_file: str=None) -> pd.DataFrame:
        '''
        Run a linear time search with mmseqs2.

        This is the mmseqs2 linear search workflow (easy-linsearch). It builds
        a linear index for the target and is faster than the standard search
        when there are few queries. A target database is generated from a
        fasta file automatically; a pre-built database can also be supplied.

        Args:
            qfile: fasta or fastq file of query sequences.
            db_folder: fasta file or mmseqs2 database prefix of the target.
            out_file: filename of the mmseqs2 results file.

        Returns:
            A pandas dataframe of alignment results.

        Raises:
            ValueError: if no database is given.
        '''
        work_dir = self.work_dir.name
        if db_folder is None:
            raise ValueError('target database is required')
        if not os.path.exists(str(db_folder)+'.dbtype'):
            db_folder = self.createdb(db_folder)
        if out_file is None:
            out_file = os.path.join(work_dir, 'results.m8')

        options = []
        if self.min_seq_id is not None:
            options += ['--min-seq-id', str(self.min_seq_id)]
        if self.sensitivity is not None:
            options += ['-s', str(self.sensitivity)]
        if self.coverage is not None:
            options += ['-c', str(self.coverage)]
        if self.cov_mode is not None:
            options += ['--cov-mode', str(self.cov_mode)]
        if self.evalue is not None:
            options += ['-e', str(self.evalue)]

        cmd = [self.binary, self.config_linsearch, qfile, str(db_folder),
               out_file, os.path.join(work_dir, 'tmp'),
               '--format-output', FORMAT_OUTPUT] + options
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 linsearch failed. '+str(ret))

        return read_MMSEQ_results(out_file)
