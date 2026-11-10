#!/usr/bin/env python

# computing libraries
import pandas as pd

# system libraries
import glob
import os
import subprocess
import tempfile

from . import io

class Aligner:
    '''
    This class holds functions pertaining to sequence clustering and
    database search with mmseqs2.

    The class takes fasta and fastq filenames as input or pandas dataframes
    with columns [name, sequence]. Clustering and search results are parsed
    into pandas dataframes.

    Attributes:
        verbose: print run settings when True.
        config: mmseqs2 subcommand for each operation.
        binary: mmseqs2 binary.
    '''
    verbose = False
    binary = 'mmseqs'
    config = {'database': 'database',
              'cluster': 'easy-cluster',
              'blastx': 'blastx',
              'createdb': 'createdb',
              'linsearch': 'easy-linsearch',
              'taxonomy': 'easy-taxonomy'}

    def __init__(self, **kwargs):
        '''
        Initialize the Aligner object.

        Args:
            kwargs: settings that override the class defaults.

        Raises:
            TypeError: if kwargs contains a name that is not an Aligner setting.
        '''
        for key, value in kwargs.items():
            if not hasattr(self, key):
                continue
            if value is not None:
                setattr(self, key, value)
        self.work_dir = tempfile.TemporaryDirectory(prefix='mmseqs2_')

    # ------------------------------------------------------------------
    # mmseqs2 database creation, clustering, and search
    # ------------------------------------------------------------------
    def get_database(self, dbfile: str, outdir: str) -> str:
        '''
        Wrapper for mmseqs2 database creation.

        Args:
            dbfile: fasta file of sequences to build a database from.
            outdir: name of the output directory of the database.

        Returns:
            A string containing the output directory.
        '''
        cmd = [self.binary, self.config['database'], str(dbfile), str(outdir), self.work_dir.name]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 database failed. '+str(ret))
        return outdir

    def createdb(self, qfile: str, qry_db: str) -> list:
        '''
        Create mmseq2 database files.

        Args:
            qfile: fasta file of amino acid reference sequences.
            qry_db: prefix for the mmseq2 query database.

        Returns:
            A list of mmseq2 database files.
        '''
        cmd = [self.binary, self.config['createdb'], qfile, qry_db]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 createdb failed. '+str(ret))
        files = glob.glob(str(qry_db)+'*')
        return files

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
            qfile = io.write_fasta(os.path.join(self.work_dir.name, 'query.fa'), qfile)
        out_dir = os.path.join(self.work_dir.name, 'results')
        cmd = [self.binary] + str(self.config['cluster']).split(' ') + [qfile, out_dir, self.work_dir.name]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 cluster failed. '+str(ret))

        # load results
        ofile = out_dir + '_cluster.tsv'
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

    def blastx(self, ref: str, qfile: str) -> pd.core.frame.DataFrame:
        '''
        Run blastx with mmseqs.

        Args:
            ref: fasta file of amino acid reference sequences.
            qfile: fasta or fastq file of reads.

        Returns:
            A pandas dataframe of read alignments.
        '''
        work_dir = self.work_dir

        # create mmseq database files
        qry_db = os.path.join(work_dir.name, 'nuc_db')
        ref_db = os.path.join(work_dir.name, 'AA_db')
        cmd1 = [self.binary, self.config['createdb'], qfile, qry_db]
        cmd2 = [self.binary, self.config['createdb'], ref, ref_db]
        for c in [cmd1, cmd2]:
            ret = subprocess.run(c)
            if ret.returncode != 0:
                raise RuntimeError('mmseqs2 createdb failed. '+str(ret))

        # run blastx search
        aln_db = os.path.join(work_dir.name, 'aln_db')
        tmp_dir = os.path.join(work_dir.name, 'tmp')
        cmd = [self.binary, 'search', qry_db, ref_db, aln_db, tmp_dir]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 search failed. '+str(ret))

        # convert aln_db to human readable format
        out_file = os.path.join(work_dir.name, 'results.csv')
        out_format = 'query,qstart,qend,qlen'
        out_format += ',target,tstart,tend,tlen'
        out_format += ',qframe,raw,evalue,fident,alnlen,mismatch'
        cmd = [self.binary, 'convertalis', '--format-output', out_format, qry_db, ref_db, aln_db, out_file]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 convertalis failed. '+str(ret))

        df = pd.read_csv(out_file, delimiter='\t', header=None)
        cols = ['q_id','q_start','q_end','q_len']
        cols += ['t_id','t_start','t_end','t_len']
        cols += ['qframe','AS','evalue','match_score','alnlen','mismatch']
        df.columns = cols
        df['match'] = df['alnlen'] - df['mismatch']

        # make things zero indexed
        cols = ['q_start','q_end','t_start','t_end']
        for col in cols:
            df[col] -= 1

        # set orientation correctly
        c = df['q_start'] > df['q_end']
        df['orientation'] = '+'
        df.loc[c, 'orientation'] = '-'
        q_start = df.loc[c, 'q_start']
        df.loc[c, 'q_start'] = df.loc[c, 'q_end']
        df.loc[c, 'q_end'] = q_start

        df = df.sort_values(by=['q_id','q_start'])
        return df

    def search(self, qfile, out_dir: str = None) -> str:
        '''
        Wrapper for mmseq2 search.

        Args:
            qfile: fasta or fastq file of sequences to search.

        Returns:
            A string containing the output directory of the search results.
        '''
        if out_dir is None:
            out_dir = os.path.join(self.work_dir.name, 'results')
        cmd = [self.binary, self.config['linsearch'], qfile, out_dir, self.work_dir.name]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 '+self.config['linsearch']+' failed. '+str(ret))
        return out_dir

    def taxon(self, qfile, out_dir: str = None) -> str:
        '''
        Wrapper for mmseqs2 taxonomy.

        Args:
            qfile: database to assign taxonomy to.

        Returns:
            A string containing the output directory of the taxonomy results.
        '''
        if out_dir is None:
            out_dir = os.path.join(self.work_dir.name, 'results')
        cmd = [self.binary, self.config['taxonomy'], str(qfile), out_dir, self.work_dir.name]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 '+self.config['taxonomy']+' failed. '+str(ret))
        return out_dir
