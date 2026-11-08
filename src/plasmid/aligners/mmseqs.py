#!/usr/bin/env python
"""Class functions for sequence clustering and database search with mmseqs2."""

# computing libraries
import pandas as pd

# system and time
import tempfile
import subprocess
import glob

# import custom libraries
from . import io

__all__ = ['Aligners']


class Aligners:
    '''
    This class holds functions pertaining to sequence clustering and
    database search with mmseqs2.

    Attributes:
        params (dict): stores clustering parameters and run settings for mmseqs2
        work_dir (tempfile.TemporaryDirectory): temporary working directory for intermediate files
    '''
    params = {'verbose': False}
    params['mmseqs2'] = {'cluster': 'easy-cluster',
                        'blastx': 'blastx',
                        'createdb': 'createdb',
                        'linsearch': 'easy-linsearch',
                        'taxonomy': 'easy-taxonomy'}
    params['work_dir'] = None
    params['binaries'] = {'mmseqs2': 'mmseqs'}

    def __init__(self, args=None):
        '''
        Initialize the Aligners object.

        Args:
            args: Command line arguments or dictionary of settings. Defaults to None.

        Raises:
            KeyError: if args contains a key that is not an Aligners setting.
        '''
        params = dict(self.params)
        if args is not None:
            values = args if isinstance(args, dict) else vars(args)
            for key, value in values.items():
                if value is None:
                    continue
                if key not in params:
                    raise KeyError('unknown Aligners parameter: '+str(key))
                params[key] = value
        self.params = params
        if self.params['verbose']:
            print(self.params)
        self.create_tempdir()

    def __del__(self):
        '''
        Destructor class to clean up temporary files
        '''
        try:
            self.params['work_dir'].cleanup()
        except Exception:
            pass

    def __enter__(self):
        '''
        Entry function for use in with statements
        '''
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        '''
        Exit function for use in with statements
        '''
        self.cleanup()
        return False

    def create_tempdir(self):
        '''
        Create a temporary working directory for the clusterer. Name of the temporary working directory is stored in self.params['work_dir']
        '''
        self.params['work_dir'] = tempfile.TemporaryDirectory(prefix='mmseqs2_')

    def cleanup(self):
        '''
        Remove the temporary working directory.
        '''
        if self.params['work_dir'] is not None:
            self.params['work_dir'].cleanup()
            self.params['work_dir'] = None

    def mmseq2_get_database(self, database: str, outdir: str) -> str:
        '''
        Wrapper for mmseqs2 database creation.

        Args:
            database: database to fetch or fasta file of sequences to build a database from.
            outdir: name of the output directory of the database.

        Returns:
            A string containing the output directory.
        '''
        mmseq = self.params['binaries']['mmseqs2']
        work_dir = self.params['work_dir']
        cmd = mmseq +' database '+str(database)+' '+str(outdir)+' '+work_dir.name
        subprocess.run(cmd.split(' '))
        return outdir

    def mmseq2_createdb(self, query: str, qry_db: str, config: str = 'createdb', mmbin: str = None) -> list:
        '''
        Create mmseq2 database files.

        Args:
            query: fasta file of amino acid reference sequences.
            qry_db: prefix for the mmseq2 query database.
            config: config options to pass to mmseq2.
            mmbin: location of the mmseq2 binary. Defaults to the mmseqs2 binary.

        Returns:
            A list of mmseq2 database files.
        '''
        if mmbin is None:
            mmbin = self.params['binaries']['mmseqs2']
        cmd = [mmbin, 'createdb', query, qry_db]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 createdb failed. '+str(ret))
        files = glob.glob(str(qry_db)+'*')
        return files

    def mmseq2_cluster(self, query, config: str = 'easy-cluster', mmbin: str = None) -> pd.core.frame.DataFrame:
        '''
        Wrapper for mmseq2 clustering.

        Args:
            query: fasta or fastq file of sequences to cluster or a dataframe with columns [name, sequence].
            config: config options to pass to mmseq2. Defaults to 'easy-cluster'.
            mmbin: location of the mmseq2 binary. Defaults to the mmseqs2 binary.

        Returns:
            A pandas dataframe of sequence clusters with columns [cluster_id, q_id, count].
        '''
        if mmbin is None:
            mmbin = self.params['binaries']['mmseqs2']
        work_dir = self.params['work_dir']

        # write query to temporary file if a dataframe is provided
        if not isinstance(query, str):
            qfile = work_dir.name + '/query'
            qfile = io.write_fastx(qfile, query)
        else:
            qfile = query

        out_dir = work_dir.name + '/results'
        cmd = [mmbin] + str(config).split(' ') + [qfile, out_dir, work_dir.name]
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

    def mmseq2_blastx(self, ref: str, query: str, config: str = 'blastx', mmbin: str = None) -> pd.core.frame.DataFrame:
        '''
        Run blastx with mmseqs.

        Args:
            ref: fasta file of amino acid reference sequences.
            query: fasta or fastq file of reads.
            config: config options to pass to mmseq2. Defaults to 'blastx'.
            mmbin: location of the mmseq2 binary. Defaults to the mmseqs2 binary.

        Returns:
            A pandas dataframe of read alignments.
        '''
        if mmbin is None:
            mmbin = self.params['binaries']['mmseqs2']
        work_dir = self.params['work_dir']

        # create mmseq database files
        qry_db = work_dir.name + '/nuc_db'
        ref_db = work_dir.name + '/AA_db'
        cmd1 = [mmbin, 'createdb', query, qry_db]
        cmd2 = [mmbin, 'createdb', ref, ref_db]
        for c in [cmd1, cmd2]:
            ret = subprocess.run(c)
            if ret.returncode != 0:
                raise RuntimeError('mmseqs2 createdb failed. '+str(ret))

        # run blastx search
        aln_db = work_dir.name + '/aln_db'
        tmp_dir = work_dir.name + '/tmp'
        cmd = [mmbin, 'search', qry_db, ref_db, aln_db, tmp_dir]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 search failed. '+str(ret))

        # convert aln_db to human readable format
        out_file = work_dir.name + '/results.csv'
        out_format = 'query,qstart,qend,qlen'
        out_format += ',target,tstart,tend,tlen'
        out_format += ',qframe,raw,evalue,fident,alnlen,mismatch'
        cmd = [mmbin, 'convertalis', '--format-output', out_format, qry_db, ref_db, aln_db, out_file]
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

    def mmseq2_search(self, query, config: str = 'easy-linsearch', mmbin: str = None) -> str:
        '''
        Wrapper for mmseq2 search.

        Args:
            query: fasta or fastq file of sequences to search.
            config: config options to pass to mmseq2. Defaults to 'easy-linsearch'.
            mmbin: location of the mmseq2 binary. Defaults to the mmseqs2 binary.

        Returns:
            A string containing the output directory of the search results.
        '''
        if mmbin is None:
            mmbin = self.params['binaries']['mmseqs2']
        work_dir = self.params['work_dir']

        # write query to temporary file if a dataframe is provided
        if not isinstance(query, str):
            qfile = work_dir.name + '/query'
            qfile = io.write_fastx(qfile, query)
        else:
            qfile = query

        out_dir = work_dir.name + '/results'
        cmd = [mmbin, config, qfile, out_dir, work_dir.name]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 '+config+' failed. '+str(ret))
        return out_dir

    def mmseq2_taxon(self, query, config: str = 'easy-taxonomy', mmbin: str = None) -> str:
        '''
        Wrapper for mmseq2 taxonomy.

        Args:
            query: database to assign taxonomy to.
            config: config options to pass to mmseq2. Defaults to 'easy-taxonomy'.
            mmbin: location of the mmseq2 binary. Defaults to the mmseqs2 binary.

        Returns:
            A string containing the output directory of the taxonomy results.
        '''
        if mmbin is None:
            mmbin = self.params['binaries']['mmseqs2']
        work_dir = self.params['work_dir']

        out_dir = work_dir.name + '/results'
        cmd = [mmbin, config, str(query), out_dir, work_dir.name]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            raise RuntimeError('mmseqs2 '+config+' failed. '+str(ret))
        return out_dir
