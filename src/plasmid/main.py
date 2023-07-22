#!/usr/bin/env python
# Main function to execute commandline utilities

# local libraries
from .aligner import Aligner
from .oligos import Oligos

# system and time
import argparse

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('-v', dest='verbose', action='store_true', help='print verbose info')
    parser.add_argument('-c', dest='config', action='store_true', help='print verbose info')
    subparser = parser.add_subparsers(title='subcommands', dest='subcommand')
    
    # parse commands related to sequence alignment
    aparser = subparser.add_parser('aligner', help='suboptions for sequence alignment')
    aparser.add_argument('-q', dest='query', nargs='+', type=str, help='query sequences as fasta, csv, genbank')
    aparser.add_argument('-d', dest='reference', nargs='+', type=str, help='reference sequences as fasta, csv, genbank')
    aparser.add_argument('-o', dest='ofile', type=str, help='output filename')
    aparser.add_argument('-sort', dest='sort', nargs=2, default=['AS','idxmax'], type=str, help='output filename')
    aparser.add_argument('-m', dest='method', type=str, default='parasail', help='method to use')
    
    # commands to generate oligos 
    oligo_parser = subparser.add('oligos', help='suboptions for oligo generation')
    oligo_parser.add_argument('-i', dest='infile', nargs='+', type=str, help='input csv files with strands to process')
    oligo_parser.add_argument('-o', dest='outfile', type=str, help='output filename for primers or oligos')
    oligo_parser.add_argument('-t', dest='strands', nargs='+', type=str, help='header of strand name to process')
    oligo_parser.add_argument('-m',dest='method', type=str, default='twist_outer', help='oligo synthesis method')
    oligo_parser.add_argument('-g', dest='ggsite', nargs='+', type=str, help='oligo_parser overhang')
    oligo_parser.add_argument('-pad', dest='padding', type=int, default=300, help='pad ends with random base pairs so total length is at least N base pairs')
    oligo_parser.add_argument('-offset', dest='offset', type=int, nargs='+', default=[], help='offset strands in build order by N')
    oligo_parser.add_argument('-noterm', dest='noterm' , action='store_true', help='strip terminator sequences')
    oligo_parser.add_argument('-mask', dest='mask' , type=str, nargs='+', help='strip custom sequences from the oligos')
    oligo_parser.add_argument('-nonum', dest='renumber', action='store_false', help='do not renumber oligos')

    # commands to generate genbank files
    cparser = subparser.add('genbank', help='suboptions for genbank file generation')
    cparser.add_argument('-i', dest='infile', nargs='+', type=str, help='input csv files with strands to process')
    cparser.add_argument('-o', dest='outfile', type=str, help='output filename for primers or oligos')
    cparser.add_argument('-t', dest='strands', nargs='+', type=str, help='header of strand name to process')
    cparser.add_argument('-p', dest='pfile', nargs='+', type=str, default=[], help='plasmids to insert sequences into')
    cparser.add_argument('-l', dest='loc', type=int, nargs=2, default=[None,None], help='plasmid bp indices to replace with insert')
    cparser.add_argument('-ori', dest='origin', type=int, default=0, help='set origin to new bp index')
    cparser.add_argument('-rs', dest='reset_colors', action='store_true', help='reset plasmid colors')

    # parse arguments
    args = parser.parse_args()
    if args.subcommand=='aligner':
        opt = Aligner(args)
        opt.run(args.qry, args.ref)
    elif args.subcommand=='oligos' or args.subcommand=='genbank':
        opt = Oligos(args)
        opt.run()
