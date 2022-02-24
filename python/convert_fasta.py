#!/usr/bin/env python

# converts fasta to oneline fasta

from Bioinfo import oneline_fasta

import argparse

# def args
def parse_args():
    parser = argparse.ArgumentParser(description='converts fasta to oneline fasta')
    parser.add_argument('--input_file',help='input fasta file')
    parser.add_argument('--oneline_fasta_file',help='desired output file name')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    oneline_fasta(args.input_file,args.oneline_fasta_file)

if __name__ == '__main__':
    main()