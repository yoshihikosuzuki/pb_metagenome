#!/usr/bin/env python
#coding: utf-8

import re
from pbcore.io import FastaReader


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="filter out subreads mapped to human genome from fasta file")
    parser.add_argument('subread_fasta', help="subread file")
    parser.add_argument('filter_header', help="headers to be filtered out")
    args = parser.parse_args()

    filter_headers = set()
    with open(args.filter_header, 'r') as f:
        for line in f:
            filter_headers.add(re.search(r'(.*)\/.*', line.strip()).groups()[0])

    for subread in FastaReader(args.subread_fasta):
        if re.search(r'(.*)\/.*', subread.name).groups()[0] not in filter_headers:
            print ">" + subread.name + "\n" + subread.sequence
            continue
