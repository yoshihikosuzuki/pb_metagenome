#!/usr/bin/env python
#coding: utf-8

from interval import interval
from collections import defaultdict
import subprocess


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="calculate cover rate of subreads against reference genome")
    parser.add_argument('ovl_file', help="output file of LA4Falcon") 
    parser.add_argument('--subread_prefix', type=str, default="subreads", help="prefix of input subread fasta used for daligner")
    parser.add_argument('--cover_rate', type=float, default=0.5, help="threshold for cover rate")
    args = parser.parse_args()

    mapped_region = defaultdict(interval)
    subread_length = {}
    with open(args.ovl_file, 'r') as f:
        for line in f:
            data = line.strip().split(' ')
            subread_id = data[0]
            start, end, length = map(int, data[5:8])
            
            mapped_region[subread_id] = mapped_region[subread_id] | interval[start, end]
            if subread_id not in subread_length:
                subread_length[subread_id] = length

    filter_id = set()
    for subread_id, region in mapped_region.items():
        mapped_length = sum([i[0][1] - i[0][0] for i in region.components])
        #print subread_id, subread_length[subread_id], mapped_length, float(mapped_length) / subread_length[subread_id]
        if float(mapped_length) / subread_length[subread_id] >= args.cover_rate:
            filter_id.add(subread_id)

    subread_header = subprocess.check_output("DBshow -n " + args.subread_prefix, shell = True).split('\n')
    for subread_id in sorted(filter_id):
        #print subread_id, subread_length[subread_id], subread_header[int(subread_id)][1:]
        print subread_header[int(subread_id)][1:]
