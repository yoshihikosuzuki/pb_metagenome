#!/usr/bin/env python
#coding: utf-8

import subprocess
import h5py
from collections import defaultdict


def replace_reads(bax_file, hns):

    with h5py.File(bax_file, 'a') as f:
        base_call = f["/PulseData/BaseCalls/Basecall"]
        hole_number = f["/PulseData/BaseCalls/ZMW/HoleNumber"]
        num_event = f["/PulseData/BaseCalls/ZMW/NumEvent"]

        for i in range(len(hole_number)):
            if hole_number[i] in hns:
                offset = sum(num_event[:i])
                length = num_event[i]
                for j in range(offset, offset + length):
                    base_call[j] = 65

                    
def exist_in_bax(bax_file, hns):

    with h5py.File(bax_file, 'r') as f:
        hole_number = f["/PulseData/BaseCalls/ZMW/HoleNumber"]

        for i in range(len(hole_number)):
            if hole_number[i] in hns:
                return True

    return False

                    
if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="replace reads to AAA... in raw .bax.h5 files")
    parser.add_argument('bax_fofn', help="file of .bax.h5 file names with absolute paths")
    parser.add_argument('subread_list', help="file of the list of subread headers to replace")
    args = parser.parse_args()

    movie_hole = defaultdict(set)
    with open(args.subread_list, 'r') as f:
        for line in f:
            movie, hole = line.strip().split('/')[0:2]
            movie_hole[movie].add(int(hole))

    with open(args.bax_fofn, 'r') as f:
        for line in f:
            fpath = line.strip()
            movie = fpath.split('/')[-1].split('.')[0]

            if movie not in movie_hole:
                continue
            
            if exist_in_bax(fpath, movie_hole[movie]):

                subprocess.call(' '.join(["cp", fpath, fpath + ".old"]), shell = True)
                replace_reads(fpath, movie_hole[movie])
