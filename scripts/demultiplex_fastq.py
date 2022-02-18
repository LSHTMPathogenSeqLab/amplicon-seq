#! /usr/bin/env python
import sys
import argparse
import csv
import subprocess
from tqdm import tqdm
import gzip

def num_mismatches(str_x,str_y):
    mismatches = 0
    for char_x,char_y in zip(str_x,str_y):
        if char_x != char_y:
            mismatches+=1
    return mismatches

def main(args):
    index = {}
    reader = csv.DictReader(open(args.index))
    if "sample" not in reader.fieldnames:
        reader = csv.DictReader(open(args.index,encoding='utf-8-sig'))
    for row in reader:
        if row["sample"]=="": continue
        index[(row["I1"],row["I2"])] = row["sample"]
        index[(row["I2"],row["I1"])] = row["sample"]
        index_size = len(row["I1"])
    R1 = gzip.open(args.R1)
    R2 = gzip.open(args.R2)
    R1_FILES = {s:gzip.open("%s_1.fastq.gz" % s,"wb") for s in index.values()}
    R2_FILES = {s:gzip.open("%s_2.fastq.gz" % s,"wb") for s in index.values()}
    UNKNOWN1 = gzip.open("Unknown_1.fastq.gz","wb")
    UNKNOWN2 = gzip.open("Unknown_2.fastq.gz","wb")
    for l in subprocess.Popen("gunzip -c %s | wc -l" % args.R1, shell=True,stdout=subprocess.PIPE).stdout:
        num_lines = int(l.strip())
    reads_assigned = {s:0 for s in index.values()}
    reads_assigned["not_assigned"] = 0
    for i in tqdm(range(int(num_lines/4))):
        data1 = [R1.readline() for _ in range(4)]
        data2 = [R2.readline() for _ in range(4)]

        indices = (data1[1].decode()[:index_size],data2[1].decode()[:index_size])
        if indices in index:
            sample = index[indices]
            reads_assigned[sample]+=1
            for i in range(4):
                R1_FILES[sample].write(data1[i])
                R2_FILES[sample].write(data2[i])

        else:
            mismatches = {}
            possible_matches = {}
            for i1,i2 in index:
                mismatch_r1 = num_mismatches(i1,indices[0])
                mismatch_r2 = num_mismatches(i2,indices[1])
                mismatches[index[i1,i2]] = (mismatch_r1,mismatch_r2)
                if mismatch_r1<=args.max_mismatches and mismatch_r2<=args.max_mismatches:
                    possible_matches[index[i1,i2]] = (i1,i2)
            if len(list(possible_matches.keys()))==1:
                sample = list(possible_matches.keys())[0]
                reads_assigned[sample]+=1
                for i in range(4):
                    R1_FILES[sample].write(data1[i])
                    R2_FILES[sample].write(data2[i])
                # R2_FILES[sample].write("".join(data2))
            else:
                reads_assigned["not_assigned"]+=1
                for i in range(4):
                    UNKNOWN1.write(data1[i])
                    UNKNOWN2.write(data2[i])

    print(reads_assigned)



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--R1','-1',type=str,help='NGS Platform',required=True)
parser.add_argument('--R2','-2',type=str,help='NGS Platform',required=True)
parser.add_argument('--index',type=str,help='NGS Platform',required=True)
parser.add_argument('--max-mismatches',default=1,type=int,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
