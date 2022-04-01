#! /usr/bin/env python
import sys
import subprocess as sp
import argparse
from uuid import uuid4 
from collections import defaultdict
import glob
import os

class hit:
    def __init__(self,seq,chrom,start,end) -> None:
        self.seq = seq
        self.chrom = chrom
        self.start = start
        self.end = end


def main(args):
    sp.call("samtools faidx %s" % args.ref,shell=True)
    tmp = str(uuid4())
    tmp_fasta = f"{tmp}.fasta"
    with open(tmp_fasta,"w") as O:
        O.write(f">forward\n{args.forward}\n")
        O.write(f">reverse\n{args.reverse}\n")

    sp.call("blastn -task blastn -subject %s -query %s -outfmt 6 > %s.result" % (args.ref,tmp_fasta,tmp),shell=True)
    
    results = defaultdict(list)
    for l in open("%s.result" % tmp):
        row = l.rstrip().split()
        if row[0] in results: continue
        start = int(row[8])
        end = int(row[9])
        if start>end:
            start,end = end,start
        results[row[0]].append(hit(row[0],row[1],start,end))

    for fp in results["forward"]:
        for rp in results["reverse"]:
            if fp.chrom != rp.chrom:
                continue
            if abs(fp.start - rp.start)>args.max_dist:
                continue
            if fp.start<rp.end:
                start = fp.start
                end = rp.end
                print(vars(fp),"<---->",vars(rp))
            else:
                start = rp.start
                end = fp.end
                print(vars(rp),"<---->",vars(fp))
            for l in sp.Popen("samtools faidx %s %s:%s-%s" % (args.ref,fp.chrom,start,end),shell=True,stdout=sp.PIPE).stdout:
                sys.stdout.write(l.decode())
    
    for f in glob.glob(tmp):
        os.remove(f)
    


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--forward','-F',type=str,help='Forward primer')
parser.add_argument('--reverse','-R',type=str,help='Reverse primer')
parser.add_argument('ref',type=str,help='Reference genome fasta file')
parser.add_argument('--max-dist',type=int,default=5000,help='Reference genome fasta file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
