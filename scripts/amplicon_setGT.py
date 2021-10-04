#! /usr/bin/env python
import sys
from tqdm import tqdm
import argparse
import re
def main(args):
    for l in tqdm(sys.stdin):
        line = l.strip().replace("|","/")
        something_changed = False
        row = line.split()
        if l[0]=="#":
            sys.stdout.write(l)
            continue
        alleles = [row[3]]+row[4].split(",")
        if len(alleles)>9: continue

        format = row[8].split(":")
        if "AD" not in format:
            continue 
        ad_i = format.index("AD")
        for i in range(9,len(row)):
            fmt = row[i].split(":")
            ad = [int(x) if x!="." else 0 for x in fmt[ad_i].split(",")]
            if sum(ad)==0:
                continue
            adf = [d/sum(ad) for d in ad]
            if max(adf)<1-args.min_adf:
                ordered_idx = sorted(sorted(range(len(adf)), key=lambda k: adf[k],reverse=True)[:2])
                # print(l)
                # print(ordered_idx)
                # quit()
                new_gt = f"{ordered_idx[0]}/{ordered_idx[1]}"
                if args.debug and fmt[0]!=new_gt:
                    sys.stderr.write(f"{row[0]}\t{row[1]}\t{fmt[0]}\t{new_gt}\t{adf}\n")
                fmt[0] = new_gt
                something_changed = True
                row[i] = ":".join(fmt)
        
        
        if something_changed:
            sys.stdout.write("\t".join(row)+"\n")
        else:
            sys.stdout.write(l)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--min-adf',default=0.05,type=float,help='Fraction of coverage to assign major')
parser.add_argument('--debug',action="store_true",help='Fraction of coverage to assign major')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
