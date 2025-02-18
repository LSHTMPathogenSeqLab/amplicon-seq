#! /usr/bin/env python
import argparse
import subprocess as sp
from uuid import uuid4
import csv
from collections import defaultdict
import sys
import gzip 
from prettytable import PrettyTable

def find_barcodes_in_read(hits,read_size,barcodes,edge_size=10):
    found_barcodes = []
    for f,r in barcodes:
        forward_found = any([f==d[2] for d in hits if is_on_edge(d,read_size,edge_size)])
        reverse_found = any([r==d[2] for d in hits if is_on_edge(d,read_size,edge_size)])
        if forward_found or reverse_found:
            found_barcodes.append(barcodes[(f,r)])
    return found_barcodes

def is_on_edge(hit,read_size,edge_size):
    if int(hit[5])<edge_size or int(hit[4])>read_size-edge_size:
        return True
    else:
        return False

def read_fastq(fastq):
    """Yield (name, and sequence) tuples from a FASTQ file."""
    is_gzip = fastq.endswith(".gz")
    with gzip.open(fastq) if is_gzip else open(fastq) as fh:
        while True:
            name = fh.readline().strip()
            if not name:
                break
            if is_gzip:
                name = name.decode()
            name = name.split()[0][1:]
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            yield name, seq


def main(args):
    if args.log_prefix:
        tmp = args.log_prefix
    else:
        tmp = "tmp" #str(uuid4())
    LOG = open(f"{tmp}.log","w")
    barcode_fasta = f"{tmp}.barcodes.fasta"
    seqkit_output = f"{tmp}.seqkit.txt"

    id2barcodes = {}
    barcodes2id = {}
    for row in csv.DictReader(open(args.barcodes)):
        id2barcodes[row['id']] = (row['forward'],row['reverse'])
        barcodes2id[(row['forward'],row['reverse'])] = row['id']
    
    with open(barcode_fasta,"w") as O:
        for p in barcodes2id:
            O.write(f">{barcodes2id[p]}_forward\n{p[0]}\n")
            O.write(f">{barcodes2id[p]}_reverse\n{p[1]}\n")

    sys.stderr.write("Searching for barcodes in reads\n")
    cmd = f"seqkit locate --max-mismatch {args.max_mismatch} -f {barcode_fasta} {args.fastq} > {seqkit_output}"
    sp.call(cmd,shell=True)

    hits = defaultdict(list)
    for l in open(seqkit_output):
        row = l.strip().split()
        if row[0]=="seqID": continue
        hits[row[0]].append(row)


    sys.stderr.write("Loading all read names and length\n")
    read_lengths = {}
    read_names = set()
    for name,seq in read_fastq(args.fastq):
        read_lengths[name] = len(seq)
        read_names.add(name)

    barcode_assignments = defaultdict(set)
    assigned_reads = set()
    for readname,hitlines in hits.items():
        barcode_matches = find_barcodes_in_read(hitlines,read_lengths[readname],barcodes2id,args.edge_size)
        if len(barcode_matches)==1:
            barcode_assignments[barcode_matches[0]].add(readname)
            assigned_reads.add(readname)
            LOG.write(f"{readname}\t{barcode_matches[0]}\t{read_lengths[readname]}\n")
        elif len(barcode_matches)==0:
            LOG.write(f"{readname}\tno_barcode_matches\t{read_lengths[readname]}\n")
        else:
            LOG.write(f"{readname}\tmultiple_barcode_matches\t{read_lengths[readname]}\n")

    unassigned_reads = read_names - assigned_reads
    for readname in unassigned_reads:
        barcode_assignments['unassigned'].add(readname)
    

    sys.stderr.write("Writing reads to files\n")
    for bid in barcode_assignments:
        tmp_barcode_id_file = f"{tmp}.{bid}.ids"
        with open(tmp_barcode_id_file,"w") as O:
            O.write("\n".join(barcode_assignments[bid]))
        sys.stderr.write(f"Writing {len(barcode_assignments[bid])} reads for {bid}\n")
        sp.call(f"seqkit grep -f {tmp_barcode_id_file} {args.fastq} > {bid}.fastq 2> /dev/null",shell=True)

    LOG.close()

    table = [[bid,len(barcode_assignments[bid]),round(len(barcode_assignments[bid])/len(read_names)*100,2)] for bid in barcode_assignments]
    tab = PrettyTable(["barcode","N","%"])
    tab.add_rows(table)
    # tab.add_column('col 5', [-123, 43], align='r', valign='t')
    print(tab)
    
    
parser = argparse.ArgumentParser(description='add required annotations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fastq',type=str,help='Reference file (lofreq required)',required = True)
parser.add_argument('--barcodes',type=str,help='Sample name (lofreq required)',required = True)
parser.add_argument('--max-mismatch',type=int,default=0,help='Sample name (lofreq required)')
parser.add_argument('--edge-size',type=int,default=10,help='Sample name (lofreq required)')
parser.add_argument('--log-prefix',type=str,default=None,help='Sample name (lofreq required)')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)