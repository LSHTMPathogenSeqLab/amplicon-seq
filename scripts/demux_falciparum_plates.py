#! /usr/bin/env python
import sys
import argparse
import subprocess as sp
import csv
from uuid import uuid4
from joblib import Parallel, delayed
from tqdm import tqdm
import re
import os


parser = argparse.ArgumentParser(description='Demultiplex plates',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p', '--plate-layout', help='Input file', required=True)
parser.add_argument('-b', '--barcodes', help='Output file', required=True)
parser.add_argument('-f', '--fastq-dir', help='Fastq firectory', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', default=1, type=int)
parser.add_argument('-j', '--jobs', help='Number of Jobs', default=4, type=int)
parser.add_argument('-d', '--malaria-profiler-db', help='Malaria profiler db', default="amplicon")
parser.add_argument('--depth',default="0,10",type=str,help='Minimum depth hard and soft cutoff specified as comma separated values')
parser.add_argument('--af',default="0.1,0.4",type=str,help='Minimum allele frequency hard and soft cutoff specified as comma separated values')
parser.add_argument('--strand',default="0,3",type=str,help='Minimum read number per strand hard and soft cutoff specified as comma separated values')


args = parser.parse_args()

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True
    return False

def run_cmd(cmd: str, desc=None, log: str=None) -> sp.CompletedProcess:
    programs = set([x.strip().split()[0] for x in re.split("[|&;]",cmd) if x!=""])
    missing = [p for p in programs if which(p)==False]
    if len(missing)>0:
        raise ValueError("Cant find programs: %s\n" % (", ".join(missing)))
    cmd = "/bin/bash -c set -o pipefail; " + cmd
    output = open(log,"w") if log else sp.PIPE
    result = sp.run(cmd,shell=True,stderr=output,stdout=output)
    if result.returncode != 0:
        raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,result.stderr.decode()))
    return result

# Read in plate layout
plate_layout = {}
for row in csv.DictReader(open(args.plate_layout,encoding='utf-8-sig')):
    plate_layout[row['Nanopore barcode ID']] = {k:f"{row['Well']}_{v}" for k,v in row.items() if k not in ('Nanopore barcode ID','Well')}

# read in barcodes
barcodes = {}
for row in csv.DictReader(open(args.barcodes,encoding='utf-8-sig')):
    barcodes[row['id']] = (row['forward'], row['reverse'])

# Demultiplex

cmds = []
for bc in plate_layout:
    # write temp file
    tmp_barcode_file = f"{uuid4()}.csv"
    with open(tmp_barcode_file, 'w') as temp:
        rows = []
        for key,val in plate_layout[bc].items():
            new_id = f"{key}_{val}"
            rows.append({'id':new_id, 'forward':barcodes[key][0], 'reverse':barcodes[key][1]})
        writer = csv.DictWriter(temp,fieldnames=['id','forward','reverse'])
        writer.writeheader()
        writer.writerows(rows)
    sp.run(f"demux_nanopore_amplicon.py --fastq {args.fastq_dir}/{bc}.fastq.gz --barcodes {tmp_barcode_file} --max-mismatch 0 --edge-size 12 --log-prefix {bc}",shell=True)
    sp.run(f"rm {tmp_barcode_file}",shell=True)

for bc in plate_layout:
    for key,val in plate_layout[bc].items():
        new_id = f"{key}_{val}"
        if os.path.isfile(f"{new_id}.fastq"):
            cmds.append(f"malaria-profiler profile --no_species --platform nanopore --txt -1 {new_id}.fastq -p {new_id} --dir malaria-profiler-results --resistance_db {args.malaria_profiler_db} --depth {args.depth} --af {args.af} --strand {args.strand} --threads {args.threads}")

sp.run("mkdir malaria-profiler-results",shell=True)
parallel = Parallel(n_jobs=args.jobs, return_as='generator')
[r for r in tqdm(parallel(delayed(run_cmd)(cmd) for cmd in cmds),total=len(cmds),desc="Running jobs")]
sp.run(f"malaria-profiler collate --dir malaria-profiler-results/",shell=True)
