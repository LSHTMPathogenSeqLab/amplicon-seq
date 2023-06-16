#! /usr/bin/env pyhton
import sys
import argparse
import subprocess as sp
import csv
from uuid import uuid4

parser = argparse.ArgumentParser(description='Demultiplex plates')
parser.add_argument('-p', '--plate-layout', help='Input file', required=True)
parser.add_argument('-b', '--barcodes', help='Output file', required=True)
parser.add_argument('-f', '--fastq-dir', help='Fastq firectory', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', default=4, type=int)
parser.add_argument('-d', '--malaria-profiler-db', help='Malaria profiler db', default="amplicon")

args = parser.parse_args()

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
            cmds.append(f"malaria-profiler profile --platform nanopore --txt -1 {new_id}.fastq -p {new_id} --dir malaria-profiler-results --resistance_db {args.malaria_profiler_db}")
        writer = csv.DictWriter(temp,fieldnames=['id','forward','reverse'])
        writer.writeheader()
        writer.writerows(rows)
    sp.run(f"demux_nanopore_amplicon.py --fastq {args.fastq_dir}/{bc}.fastq.gz --barcodes {tmp_barcode_file} --max-mismatch 0 --edge-size 12 --log-prefix {bc}",shell=True)
    sp.run(f"rm {tmp_barcode_file}",shell=True)

run_file = f"{uuid4()}.sh"
with open(run_file,'w') as O:
    O.write("\n".join(cmds))
sp.run("mkdir malaria-profiler-results",shell=True)
sp.run(f"cat {run_file} | parallel -j {args.threads} --bar",shell=True)
sp.run(f"malaria-profiler collate",shell=True)
