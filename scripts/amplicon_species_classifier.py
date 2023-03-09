#! /usr/bin/env python
import sys
import argparse

import os
import statistics as stats
from tqdm import tqdm
from uuid import uuid4
import subprocess
from os.path import expanduser
import csv 

_VERSION = "0.0.1"

def filecheck(filename):
    """
    Check if file is there and quit if it isn't
    """
    if filename=="/dev/null":
        return filename
    elif not os.path.isfile(filename):
        sys.stderr.write("Can't find %s\n" % filename)
        exit(1)
    else:
        return filename

def run_cmd(cmd,verbose=1,target=None,terminate_on_error=True):
    """
    Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
    """
    if target and filecheck(target): return True
    cmd = "set -u pipefail; " + cmd
    if verbose>0:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)

    p = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()

    if terminate_on_error is True and p.returncode!=0:
        raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,stderr.decode()))

    if verbose>1:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)

    return (stdout.decode(),stderr.decode())


def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])



def get_canonical_kmer(kmer):
    t = {"A":"1","C":"2","T":"3","G":"4"}
    rkmer = revcom(kmer)
    nkmer = int("".join([t[x] for x in list(kmer)]))
    nrkmer = int("".join([t[x] for x in list(rkmer)]))
    return kmer if nkmer<nrkmer else rkmer


def check_for_kmers(kmer_list_file,read1,read2=None,kmer_table=None):

    kmer_dict = {}
    for l in open(kmer_list_file):
        row  = l.strip().split("\t")
        kmer_dict[get_canonical_kmer(row[0])] = row[1]


    tmp = str(uuid4())
    fastq_files = f"-file {read1}"
    if read2:
        fastq_files = fastq_files + f",{read2}"
    run_cmd("dsk %s -out %s" % (fastq_files,tmp))
    run_cmd("dsk2ascii -file %s.h5 -out %s.kmers.txt" % (tmp,tmp))

    file_kmers = {}
    for l in tqdm(open(tmp+".kmers.txt")):
        row = l.strip().split()
        if row[0] in kmer_dict:
            file_kmers[row[0]] = int(row[1])

    os.remove("%s.h5" % tmp)
    os.remove("%s.kmers.txt" % tmp)

    kmer_support = []
    species_set = set()
    for kmer,species in kmer_dict.items():
        num = file_kmers.get(kmer,0)
        kmer_support.append({"kmer":kmer,"species":species,"num":num})
        species_set.add(species)

    if kmer_table:
        with open(kmer_table,"w") as O:
            for k in kmer_support:
                O.write("%(kmer)s\t%(species)s\t%(num)s\n" % k)
    print(kmer_support)
    #### Test for coluzzi ###
    if len([x for x in kmer_support if x["species"]=="coluzzi" and x["num"]>0])>0:
        kmer_support = [x for x in kmer_support if x["species"]!="coluzzi"]
        for k in kmer_support:
            if k["species"]=="gambiae_complex":
                k["species"] = "coluzzi"
        
    else:
        kmer_support = [x for x in kmer_support if x["species"]!="gambiae"]
        for k in kmer_support:
            if k["species"]=="gambiae_complex":
                k["species"] = "gambiae"
    ###       End test     ###   

    species_support = []
    for s in set([x["species"] for x in kmer_support]):
        support = [x["num"] for x in kmer_support if x["species"]==s]

        if len(support)==0:
            continue
        if len([x for x in support if x!=0])<len(support)/2:
            continue
        mean = stats.mean(support)
        # std = stats.stdev(support)
        # species_support.append({"species":s,"mean":mean,"std":std})
        species_support.append({"species":s,"mean":mean})

    return species_support


def get_db(update = False):
    db_dir = f"{expanduser('~')}/.vector-profiler"
    kmer_file = f"{db_dir}/kmers.txt"
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)
    if not os.path.isfile(kmer_file) or update==True:
        import urllib.request
        sys.stderr.write("Downloading kmers\n")
        urllib.request.urlretrieve("https://raw.githubusercontent.com/LSHTMPathogenSeqLab/amplicon-seq/main/db/kmers.txt",kmer_file)
    
    return kmer_file

def main_update_db(args):
    get_db(update = True)

def stringify(x):
    try:
        return "%.2f" % x
    except:
        return "%s" % x

def main_collate(args):
        if args.samples:
            samples = [x.rstrip() for x in open(args.samples).readlines()]
        else:
            samples = [x.replace(".species.txt","") for x in os.listdir(args.dir) if x[-len(".species.txt"):]==".species.txt"]

        rows = []
        for s in samples:
            tmp = []
            for r in csv.DictReader(open(filecheck(f"{args.dir}/{s}.species.txt")),delimiter="\t"):
                tmp.append(r)
            combined_row = {"sample":s}
            for col in r:
                combined_row[col] = ";".join([str(x[col]) for x in tmp])
            rows.append(combined_row)
                

        with open(args.outfile,"w") as O:
            writer = csv.DictWriter(O,fieldnames = list(rows[0]),delimiter=args.sep)
            writer.writeheader()
            writer.writerows(rows)

def main_speciate(args):
    kmer_file  = get_db()

    kmer_table_file = f"{args.prefix}.kmer_table.txt"
    results = check_for_kmers(kmer_file,args.read1,args.read2,kmer_table=kmer_table_file)
    results = sorted(results,key=lambda x: x["mean"],reverse=True)
    with open(args.prefix+".species.txt","w") as O:
        O.write("\t".join(results[0].keys())+"\n")
        for res in results:
            O.write("\t".join([stringify(x) for x in res.values()])+"\n")



parser = argparse.ArgumentParser(description='NTM-Profiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--version', action='version', version="NTM-Profiler version %s" % _VERSION)
subparsers = parser.add_subparsers(help="Task to perform")

# Profile #
parser_sub = subparsers.add_parser('speciate', help='Run whole profiling pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
input=parser_sub.add_argument_group("Input options")
input.add_argument('--read1','-1',help='First read file')
input.add_argument('--read2','-2',help='Second read file')


output=parser_sub.add_argument_group("Output options")

output.add_argument('--prefix','-p',help='Output file with results generated',required = True)
parser_sub.add_argument('--version', action='version', version="NTM-Profiler version %s" % _VERSION)
parser_sub.set_defaults(func=main_speciate)


# Collate results #
parser_sub = subparsers.add_parser('collate', help='Collate results form multiple samples together', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--outfile','-o',help='Output file name',required = True)
parser_sub.add_argument('--samples',help='File with samples (one per line)')
parser_sub.add_argument('--dir','-d',default=".",help='Storage directory')
parser_sub.add_argument('--sep',default="\t",help='Field seperator')
parser_sub.add_argument('--version', action='version', version="NTM-Profiler version %s" % _VERSION)
parser_sub.set_defaults(func=main_collate)

# Update database #
parser_sub = subparsers.add_parser('update_db', help='Collate results form multiple samples together', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--version', action='version', version="NTM-Profiler version %s" % _VERSION)
parser_sub.set_defaults(func=main_update_db)




args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)