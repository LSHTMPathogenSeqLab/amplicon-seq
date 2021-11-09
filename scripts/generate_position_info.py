import sys
import argparse
import pyfastx
import subprocess as sp
import re
import os
import shutil
import glob
# gene_lookup = {
#     "transcript:CZT98080":"MSP2", 
#     "transcript:CAD52497":"TRAP"
# }

def crate_snpeff_db(name,ref,gff):
    tmp = glob.glob(f"{sys.base_prefix}/share/snpeff*")[0]
    if not os.path.isdir(f"{tmp}/data/"):
        os.mkdir(f"{tmp}/data/")
    d = f"{tmp}/data/{name}"
    if os.path.isdir(d):
        return
    os.mkdir(d)
    shutil.copy(ref,f"{d}/sequences.fa")
    shutil.copy(gff,f"{d}/genes.gff")
    with open(f"{tmp}/snpEff.config","a") as O:
        O.write(f"{name}.genome : {name}\n")
    sp.call(f"snpEff build -v -gff3 {name}",shell=True)

def get_user_input(x):
    i = input(f"I've encountered a gene called {x}. What should I call this?\n")
    return i

def main(args):
    gene_lookup = {}
    
    # for l in open(args.gene_lookup_file):
    #     row = l.strip().split()
    #     gene_lookup[row[0]] = row[1]
    args.species = "_".join(args.species).replace(" ","_")
    crate_snpeff_db("custom_"+args.species,args.ref,args.gff)
    fa = pyfastx.Fasta(args.ref)
    O = open("tmp.vcf","w")
    O.write("""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20210414
##source=freeBayes v1.3.5
##reference=../Reference_files/Aedes_aegypti_lvpagwg.AaegL5.fa
##contig=<ID=1,length=310827022>
##contig=<ID=2,length=474425716>
##contig=<ID=3,length=409777670>
##contig=<ID=MT,length=16790>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO    FORMAT  SAMPLE
""")
    for l in open(args.bed):
        row = l.strip().split()
        for i in range(int(row[1]),int(row[2])+1):
            ref_allele = fa.fetch(row[0],(i,i))
            n = ["A","C","G","T"]
            n.remove(ref_allele.upper())
            O.write(f"{row[0]}\t{i}\t.\t{ref_allele}\t{n[0]}\t.\t.\tDP=100\tGT\t1/1\n")
    sp.call(f"snpEff ann custom_{args.species} -no-upstream -no-downstream tmp.vcf | bcftools query -f '%CHROM\t%POS\t%ANN\n' > info.txt",shell=True)
    with open(args.outfile,"w") as O:
        for l in open("info.txt"):
            row = l.strip().split()
            ann_field = [a.split("|") for a in row[2].split(",")]

            if len(ann_field)==1 and ann_field[0][1]=="intergenic_region":
                ann = ann_field[0]
            else:
                possible_transcripts = set([a[6] for a in ann_field])
                if len(possible_transcripts.intersection(set(gene_lookup)))==0:
                    print(f"Position {row[0]}:{row[1]} is associated with {len(possible_transcripts)} transcripts")
                    tr = int(input(
                        "Which transcript shall I choose?\n" +  
                        "\n".join([f"{i}:\t{a[6]}\t{a[9]}\t{a[10]}" for i,a in enumerate(ann_field)]) + 
                        "\nPlease enter number: "

                        ))
                    tr = list(ann_field)[tr][6]
                    gene_name = input("And what should I call this gene?\n").strip()
                    gene_lookup[tr] = gene_name
                    print(gene_lookup)
                ann = [a for a in ann_field if a[6] in gene_lookup][0]
            if "intergenic_region" in ann[1]:
                csq = ("-",ann[1])
            elif "missense" not in ann[1] and "stop" not in ann[1] and "synonymous" not in ann[1]:
                if ann[6] not in gene_lookup:
                    gene_lookup[ann[6]] = get_user_input(ann[6])
                csq = (gene_lookup[ann[6]],ann[1])
            else:
                if ann[6] not in gene_lookup:
                    gene_lookup[ann[6]] = get_user_input(ann[6])
                re_obj = re.search("p.([A-Za-z]+)([0-9]+)([A-Za-z\*]+)",ann[10])
                m = "%s%s" % (re_obj.group(1),re_obj.group(2))
                csq = (gene_lookup[ann[6]],m)
                # print(ann[6])
            O.write("\t".join([row[0],row[1],csq[0],csq[1]])+"\n")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--species',nargs="+",type=str,help='',required=True)
parser.add_argument('--bed',type=str,help='',required=True)
parser.add_argument('--ref',type=str,help='',required=True)
parser.add_argument('--gff',type=str,help='',required=True)
parser.add_argument('--outfile',type=str,help='',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)