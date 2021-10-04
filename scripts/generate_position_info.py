import sys
import argparse
import pyfastx
import subprocess as sp
import re

gene_lookup = {
    "transcript:CZT98080":"MSP2", 
    "transcript:CAD52497":"TRAP"
}

def main(args):
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
    sp.call(r"snpEff ann Plasmodium_falciparum_3D7 -no-upstream -no-downstream tmp.vcf | bcftools query -f '%CHROM\t%POS\t%ANN\n' > info.txt",shell=True)
    for l in open("info.txt"):
        row = l.strip().split()
        ann = row[2].split(",")[0].split("|")
        if "intergenic_region" in ann[1]:
            csq = ("-",ann[1])
        elif "missense" not in ann[1] and "stop" not in ann[1] and "synonymous" not in ann[1]:
            csq = (gene_lookup[ann[6]],ann[1])
        else:
            re_obj = re.search("p.([A-Za-z]+)([0-9]+)([A-Za-z\*]+)",ann[10])
            m = "%s%s" % (re_obj.group(1),re_obj.group(2))
            csq = (gene_lookup[ann[6]],m)
            # print(ann[6])
        print("\t".join([row[0],row[1],csq[0],csq[1]]))

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bed',type=str,help='')
parser.add_argument('--ref',type=str,help='')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)