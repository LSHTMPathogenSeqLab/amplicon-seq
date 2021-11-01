import sys
import argparse
import subprocess as sp
import csv
import fastq2matrix as fm
from fastq2matrix import run_cmd





def main(args):
    samples = []
    for f in args.index_files:
        for row in csv.DictReader(open(f)):
            if row["sample"] in samples:
                sys.stderr.write(f"Warning! You have a duplicate sample name: {row['sample']}\n")
            samples.append(row["sample"])

    with open("vcf_files.txt","w") as O:
        for s in samples:
            O.write("%s.freebayes.vcf")
            O.write("%s.gatk.vcf")
    for sample in samples:
        args.sample = sample
        run_cmd("naive_variant_caller.py --ref %(ref)s --bam %(sample)s.bam --sample %(sample)s --min-af %(min_sample_af)s --vcf-file-list vcf_list.txt | bcftools view -Oz -o %(sample)s.vcf.gz" % vars(args))
        run_cmd("tabix -f %(sample)s.vcf.gz" % vars(args))
    with open("vcf_list.txt","w") as O:
        for s in samples:
            O.write("%s.vcf.gz\n" % (s))
    run_cmd("bcftools merge -l vcf_list.txt -Oz -o combined.vcf.gz" )
    run_cmd(r"bcftools query -f '%CHROM\t%POS[\t%DP]\n' combined.vcf.gz > tmp.txt")

    run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined.vcf.gz | bcftools sort -Oz -o tmp.vcf.gz" % vars(args))
    run_cmd("bcftools view -v snps tmp.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o snps.vcf.gz" % vars(args))
    run_cmd("tabix snps.vcf.gz" % vars(args))
    run_cmd("bcftools view -v indels tmp.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o indels.vcf.gz" % vars(args))
    run_cmd("tabix indels.vcf.gz" % vars(args))
    run_cmd(r"bcftools query snps.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\n]' > combined_genotyped_filtered_formatted.snps.txt")    
    run_cmd(r"bcftools query snps.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%TBCSQ\n]' > combined_genotyped_filtered_formatted.snps.trans.txt")
    run_cmd(r"bcftools query indels.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\n]' > combined_genotyped_filtered_formatted.indels.txt")
    run_cmd(r"bcftools query indels.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%TBCSQ\n]' > combined_genotyped_filtered_formatted.indels.trans.txt")


    
    pos_info = {}
    for l in open(args.position_info):
        row = l.strip().split()
        pos_info[(row[0],row[1])] = (row[2],row[3])
    with open("depth_info.txt", "w") as O:
        O.write("chrom\tpos\tgene\tcsq\t%s\n" % "\t".join(samples))
        for l in open("tmp.txt"):
            row = l.strip().split()
            if (row[0],row[1]) in pos_info:
                p = pos_info[(row[0],row[1])]
                row.insert(2,p[1])
                row.insert(2,p[0])
            else:
                row.insert(2,"NA")
                row.insert(2,"NA")
            O.write("%s\n" % ("\t".join(row)))




# Set up the parser
parser = argparse.ArgumentParser(description='Amplicon sequencing analysis script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-files',nargs="+",type=str,help='CSV files containing fields "Sample,I1,I2"',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
parser.add_argument('--gff',type=str,help='GFF file',required=True)
parser.add_argument('--bed',type=str,help='BED file with genes/amplicon locations',required=True)
parser.add_argument('--position-info',type=str,help='Position info file',required=True)
parser.add_argument('--trim',action="store_true",help='Perform triming')
parser.add_argument('--trim-qv',default=20,type=int,help='Quality value to use in the sliding window analysis')
parser.add_argument('--min-base-qual',default=30,type=int,help='Minimum base quality to use by freebayes')
parser.add_argument('--min-adf',type=float,help='Set a minimum frequency for a mixed call')
parser.add_argument('--sample-prefix',default="",help=argparse.SUPPRESS)
parser.add_argument('--vcf-files',type=str,nargs="+",help='VCF files with positions to include (optional)')
parser.add_argument('--min-variant-qual',default=30,type=int,help='Quality value to use in the sliding window analysis')
parser.add_argument('--min-sample-af',default=0.05,type=float,help='Quality value to use in the sliding window analysis')
parser.add_argument('--per-sample-only',action="store_true",help='Perform triming')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
