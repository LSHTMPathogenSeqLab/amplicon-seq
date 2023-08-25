#! /usr/bin/env python
import sys
import argparse
import subprocess as sp
import csv
import fastq2matrix as fm
from fastq2matrix import run_cmd
from collections import defaultdict
import gzip



def main(args):

    samples = []
    reader = csv.DictReader(open(args.index_file))
    if "sample" not in reader.fieldnames:
        reader = csv.DictReader(open(args.index_file,encoding='utf-8-sig'))
        for row in reader:
            if row["sample"]=="": continue
            samples.append(row["sample"])

    fm.bwa_index(args.ref)
    fm.create_seq_dict(args.ref)
    fm.faidx(args.ref)
    
    for sample in samples:
        args.sample = sample
        run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:nanopore\" %(ref)s %(sample)s.fastq | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_coverage_pf_mean.txt" % vars(args))
        
        with open("bam_list.txt","w") as O:
            for s in samples:
                O.write("%s.bam\n" % (s))
                
        run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 --min-coverage 50 --min-base-quality %(min_base_qual)s --gvcf --gvcf-dont-use-chunk true | bcftools view -T %(bed)s | bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))
        run_cmd(r"bcftools query -f '%CHROM\t%POS[\t%DP]\n' combined.genotyped.vcf.gz > tmp.txt")
        run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined.genotyped.vcf.gz | bcftools view -i 'QUAL>30' | bcftools sort | bcftools norm -m - -Oz -o tmp.vcf.gz" % vars(args))
        run_cmd("bcftools view -v snps tmp.vcf.gz > snps.vcf" % vars(args))
        run_cmd("bgzip -c snps.vcf > snps.vcf.gz" % vars(args))
        run_cmd("tabix -f snps.vcf.gz" % vars(args))
        run_cmd("bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o snps.vcf.gz" % vars(args))
        run_cmd("tabix snps.vcf.gz" % vars(args))
        run_cmd(r"bcftools query snps.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%TBCSQ\n]' > combined_genotyped_filtered_formatted.snps.trans.txt")    

        
        bedlines = []
        amplicon_positions = []
        for l in open(args.bed):
            row = l.strip().split()
            bedlines.append(row)        
            for p in range(int(row[1]),int(row[2])):
                amplicon_positions.append((row[0],p))

        def overlap_bedlines(a,bedlines):
            overlaps = []
            for b in bedlines:
                if b[0]==a[0]:
                    overlap = max(0, min(int(a[2]), int(b[2])) - max(int(a[1]), int(b[1])))
                    if overlap>0:
                        overlaps.append([b[0],max(int(a[1]),int(b[1])),min(int(a[2]),int(b[2]))])
            return overlaps

        dp = defaultdict(dict)
        for s in samples:
            for l in gzip.open(f"{s}.per-base.bed.gz"):
                row = l.decode().strip().split()
                overlaps = overlap_bedlines(row,bedlines)
                if len(overlaps)>0:
                    for overlap in overlaps:
                        for pos in range(int(overlap[1]),int(overlap[2])):
                            dp[s][(row[0],pos)] = int(row[3])


# Set up the parser
parser = argparse.ArgumentParser(description='Amplicon sequencing analysis script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='samples.csv with the "sample" column for sample IDs; created by "demux_nanopore_plates.py" if using 96-well plates',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
parser.add_argument('--gff',type=str,help='GFF file',required=True)
parser.add_argument('--bed',type=str,help='BED file with genes/amplicon locations',required=True)
parser.add_argument('--min-base-qual',default=30,type=int,help='Minimum base quality to use by freebayes')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
