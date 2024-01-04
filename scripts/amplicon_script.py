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
    cmd = "demultiplex_fastq.py --R1 %(read1)s --R2 %(read2)s --index %(index_file)s" % vars(args)
    if args.search_flipped_index:
        cmd += " --search-flipped-index"
    run_cmd(cmd)

    for sample in samples:
        args.sample = sample
        run_cmd("fastqc %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz" % vars(args))
        
        if(args.trim = True and args.samclip = True):
            run_cmd("trimmomatic PE %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz %(sample)s_1.trimmed.fastq.gz %(sample)s_1.untrimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz %(sample)s_2.untrimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:36 2> %(sample)s.trimlog" % vars(args))
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        if(args.trim = True and args.samclip = False):
            run_cmd("trimmomatic PE %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz %(sample)s_1.trimmed.fastq.gz %(sample)s_1.untrimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz %(sample)s_2.untrimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:36 2> %(sample)s.trimlog" % vars(args))
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samtools sort -o %(sample)s.bam -" % vars(args))
        if(args.trim = False and args.samclip = False):
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samtools sort -o %(sample)s.bam -" % vars(args))
        if(args.trim = False and args.samclip = True):
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        
        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_region_coverage.txt" % vars(args))
        run_cmd("sambamba depth base %(sample)s.bam > %(sample)s.coverage.txt" % vars(args))
        run_cmd("freebayes -f %(ref)s -t %(bed)s %(sample)s.bam --haplotype-length -1> %(sample)s.freebayes.vcf" % vars(args))
        run_cmd("gatk HaplotypeCaller -R %(ref)s -L %(bed)s  -I %(sample)s.bam -O %(sample)s.gatk.vcf" % vars(args))

    if not args.per_sample_only:
        with open("vcf_list.txt","w") as O:
            for s in samples:
                O.write("%s.freebayes.vcf\n" % s)
                O.write("%s.gatk.vcf\n" % s)
        for sample in samples:
            args.sample = sample
            run_cmd("naive_variant_caller.py --ref %(ref)s --bam %(sample)s.bam --sample %(sample)s --min-af %(min_sample_af)s --vcf-file-list vcf_list.txt | bcftools view -Oz -o %(sample)s.vcf.gz" % vars(args))
            run_cmd("tabix -f %(sample)s.vcf.gz" % vars(args))
        with open("sample_vcf_list.txt","w") as O:
            for s in samples:
                O.write("%s.vcf.gz\n" % (s))
        run_cmd("bcftools merge -l sample_vcf_list.txt -Oz -o combined.vcf.gz" )
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

        pos_info = {}
        for l in open(args.position_info):
            row = l.strip().split()
            pos_info[(row[0],int(row[1]))] = (row[2],row[3])
        
        with open("depth_info.txt", "w") as O:
            O.write("chrom\tpos\tgene\tcsq\t%s\n" % "\t".join(samples))
            for chrom,pos in amplicon_positions:
                if (chrom,pos) in pos_info:
                    d = pos_info[(chrom,pos)]
                    O.write("%s\t%s\t%s\t%s\t%s\n" % (
                        chrom,pos,d[0],d[1],
                        "\t".join([str(dp[s].get((chrom,pos),0)) for s in samples])
                    ))

# Set up the parser
parser = argparse.ArgumentParser(description='Amplicon sequencing analysis script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='CSV file containing fields "Sample,I1,I2"',required=True)
parser.add_argument('--read1',type=str,help='Forward read file',required=True)
parser.add_argument('--read2',type=str,help='Reverse read file',required=True)
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
parser.add_argument('--search-flipped-index',action='store_true',help='NGS Platform')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
pasrer.add_argument('--samclip', action="store_true", help = 'Use samclip or not')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
