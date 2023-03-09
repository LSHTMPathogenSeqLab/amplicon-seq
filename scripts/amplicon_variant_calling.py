#! /usr/bin/env python
import sys
import argparse
import fastq2matrix as fm
import csv 

def main(args):
    samples = []
    reader = csv.DictReader(open(args.index_file))
    if "sample" not in reader.fieldnames:
        reader = csv.DictReader(open(args.index_file,encoding='utf-8-sig'))
    for row in reader:
        if row["sample"]=="": continue
        samples.append(row["sample"])

    fm.bwa_index(args.ref)
    fm.faidx(args.ref)

    cmd = "demultiplex_fastq.py --R1 %(read1)s --R2 %(read2)s --index %(index_file)s" % vars(args)
    if args.search_flipped_index:
        cmd += " --search-flipped-index"
    fm.run_cmd(cmd)

    for sample in samples:
        args.sample = sample
        fm.run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))

        fm.run_cmd("samtools index %(sample)s.bam" % vars(args))

    with open("bam_list.txt","w") as O:
        for s in samples:
            O.write("%s.bam\n" % (s))

    fm.run_cmd("freebayes -f %(ref)s -L bam_list.txt --haplotype-length -1 --min-coverage 50 --min-base-quality %(min_base_qual)s  --gvcf --gvcf-dont-use-chunk true | bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))
    fm.run_cmd(r"bcftools view -c 1 combined.genotyped.vcf.gz | bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\n]' > combined_genotyped_filtered_formatted.snps.txt")    
    fm.run_cmd(r"bcftools query -f '%CHROM\t%POS[\t%DP]\n' combined.genotyped.vcf.gz > depth_info.txt")

parser = argparse.ArgumentParser(description='Amplicon variant calling script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='CSV file containing fields "Sample,I1,I2"',required=True)
parser.add_argument('--read1',type=str,help='Forward read file',required=True)
parser.add_argument('--read2',type=str,help='Reverse read file',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
parser.add_argument('--min-base-qual',default=30,type=int,help='Minimum base quality to use by freebayes')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)