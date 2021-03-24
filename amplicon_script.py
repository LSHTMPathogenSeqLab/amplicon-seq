import sys
import argparse
import subprocess as sp
import csv
import fastq2matrix as fm

def run_cmd(cmd):
    sys.stderr.write("Running command:\n%s\n\n" % cmd)
    with open("/dev/null","w") as O:
        res = sp.call(cmd,shell=True,stderr=O,stdout=O)
    if res!=0:
        sys.exit("Error running last command, please check!\n")


def main(args):

    samples = []
    for row in csv.DictReader(open(args.index_file)):
        samples.append(row["sample"])

    fm.bwa_index(args.ref)
    fm.create_seq_dict(args.ref)
    fm.faidx(args.ref)

    run_cmd("demultiplex_fastq.py --R1 %(read1)s --R2 %(read2)s --index %(index_file)s" % vars(args))
    
    for sample in samples:
        args.sample = sample
        run_cmd("fastqc %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz" % vars(args))
        if args.trim:
            run_cmd("trimmomatic PE %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz %(sample)s_1.trimmed.fastq.gz %(sample)s_1.untrimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz %(sample)s_2.untrimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:36 2> %(sample)s.trimlog" % vars(args))
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samtools sort -o %(sample)s.bam -" % vars(args))
        else:
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz | samtools sort -o %(sample)s.bam -" % vars(args))

        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))

    with open("bam_list.txt","w") as O:
        for s in samples:
            O.write("%s.bam\n" % (s))

    run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 | bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))

    run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined.genotyped.vcf.gz | bcftools view -i 'QUAL>30' | bcftools sort | bcftools norm -m - -Oz -o tmp.vcf.gz" % vars(args))
    run_cmd("bcftools view -v snps tmp.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o tmp.snps.vcf.gz" % vars(args))
    run_cmd("tabix tmp.snps.vcf.gz" % vars(args))
    run_cmd("bcftools view -v indels tmp.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o tmp.indels.vcf.gz" % vars(args))
    run_cmd("tabix tmp.indels.vcf.gz" % vars(args))
    run_cmd("bcftools concat -a tmp.snps.vcf.gz tmp.indels.vcf.gz | bcftools sort -Oz -o combined.genotyped.csq.vcf.gz" % vars(args))
    run_cmd("bcftools index combined.genotyped.csq.vcf.gz" % vars(args))
    run_cmd(r"bcftools query combined.genotyped.csq.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%INFO/BCSQ\n]' > combined_genotyped_filtered_formatted.trans.gatk.txt")



# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='File suffix',required=True)
parser.add_argument('--read1',type=str,help='File suffix',required=True)
parser.add_argument('--read2',type=str,help='File suffix',required=True)
parser.add_argument('--ref',type=str,help='File suffix',required=True)
parser.add_argument('--gff',type=str,help='File suffix',required=True)
parser.add_argument('--bed',type=str,help='File suffix',required=True)
parser.add_argument('--trim',action="store_true",help='File suffix')
parser.add_argument('--sample-prefix',default="",help='File suffix')
parser.add_argument('--trim-qv',default=20,type=int,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
