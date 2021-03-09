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
        run_cmd("trimmomatic PE %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz %(sample)s_1.trimmed.fastq.gz %(sample)s_1.untrimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz %(sample)s_2.untrimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36" % vars(args))
        run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samtools sort -o %(sample)s.bam -" % vars(args))
        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s_flagstat.txt" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_coverage_aedes_mean.txt" % vars(args))
        run_cmd("samtools stats %(sample)s.bam > %(sample)s.stats" % vars(args)) # needed?
        run_cmd("gatk HaplotypeCaller -R %(ref)s -L %(bed)s -I %(sample)s.bam -O %(sample)s.gatk.vcf.gz -ERC BP_RESOLUTION -mbq 30" % vars(args)) # gatk -> g?

    with open("vcfout.list","w") as O:
        for s in samples:
            O.write("%s.gatk.vcf.gz\n" % (s))

    run_cmd("gatk CombineGVCFs -R %(ref)s -L %(bed)s --variant vcfout.list -O combined.gatk.vcf" % vars(args))
    run_cmd("gatk GenotypeGVCFs -R %(ref)s -L %(bed)s -V combined.gatk.vcf -O combined_genotyped.gatk.vcf.gz" % vars(args))

    run_cmd("bcftools filter -i 'FMT/DP>10 & FMT/GQ>30' -S . combined_genotyped.gatk.vcf.gz -Oz -o combined_genotyped.gatk.filtered.vcf.gz")
    run_cmd("bcftools index combined_genotyped.gatk.filtered.vcf.gz" % vars(args))

    # tranlation
    run_cmd("bcftools csq -p a -f %(ref)s -g %(gff)s combined_genotyped.gatk.filtered.vcf.gz -Oz -o combined_genotyped.gatk.filtered.trans.vcf.gz" % vars(args))
    run_cmd("bcftools index combined_genotyped.gatk.filtered.trans.vcf.gz" % vars(args))

    run_cmd("bcftools query combined_genotyped.gatk.filtered.vcf.gz -f '[%SAMPLE\t%POS\t%REF\t%ALT\t%QUAL\t%GQ\t%GT\t%TGT\t%DP\t%AD\n]' > combined_genotyped_filtered_formatted.gatk.txt")

    run_cmd("bcftools query combined_genotyped.gatk.filtered.trans.vcf.gz -f '[%SAMPLE\t%POS\t%REF\t%ALT\t%QUAL\t%GQ\t%GT\t%TGT\t%DP\t%AD\t%INFO/BCSQ\n]' > combined_genotyped_filtered_formatted.trans.gatk.txt")



# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='File suffix',required=True)
parser.add_argument('--read1',type=str,help='File suffix',required=True)
parser.add_argument('--read2',type=str,help='File suffix',required=True)
parser.add_argument('--ref',type=str,help='File suffix',required=True)
parser.add_argument('--gff',type=str,help='File suffix',required=True)
parser.add_argument('--bed',type=str,help='File suffix',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
