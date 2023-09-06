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
        run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:nanopore\" %(ref)s %(sample)s.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        #run_cmd("minimap2 -ax map-ont %(ref)s %(sample)s.group.fastq | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_coverage_pf_mean.txt" % vars(args))
        
    
    with open("bam_list.txt","w") as O:
        for s in samples:
            O.write("%s.bam\n" % (s))

    run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 --min-coverage 50 --min-base-quality %(min_base_qual)s --gvcf --gvcf-dont-use-chunk true | bcftools view -T %(bed)s | bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))
    run_cmd("bcftools query -f '%CHROM\t%POS[\t%DP]\n' combined.genotyped.vcf.gz > tmp.txt")

    run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined.genotyped.vcf.gz | bcftools view -i 'QUAL>30' | bcftools sort | bcftools norm -m - -Oz -o tmp.vcf.gz" % vars(args))
    run_cmd("bcftools view -v snps tmp.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o snps.vcf.gz" % vars(args))
    run_cmd("tabix snps.vcf.gz" % vars(args))
    run_cmd("gunzip -c snps.vcf.gz > snps.vcf" % vars(args))
    run_cmd("bcftools view -v indels tmp.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o indels.vcf.gz" % vars(args))
    run_cmd("tabix indels.vcf.gz" % vars(args))

    run_cmd("awk '{gsub(/NC_000001.11/, \"1\"); gsub(/NC_000003.12/, \"3\"); gsub(/NC_000004.12/, \"4\"); gsub(/NC_000009.12/, \"9\"); gsub(/NC_000011.10/, \"11\"); gsub(/NC_000023.11/, \"X\"); print;}' snps.vcf > snps_num.vcf" % vars(args))
    run_cmd("SnpSift annotate %(clinVar)s snps_num.vcf > snps_num_clinvar_GRCh38.vcf" % vars(args))
    run_cmd(r"bcftools query snps_num_clinvar_GRCh38.vcf -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%RS\t%CLNDN\n]' > combined_genotyped_filtered_formatted.snps.txt")   



# Set up the parser
parser = argparse.ArgumentParser(description='Amplicon sequencing analysis script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='CSV file with the "sample" column for sample IDs',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
parser.add_argument('--gff',type=str,help='GFF file',required=True)
parser.add_argument('--bed',type=str,help='BED file with genes/amplicon locations',required=True)
parser.add_argument('--clinVar',type=str,help='ClinVar SNP annotation file',required=True)
parser.add_argument('--min-base-qual',default=30,type=int,help='Minimum base quality to use by freebayes')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
