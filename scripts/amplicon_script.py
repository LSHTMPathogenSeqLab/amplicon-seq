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

    run_cmd("demultiplex_fastq.py --R1 %(read1)s --R2 %(read2)s --index %(index_file)s" % vars(args))

    for sample in samples:
        args.sample = sample
        run_cmd("fastqc %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz" % vars(args))
        if args.trim:
            run_cmd("trimmomatic PE %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz %(sample)s_1.trimmed.fastq.gz %(sample)s_1.untrimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz %(sample)s_2.untrimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:36 2> %(sample)s.trimlog" % vars(args))
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.trimmed.fastq.gz %(sample)s_2.trimmed.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        else:
            run_cmd("bwa mem -t 10 -R \"@RG\\tID:%(sample_prefix)s%(sample)s\\tSM:%(sample_prefix)s%(sample)s\\tPL:Illumina\" %(ref)s %(sample)s_1.fastq.gz %(sample)s_2.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))

        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_region_coverage.txt" % vars(args))
        run_cmd("sambamba depth base %(sample)s.bam > %(sample)s.coverage.txt" % vars(args))
        # run_cmd("freebayes -f %(ref)s -t %(bed)s --haplotype-length -1 %(sample)s.bam --min-base-quality %(min_base_qual)s | bcftools view -Oz -o %(sample)s.vcf.gz" % vars(args))
    
    with open("bam_list.txt","w") as O:
        for s in samples:
            O.write("%s.bam\n" % (s))

    
    args.min_adf_cmd = f" | amplicon_setGT.py --min-adf {args.min_adf}" if args.min_adf else ""
    run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 --min-coverage 50 --min-base-quality %(min_base_qual)s --gvcf --gvcf-dont-use-chunk true | bcftools view -T %(bed)s %(min_adf_cmd)s | bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))
    run_cmd(r"bcftools query -f '%CHROM\t%POS[\t%DP]\n' combined.genotyped.vcf.gz > tmp.txt")

    run_cmd("bcftools filter -i 'FMT/DP>10' -S . combined.genotyped.vcf.gz | bcftools view -i 'QUAL>30' | bcftools sort | bcftools norm -m - -Oz -o tmp.vcf.gz" % vars(args))
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
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
