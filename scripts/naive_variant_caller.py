import sys
import argparse
import subprocess as sp
from collections import defaultdict
import re

def tryint(x):
    try:
        return int(x)
    except:
        return x

def recode_indels(ref,filtered_alleles):
    
    size = [int(re.search("([\+\-0-9]+)",a["allele"]).group(1)) if ("-" in a["allele"] or "+" in a["allele"]) else 0 for a in filtered_alleles]
    i = size.index(min(size))
    smallest = filtered_alleles[i]
    if "+" in smallest["allele"]:
        for a in filtered_alleles:
            r = re.match("([A-Z])",a["allele"])
            if r:
                pass
            r = re.search("([A-Z]+)\+[0-9]+([A-Z]+)",a["allele"])
            if r:
                a["allele"] = r.group(1)+r.group(2)
    elif "-" in smallest["allele"]:
        r = re.search("([A-Z]+)\-[0-9]+([A-Z]+)",smallest["allele"])
        extra_seq = r.group(2)
        ref = r.group(1)+r.group(2)
        for a in filtered_alleles:
            r = re.search("([A-Z])",a["allele"])
            if r:
                a["allele"] = a["allele"]+extra_seq
            r = re.search("([A-Z]+)\+[0-9]+([A-Z]+)",a["allele"])
            if r:
                a["allele"] = ref+r.group(2)
            r = re.search("([A-Z]+)(\-[0-9]+)([A-Z]+)",a["allele"])
            if r:
                a["allele"] = ref[0:int(r.group(2))]
    else:
        for a in filtered_alleles:
            r = re.match("([A-Z])",a["allele"])
            if r:
                pass
            r = re.search("([A-Z]+)\+[0-9]+([A-Z]+)",a["allele"])
            if r:
                a["allele"] = ref+r.group(2)
    
    return ref,filtered_alleles 

def write_vcf(ref,variants,sample):
    chrom_lengths = {}
    for l in open(ref+".fai"):
        row = l.strip().split()
        chrom_lengths[row[0]] = row[1]
    sys.stdout.write('##fileformat=VCFv4.2\n')
    sys.stdout.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">\n')
    sys.stdout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    sys.stdout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    for chrom in set([x["chrom"] for x in variants]):
        sys.stdout.write('##contig=<ID=%s,length=%s>\n' % (chrom,chrom_lengths[chrom]))
    sys.stdout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample)

    for v in variants:
        alts = [x["allele"] for x in v["alleles"] if x["allele"]!=v["ref"]]
        alt_ads = [str(x["ad"]) for x in v["alleles"] if x["allele"]!=v["ref"]]
        ordered_alleles = [v["ref"]] + alts
        for a in v["alleles"]:
            a["gt"] = ordered_alleles.index(a["allele"])
        v["alleles"] = sorted(v["alleles"],key=lambda x:x["gt"])
        if len(v["alleles"])==1:
            v["gt"] = f'{v["alleles"][0]["gt"]}/{v["alleles"][0]["gt"]}'
        elif len(v["alleles"])==2:
            v["gt"] = f'{v["alleles"][0]["gt"]}/{v["alleles"][1]["gt"]}'
        else:
            quit("Error!!!")
        v["alts"] = ",".join(alts) if len(alts)>0 else "."
        v["AD"] = ",".join([str(v["ref_ad"])] + alt_ads) if len(alts)>0 else str(v["ref_ad"])
        sys.stdout.write("%(chrom)s\t%(pos)s\t.\t%(ref)s\t%(alts)s\t.\t.\tDP=%(dp)s\tGT:AD:DP\t%(gt)s:%(AD)s:%(dp)s\n" % v)

def main(args):
    variant_positions = defaultdict(set)
    if not args.vcf_files and args.vcf_file_list:
        args.vcf_files = [x.strip() for x in open(args.vcf_file_list)]
    for f in args.vcf_files:
        for l in sp.Popen(f"bcftools filter -e 'QUAL<{args.min_variant_qual}' {f} | bcftools view -c 1 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ",shell=True,stdout=sp.PIPE).stdout:
            row = l.decode().strip().split()
            variant_positions[(row[0],row[1])].add((row[2],row[3]))

    variants = []
    for l in sp.Popen(f"htsbox pileup -f {args.ref} {args.bam}",shell=True,stdout=sp.PIPE).stdout:
        chrom,pos,ref,allele_str,info = row = l.decode().strip().split()
        ref = ref.upper()
        if args.vcf_files and (chrom,pos) not in variant_positions: continue
        alleles = allele_str.split(",")
        ads = [int(x) for x in info.split(":")[1].split(",")]
        filtered_alleles = []
        ref_ad = 0
        for all,ad in zip(alleles,ads):
            if all==ref: 
                ref_ad = ad
            if ad/sum(ads)>args.min_af:
                filtered_alleles.append({"allele":all,"ad":ad})
        filtered_alleles = sorted(filtered_alleles,key=lambda x: x["ad"],reverse=True)
        if not args.no_ploidy:
            filtered_alleles = filtered_alleles[:2]
        ref,filtered_alleles = recode_indels(ref,filtered_alleles)
        variants.append({
            "chrom":chrom,
            "pos":pos,
            "ref":ref,
            "alleles":filtered_alleles,
            "dp":sum(ads),
            "ref_ad": ref_ad
        })
        
    write_vcf(args.ref,variants,args.sample)
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',type=str,help='Sambamba coverage file')
parser.add_argument('--ref',type=str,help='Sambamba coverage file')
parser.add_argument('--sample',type=str,help='Sambamba coverage file')
parser.add_argument('--vcf-files',type=str,nargs="+",help='Sambamba coverage file')
parser.add_argument('--vcf-file-list',type=str,help='Sambamba coverage file')
parser.add_argument('--min-af',type=float,default=0.05,help='Sambamba coverage file')
parser.add_argument('--no-ploidy',action="store_true",help='Sambamba coverage file')
parser.add_argument('--min-variant-qual',default=30,type=int,help='Quality value to use as a minimum threshold')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)