# amplicon-seq

## Install 

Use these commands to install the scripts
```
git clone https://github.com/LSHTMPathogenSeqLab/amplicon-seq.git
cd amplicon-seq
python setup.py install
```

## Updating
First change directory to the repo on your computer and then run the following commands
```
cd /path/to/amplicon-seq
git pull
python setup.py install
```


## Usage
```
usage: amplicon_script.py [-h] --index-file INDEX_FILE --read1 READ1 --read2
                          READ2 --ref REF --gff GFF --bed BED [--trim]
                          [--trim-qv TRIM_QV] [--min-base-qual MIN_BASE_QUAL]

Amplicon sequencing analysis script

optional arguments:
  -h, --help            show this help message and exit
  --index-file INDEX_FILE
                        CSV file containing fields "Sample,I1,I2" (default:
                        None)
  --read1 READ1         Forward read file (default: None)
  --read2 READ2         Reverse read file (default: None)
  --ref REF             Reference fasta (default: None)
  --gff GFF             GFF file (default: None)
  --bed BED             BED file with genes/amplicon locations (default: None)
  --trim                Perform triming (default: False)
  --trim-qv TRIM_QV     Quality value to use in the sliding window analysis
                        (default: 20)
  --min-base-qual MIN_BASE_QUAL
                        Minimum base quality to use by freebayes (default: 30)
```

## Nanopore Sequencing - Mapping, Variant Calling & Annotation (Run AFTER Demultiplexing by nanopore and primer barcodes)
Place all demultiplexed fastq files into a directory and use these commands to generate a sample ID list in .csv format with 'sample' as the column header.
```
ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv
```
Compress all FASTQ files.
```
for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done
```

### Usage
```
nanopore_amplicon_script.py [-h] --index-file INDEX_FILE --ref REF --gff GFF --bed BED
                                             [--min-base-qual MIN_BASE_QUAL]
                                             [--version]

Nanopore amplicon sequencing analysis script

optional arguments:
  -h, --help            show this help message and exit
  --index-file INDEX_FILE
                        samples.csv with the "sample" column for sample IDs;
                        created by "demux_nanopore_plates.py" if using 96-well
                        plates (default: None)
  --ref REF             Reference fasta (default: None)
  --gff GFF             GFF file (default: None)
  --bed BED             BED file with genes/amplicon locations (default: None)
  --min-base-qual MIN_BASE_QUAL
                        Minimum base quality to use by freebayes (default: 30)
  --version             show program's version number and exit

```

## Human genotyping amplicon pipeline

Download the latest clinvar vcf to annotate variations with rsIDs and clinical significance (now required).

```
curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz > clinvar_GRCh38.vcf.gz
tabix -f clinvar_GRCh38.vcf.gz

```

### Usage
```
human_amplicon_script.py [-h] --index-file INDEX_FILE --read1 READ1 --read2 READ2 
                              --ref REF --gff GFF --bed BED --clinVar CLINVAR 
                                [--trim]
                                [--trim-qv TRIM_QV]
                                [--min-base-qual MIN_BASE_QUAL]
                                [--min-adf MIN_ADF]
                                [--vcf-files VCF_FILES [VCF_FILES ...]]
                                [--min-variant-qual MIN_VARIANT_QUAL]
                                [--min-sample-af MIN_SAMPLE_AF]
                                [--per-sample-only] [--version]

Amplicon sequencing analysis script

optional arguments:
  -h, --help            show this help message and exit
  --index-file          INDEX_FILE
                        CSV file containing fields "Sample,I1,I2" (default:
                        None)
  --read1 READ1         Forward read file (default: None)
  --read2 READ2         Reverse read file (default: None)
  --ref REF             Reference fasta (default: None)
  --gff GFF             GFF file (default: None)
  --bed BED             BED file with genes/amplicon locations (default: None)
  --clinVar CLINVAR     ClinVar SNP annotation file (default: None)
                        Position info file (default: None)
  --trim                Perform triming (default: False)
  --trim-qv TRIM_QV     Quality value to use in the sliding window analysis
                        (default: 20)
  --min-base-qual       MIN_BASE_QUAL
                        Minimum base quality to use by freebayes (default: 30)
  --min-adf MIN_ADF     Set a minimum frequency for a mixed call (default:
                        None)
  --vcf-files           VCF_FILES [VCF_FILES ...]
                        VCF files with positions to include (optional)
                        (default: None)
  --min-variant-qual    MIN_VARIANT_QUAL
                        Quality value to use in the sliding window analysis
                        (default: 30)
  --min-sample-af       MIN_SAMPLE_AF
                        Quality value to use in the sliding window analysis
                        (default: 0.05)
  --per-sample-only     Perform triming (default: False)
  --version             show program's version number and exit


```

## Plot read depth for variants of interest
Create box plots using the plot_read_depth_human.r or plot_read_depth_Pfalciparum.r in the plots folder to show the read depth at positions of interest using the depth_info.txt output file from the amplicon sequencing pipeline.

Example using loci in humans associated with severe malarial disease

![alt text](https://github.com/LSHTMPathogenSeqLab/amplicon-seq/blob/main/plots/read_depth_plot.png)

