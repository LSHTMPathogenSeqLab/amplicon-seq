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