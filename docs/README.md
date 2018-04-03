# HOWTO run generate_ld_counts.py

## Requirements

  All tools should be included in the ${PATH} variable. Support for custom tool location and auto-discovery will be added at a later date

  * [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

### If using SRA data as input:

  * [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
  * [samtools](http://www.htslib.org/download/)
  * [gatk](https://software.broadinstitute.org/gatk/download/)

## Usage

### Using SRA data as input

#### Run with one SRA sample
  * `python generate_ld_counts.py -s SRR#######`

#### Run with many SRA samples
  * `python generate_ld_counts.py -b input.list`
  * input.list should have 1 SRA ID per line

#### Run with an existing .vcf or .vcf.gz file
  * `python generate_ld_counts.py -v file.vcf`
