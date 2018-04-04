# HOWTO run generate_ld_counts.py

## Requirements

  All tools should be included in the ${PATH} variable. Support for custom tool location and auto-discovery will be added at a later date

  * [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

### If using SRA data as input:

  * [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
  * [samtools](http://www.htslib.org/download/)
  * [gatk](https://software.broadinstitute.org/gatk/download/)
  * [gatk resource bundle for GRCh38](https://software.broadinstitute.org/gatk/download/bundle)
    * `Homo_sapiens_assembly38.fasta.gz` must be unzipped
    * .hts files must be generated with `hisat2-build Homo_sapiens_assembly38.fasta`

## Usage
  * `generate_ld_counts.py` must be run from the `HapPyNet/src` directory, or things will break. We will fix this later, I promise

### Using SRA data as input

#### Run with one SRA sample
  * `python generate_ld_counts.py -s SRR#######`

#### Run with many SRA samples
  * `python generate_ld_counts.py -b input.list`
    * input.list should have 1 SRA ID (e.g. SRR1608593) per line

#### Run with an existing .vcf or .vcf.gz file
  * `python generate_ld_counts.py -v file.vcf`

#### Additional options
  * `-o <dirname>` set output directory (default ../results)
  * `-i <prefix>` add prefix to output files
  
## Expected Output
  * One directory per sample (e.g. `../results/SRR#######`) containing all files from alignments to .vcf.gz for that sample
  * One directory (`../results/counts`) with .count files containing variant counts per haplotype block, before filtering (useful for QC)
  
## Additional Steps
  * Run `filter_vcfs.sh` from the `../results/` directory to perform stricter filtering of VCFs, and to generate .count files for use with Step 2

