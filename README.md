# HapPyNet

Association of **Hap**lotype blocks to phenotypes using a neural **net**work machine learning method (in **Py**thon)

A tool to test the association of variants in haplotype blocks to phenotypes.
This tool takes variants (VCF format) called by any technolgy like Exome, WGS, RNASeq or SNP Genotyping Arrays and generates association test results

![alt text](docs/images/concept.png)

<p align="center">
<img src="./docs/images/pipeline_0.png" height="300">
</p>

## Usage

   * Step 1 : Generate SNP count matrix (Number of SNPs per LD block) \
     See README [here](docs/README.md)
   * Step 2 : Differential LD block SNP load between two classes (ex: disease vs. normal)

## Method
   * Call variants using any platform (RNASeq, Exome, Whole Genome or SNP Arrays)
   * Group variants by haplotype blocks to compute SNP load in each haplotype block
   * Classify samples into disease vs normal, based on SNP load(number of SNPs per LD block) using a TensorFlow classifier
   * Associate haplotypes with phenotype. As of Apr 2018, this is NOT implemented

![alt text](docs/images/flow.png)

## Data sources

   * LD BLocks : Non-overlapping LD blocks derived from 1KG data (hg19) were obtained from : *Approximately independent linkage disequilibrium blocks in human populations,Bioinformatics. 2016 Jan 15; 32(2): 283–285 [doi:  10.1093/bioinformatics/btv546]*. Using NCBI's online remapping tool these regions were mapped to GRCh38 with merge fragments turned ON to make sure each LD block is not fragmented

   * RNASeq samples: Initial training set from healthy and disease samples were obtained from SRA. The disease sample selection query was: `(AML) AND "Homo sapiens"[orgn:__txid9606] NOT ChIP-Seq`
   * The specific list of Sample IDs for AML and normal samples are in test_data

## RNASeq Variant Calling Pipeline

   * RNASeq sample reads were aligned using HiSat2
   * Variants were called using [GATK version 4.0.3.0](https://github.com/broadinstitute/gatk/releases/download/4.0.3.0/gatk-4.0.3.0.zip)
   
## Machine Learning

   * We trained a classifier using TensorFlow 
   * We are exploring standard differential gene expression methods from Bioconductor
![alt text](docs/images/eset.png)
