# HaploPhen

Association of haplotype blocks to phenotypes using machine learning methods

A tool to test the association of variants in haplotype blocks to phenotypes.
This tool takes variants called by any technolgy like Exome, WGS, RNASeq or SNPArrays in VCF format and generates association test results

![alt text](docs/images/concept.png)

<p align="center">
<img src="./docs/images/pipeline_0.png" height="300">
</p>

## Usage
   * Step 1 : Generate SNP count matrix (Number of SNPs per LD block) \
     See README [here](docs/README.md)
   * Step 2 : Run a neural net to classify samples into disease vs. normal
     jupyter notebook training.ipynb

## Method
   * Call variants using any platform (RNASeq, Exome, Whole Genome or SNP Arrays)
   * Group variants by haplotype blocks to compute SNP load in each haplotype block
   * Classify samples into disease vs normal, based on SNP load(number of SNPs per LD block) using a TensorFlow classifier
   * Associate haplotypes with phenotype. As of Apr 2018, this is NOT implemented

![alt text](docs/images/flow.png)

## Data sources

   * LD BLocks : Non-overlapping LD blocks derived from 1KG data (hg19) were obtained from : *Approximately independent linkage disequilibrium blocks in human populations,Bioinformatics. 2016 Jan 15; 32(2): 283–285 [doi:  10.1093/bioinformatics/btv546]*. Using NCBI's online remapping tool these regions were mapped to GRCh38 with merge fragments turned ON to make sure each LD block is not fragmented

   * RNASeq samples: Initial training set from healthy and disease samples were obtained from SRA. The disease sample selection query was: (AML) AND "Homo sapiens"[orgn:__txid9606] NOT ChIP-Seq. Here are the specific list of Sample IDs for AML (ref_data/SRR_50aml_small.list) and normal samples (ref_data/SRR_47normal_small.list)

## RNASeq Variant Calling Pipeline

   * RNASeq sample reads were aligned using HiSat2
   * Variants were called using GATK and quality filtered at read depth of 50 and genotype quality of 90
![alt text](docs/images/VariantsvsCoverageDP50_GQ90.png)

## Machine Learning

   * We trained a classifier with a 4 layer NeuralNet using TensorFlow with leave-one-out cross validation resulting in 99% accuracy of our model!

   ![alt text](docs/images/simple_neural_network_header.jpg)

   Here' the confusion matrix:
   * We are exploring standard differential gene expression methods from Bioconductor
![alt text](docs/images/eset.png)
