# HapPyNet

## Test Data

This folder contains information about the test data set used to validate the HapPyNet program

## SRR_lists
RNAseq data (from [SRA](https://www.ncbi.nlm.nih.gov/sra)) for [50 AML](SRR_lists/SRR_aml_candidate50.list) and [47 normal](SRR_lists/SRR_normal_candidate47.list) samples were considered for use when developing the method. Resource availability only allowed the complete pipeline to br run on [47 AML](SRR_lists/SRR_aml_test47.list) and [30 normal](SRR_lists/SRR_normal_test30.list) samples during the NCBI Hackathon,
  
## SRR_runtables
Complete runtable (from [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) containing metadata for all 97 candidate samples considered for use. The relationship of total variant load per sample was compared to the MBases of sequence data available, and this is an important QC step.

## test_counts
Count files generated from strictly filtered VCFs for all 77 samples analyzed are provided

## test_vcfs
Examples VCFs are provided for 10 samples (5 AML, 5 normal) for use in validating the upstream variant-calling 
