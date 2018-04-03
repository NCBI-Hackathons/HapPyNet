#!/bin/bash

if [ $# -lt 1 ]; then
    exit 1
fi

ALIGN_THREADS=12
SORT_THREADS=3

SRR="${1}"
mkdir -p ${SRR}
cd ${SRR}

# Align
TMPDIR=$(mktemp -d -p $(pwd))
hisat2 -p ${ALIGN_THREADS} -x /opt/grch38/Homo_sapiens_assembly38 --sra-acc ${SRR} --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA| samtools sort -@ ${SORT_THREADS} > ${SRR}.sort.bam -T ${TMPDIR}
rmdir ${TMPDIR}

# MarkDuplicates
/opt/gatk-4.0.3.0/gatk MarkDuplicates --INPUT ${SRR}.sort.bam --OUTPUT ${SRR}.sort.markd.bam --METRICS_FILE ${SRR}.sort.markd.metrics.bam

# BaseRecalibrator
/opt/gatk-4.0.3.0/gatk BaseRecalibrator -R /opt/grch38/Homo_sapiens_assembly38.fasta --known-sites /opt/grch38/dbsnp_138.hg38.vcf.gz --known-sites /opt/grch38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -I ${SRR}.sort.markd.bam -O ${SRR}.sort.markd.recal.table

# PrintReads
/opt/gatk-4.0.3.0/gatk ApplyBQSR -R /opt/grch38/Homo_sapiens_assembly38.fasta -I ${SRR}.sort.markd.bam --bqsr-recal-file ${SRR}.sort.markd.recal.table -O ${SRR}.sort.markd.recal.bam

# HaplotypeCaller
/opt/gatk-4.0.3.0/gatk HaplotypeCaller -R /opt/grch38/Homo_sapiens_assembly38.fasta --dbsnp /opt/grch38/dbsnp_138.hg38.vcf.gz -I ${SRR}.sort.markd.recal.bam -O ${SRR}.sort.markd.recal.vcf.gz

# Bedtools intersect
/home/ubuntu/bin/bedtools intersect -a /opt/ldetect_GRCh38/EUR_ldetect.bed -b ${SRR}.sort.markd.recal.vcf.gz -c | sort -k1,1V -k2,2n > ${SRR}.count

