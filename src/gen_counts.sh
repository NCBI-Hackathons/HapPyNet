#!/bin/bash

# run from /opt/SRR_done
# run as parallel -j12 './gen_counts.sh {}' ::: $(ls SRR* ERR* -d)

if [ $# -lt 1 ]; then
    exit 1
fi

if [ ! -d ${1} ]; then
    exit 1
fi

cd ${1}

# Keep site depth >=8 (also individual depth >=8, these are single-sample .vcf.gz files)
gatk SelectVariants --java-options "-Xmx4G" -select "DP>=8" -V ${1}.sort.markd.recal.vcf.gz -O ${1}.sort.markd.recal.DP8.vcf.gz

# Keep only dbSNP
gatk SelectVariants --java-options "-Xmx4G" -select "(vc.hasID())" -V ${1}.sort.markd.recal.DP8.vcf.gz -O ${1}.sort.markd.recal.DP8.dbsnp.vcf.gz

# Add filter for GQ20
gatk VariantFiltration --java-options "-Xmx4G" --genotype-filter-expression "GQ < 20" --genotype-filter-name "GQ20" -V ${1}.sort.markd.recal.DP8.dbsnp.vcf.gz -O ${1}.sort.markd.recal.DP8.dbsnp.GQ20.vcf.gz

# Apply GQ20 Filter
gatk SelectVariants --java-options "-Xmx4G" --set-filtered-gt-to-nocall --exclude-non-variants -V ${1}.sort.markd.recal.DP8.dbsnp.GQ20.vcf.gz -O ${1}.sort.markd.recal.DP8.dbsnp.GQ20.filt.vcf.gz

# Add filter for GQ90
gatk VariantFiltration --java-options "-Xmx4G" --genotype-filter-expression "GQ < 90" --genotype-filter-name "GQ90" -V ${1}.sort.markd.recal.DP8.dbsnp.vcf.gz -O ${1}.sort.markd.recal.DP8.dbsnp.GQ90.vcf.gz

# Apply Filter for GQ90 and keep site depth >=50
gatk SelectVariants --java-options "-Xmx4G" -select "DP>=50" --set-filtered-gt-to-nocall --exclude-non-variants -V ${1}.sort.markd.recal.DP8.dbsnp.GQ90.vcf.gz -O ${1}.sort.markd.recal.DP8.dbsnp.DP50.GQ90.filt.vcf.gz



# Generate counts for the various .vcf.gz generated above
cd ../
mkdir -p \
      counts \
      counts_DP8 \
      counts_DP8_dbsnp \
      counts_DP8_dbsnp_GQ20 \
      counts_DP8_dbsnp_GQ20_DP50_GQ90

bedtools intersect \
	 -a /opt/ldetect_GRCh38/EUR_ldetect.bed \
	 -b ${1}/${1}.sort.markd.recal.vcf.gz \
	 -c \
    | sort -k1,1V -k2,2n > \
	   counts/${1}.count

bedtools intersect \
	 -a /opt/ldetect_GRCh38/EUR_ldetect.bed \
	 -b ${1}/${1}.sort.markd.recal.DP8.vcf.gz \
	 -c \
    | sort -k1,1V -k2,2n \
	   > counts_DP8/${1}.count

bedtools intersect \
	 -a /opt/ldetect_GRCh38/EUR_ldetect.bed \
	 -b ${1}/${1}.sort.markd.recal.DP8.dbsnp.vcf.gz \
	 -c \
    | sort -k1,1V -k2,2n \
	   > counts_DP8_dbsnp/${1}.count

bedtools intersect \
	 -a /opt/ldetect_GRCh38/EUR_ldetect.bed \
	 -b ${1}/${1}.sort.markd.recal.DP8.dbsnp.GQ20.filt.vcf.gz \
	 -c \
    | sort -k1,1V -k2,2n \
	   > counts_DP8_dbsnp_GQ20/${1}.count

bedtools intersect \
	 -a /opt/ldetect_GRCh38/EUR_ldetect.bed \
	 -b ${1}/${1}.sort.markd.recal.DP8.dbsnp.DP50.GQ90.filt.vcf.gz \
	 -c \
    | sort -k1,1V -k2,2n \
	   > counts_DP8_dbsnp_GQ20_DP50_GQ90/${1}.count

exit 0