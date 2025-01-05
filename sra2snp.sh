#!/bin/bash

# sra2snp.sh
# by Wenshu Chen, 2022

# This script is a workflow of SNP VCF file preparation from paired-end reads:
# SRA read retrival, fastq quality check, adapter removal, 
# alignment to reference genome, duplicate removal, 
# variant calling, and initial filtering for quality SNP.

# Run with a singularity container with all the programs, 
# or have all the programs in the PATH.
# For multiple samples, download SRA run list
# and run "while read" loop with the list on command line for all samples.

# GENOME: reference genome fasta, bowtie2 indexed, otherwise run bowtie2-build
# to index the reference sequence.

# Test reference genome: Apis mellifera (honey bee), Amel_HAv3.1, GCF_003254395.2,
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003254395.2/
# Test sample: Africanized honey bee whole genome sequencing, SRR3654792, 
# https://www.ncbi.nlm.nih.gov/sra/?term=SRR3654792

# Usage: ./sra2snp.sh [ GENOME ] [ SRR number ] [ minDP ] [ threads ]

# Output: a VCF file and a report of SNP count.

GENOME=$1   # Reference genome
SRRn=$2     # SRR number
minDP=$3    # minimum coverage
threads=$4

# Retrieve SRA run using SRA Toolkit and run number.
prefetch ${SRRn}   	# Create a directory named ${SRRn} 
cd ${SRRn}			

# Convert the run to fastq.
fasterq-dump ${SRRn}.sra	# ${SRRn}_1.fastq, ${SRRn}_1.fastq
gzip *.fastq		    	# ${SRRn}_1.fastq.gz, ${SRRn}_1.fastq.gz

# Check read quality.
mkdir FastQC
fastqc ${SRRn}_1.fastq.gz ${SRRn}_2.fastq.gz -o FastQC

# Trim adapters and run quality check again. Output read file now has 
# "_val_1" or "_val_2" added before ".fq.gz", changed from ".fastq.gz".
trim_galore --paired --fastqc -j 4 -o FastQC \
	${SRRn}_1.fastq.gz ${SRRn}_2.fastq.gz
mv FastQC/*.fq.gz ./

# Align reads to reference genome using bowtie2. 
# --score-min is set according to BWASP (https://github.com/BrendelGroup/BWASP).
echo -e "\nMapping reads to reference genome.\n"

bowtie2 --score-min L,0,-0.6 \
	-x ${GENOME} \
	-1 ${SRRn}_1_val_1.fq.gz \
	-2 ${SRRn}_2_val_2.fq.gz \
	-p ${threads} \
	-S ${SRRn}.sam \
	2> ${SRRn}_align_log
samtools view -b ${SRRn}.sam |\
samtools sort -@ ${threads} -o ${SRRn}.bam
rm *.sam

# Remove duplicates and index the file.
java -jar picard.jar MarkDuplicates \
	-I ${SRRn}.bam \
	-O ${SRRn}_markdup.bam \
	--REMOVE_DUPLICATES \
	--ASSUME_SORTED true \
	-M ${SRRn}_metrics_markdup \
	
samtools index -@ ${threads} ${SRRn}_markdup.bam

echo -e "\nDone with processing ${SRRn} reads to bam.\n"

# Call variants with bcftools mpileup and call commands.
# Some parameters (BQ, MQ, and P) were set according to
# biscuit (BISulfite-seq CUI Toolkit, https://github.com/zhou-lab/biscuit)
# pileup for consistence with our other datasets.
echo -e "\nGenerating VCF file. \n"

bcftools mpileup \
	--min-BQ 20	--min-MQ 40 \
	--annotate FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP \
	--fasta-ref ${GENOME} \
	${SRRn}_markdup.bam |\
bcftools call \
	-mv -P 1e-3 \
	-f FORMAT/GQ \
	--threads ${threads} \
	-Oz -o ${SRRn}.vcf.gz

bcftools index ${SRRn}.vcf.gz

# Filter the VCF file.
# FILTER = "." = PASS.
echo -e "\nFiltering the VCF. \n"

# Calculate maxDP. Set according to 
# Heng Li, https://academic.oup.com/bioinformatics/article/30/20/2843/2422145
bcftools query -f '%INFO/DP\n' ${SRRn}.vcf.gz |\
awk '{DP += $1; meanDP = DP/NR; maxDP = meanDP + 4*sqrt(meanDP)}
	 END {printf "%.0f", maxDP}' \
	 > maxDP
maxDP=$(cat maxDP)

bcftools filter \
         -e "TYPE != 'snp' | ALT = 'N' | QUAL < 30 | FILTER != '.' | GT = '0/0' | \
        	MIN(GQ) < 20 | MIN(INFO/DP) < ${minDP} | MAX(INFO/DP) > ${maxDP}" \
         ${SRRn}.vcf.gz --threads ${threads} \
		 -Oz -o ${SRRn}_filtered_DP10_GT0111.vcf.gz
                
bcftools index ${SRRn}_filtered_DP10_GT0111.vcf.gz
rm maxDP

# Count the sites.
echo -e "${SRRn}.vcf.gz SNP sites (DP >= 10): \n" >> ${SRRn}_vcf_site_count
bcftools view -H -i 'INFO/DP >= 10' ${SRRn}.vcf.gz | wc -l >> ${SRRn}_vcf_site_count

echo -e "\n${SRRn}_filtered_DP10_GT0111.vcf.gz SNP sites (DP >= 10): \n" >> ${SRRn}_vcf_site_count
bcftools view -H ${SRRn}_filtered_DP10_GT0111.vcf.gz | wc -l >> ${SRRn}_vcf_site_count

cd ..
echo -e "\nDone with SRA read to SNP processing. \n"
