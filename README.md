Some scripts for projects

_summarizeFeatures.py_ and _getFeatureBed.R_ are for tab-delimited data (data frame or CSV) generated with [GAinSAW](https://github.com/BrendelGroup/GAinSAW) package. The Python script summarizes the count and percentage of lifted and unlifted points according to feature types; the R script turns a data frame into a simplified BED file with all entries or a selected feature type. 

_getChainGenomeCoverage.sh_ and _getOverlapChainBlock.sh_ are two shell scripts with functions of what the names say for UCSC chain file.

_sra2snp.sh_ is a workflow of SNP calling from paired-end reads. It takes a SRA run number and a reference genome, going through read download to VCF output.

Please see the scripts for detailed usages.
