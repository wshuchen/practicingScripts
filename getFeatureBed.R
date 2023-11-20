#!/usr/bin/env Rscript

# getFeatureBed.R 
# by Wenshu Chen, 2023

# Process a GAinSAW annotation CSV of lifted points
# into a BED with a selected feature.
# Valid features: pseudogene, five_prime_UTR, three_prime_UTR, CDS, intron,
# RNA, intergenic, multiple, all; "multiple" means more than one featuere.
# If feature == "all", generate a BED without selection.

# The CSV is tab-delimited.
# header and a line (broken up into short lines):
# qkey  qftype  tkey  tftype    
# qgann_chr  qgann_from  qgann_to  qgann_strand  qgann_ftype  qgann_gene  
# tgann_chr  tgann_from  tgann_to  tgann_strand  tgann_ftype  tgann_gene
# hg38_chr3_185802277+_1   intron  canfam3_chr34_18502545+  intron
# chr3  185698347  185821011  -  intron  IGF2BP2
# chr34 18410487   18522095   -  intron  IGF2BP2

# The output (one line):
# chr3  185802277  185802278  +  intron  IGF2BP2  hg38
# chr34 18502545   18502546   +  intron  IGF2BP2  canfam3

# Usage:
# ./getFeature.R [ lifted point CSV ] [ feature ]

args <- commandArgs(trailingOnly=TRUE)
annotation_csv  = args[1]
feature  = args[2]

library(stringr)

getFeatureBed = function(annotation_csv, feature)
    {
    df = read.csv(annotation_csv, sep = "\t")
    df = df[c("qkey", "qftype", "qgann_gene", "tkey", "tftype", "tgann_gene")]
    
    # Expand the keys to coordinate columns and then select and rearrange them.
    # qkey: hg38_chr3_185802277+_1
    # tkey: canfam3_chr34_18502545+
    # Note: pos is the start position in a 0-based format.
    df[c("query", "qchr", "qpos", "qstrand")] = str_split_fixed(df$qkey, "_", 4)
    df$qstart = sub("\\+", "", df$qpos)
    df$qend = as.integer(df$qstart) + 1
    df$qstrand = c("+")
    
    df[c("target", "tchr", "tpos")] = str_split_fixed(df$tkey, "_", 3)
    df$tstart = str_sub(df$tpos, 1, -2)
    df$tend = as.integer(df$tstart) + 1
    df$tstrand = str_sub(df$tpos, -1, -1)

    df = df[c("qchr", "qstart", "qend", "qstrand", 
            "qftype", "qgann_gene", "query", 
            "tchr", "tstart", "tend", "tstrand", 
            "tftype", "tgann_gene", "target")]
    df = df[order(df$qchr, df$qstart), ]

    if (feature == "all") {
        return(df)
        }
    else if (feature == "multiple") {
        return(df[grepl(",", df$qftype), ])
        }
    else {
        # Get the subset with single feature type and then feature subset.
        df = df[!grepl(",", df$qftype), ]
        feature_df = df[df$qftype == feature, ]
        return(feature_df)
        } 
    }

feature_df = getFeatureBed(annotation_csv, feature)
out_file = paste0(sub(".csv", "_", annotation_csv), feature, ".bed")
write.table(feature_df, out_file, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = "\t")
