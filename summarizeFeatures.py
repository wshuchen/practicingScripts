#!/usr/bin/env python

# summarizeFeatures.py
# by Wenshu Chen, 2023

# Get the count and percentage of lifted and unlifted points
# in each genome feature type of a GAinSAW CSV.

# The CSV is tab-delimited.
# header and a line (broken up into short lines):
# qkey  qftype  tkey  tftype    
# qgann_chr  qgann_from  qgann_to  qgann_strand  qgann_ftype  qgann_gene  
# tgann_chr  tgann_from  tgann_to  tgann_strand  tgann_ftype  tgann_gene
# hg38_chr3_185802277+_1   intron  canfam3_chr34_18502545+  intron
# chr3  185698347  185821011  -  intron  IGF2BP2
# chr34 18410487   18522095   -  intron  IGF2BP2

# Output of one sample:
# feature,lifted,lifted_pct,unlifted,unlifted_pct, total,total_pct,
# lifted_feature_pct,unlifted_feature_pct
# pseudogene,1398,1.39,1181,1.19,2579,1.29,54.21,45.79
# CDS,1951,1.94,64,0.06,2015,1.01,96.82,3.18
# 5'-UTR,175,0.17,65,0.07,240,0.12,72.92,27.08
# 3'-UTR,1463,1.45,463,0.47,1926,0.96,75.96,24.04
# RNA,12360,12.29,12165,12.23,24525,12.26,50.4,49.6
# intron,40678,40.45,32522,32.71,73200,36.6,55.57,44.43
# intergenic,36949,36.74,49825,50.11,86774,43.39,42.58,57.42
# multiple,5591,5.56,3150,3.17,8741,4.37,63.96,36.04

# If Excel:
# Copy and paste into one cell and use Data > Text to Columns to expand.

# Usage: ./summarizeFeatures.py [ lifted_point_csv ] [ unlifted_point_csv ]

import os
import sys
import pandas as pd

def count_feature_points(annotation_csv):
    """
    Count the points with a feature in annotation CSVs.
    Return a data frame with count and percentage.
    "multiple" means more than one feature.
    """
    df = pd.read_csv(annotation_csv, sep="\t", dtype=object)
    # Get the point numbers in feature.
    if "qftype" in df.columns:  # lifted points
        print("\nLifted points:", len(df))
        feature_count = df.groupby("qftype").describe()["qkey"].drop(["top", "freq"], 
                        axis=1).sort_values(by="count").drop(["unique"], axis=1)
    elif "uftype" in df.columns:   # unlifted points
        print("\nUnLifted points:", len(df))
        feature_count = df.groupby("uftype").describe()["ukey"].drop(["top", "freq"], 
                        axis=1).sort_values(by="count").drop(["unique"], axis=1)
    single_feature_count = feature_count[~feature_count.index.str.contains(",")]
    multi_feature_count = feature_count[feature_count.index.str.contains(",")]
    multi_feature_df = pd.DataFrame({"multiple": sum(multi_feature_count["count"])}, 
                        index=["count"]).transpose()
    feature_count_df = pd.concat([single_feature_count, multi_feature_df])
    # Calculate the percentage.
    feature_count_df["percentage"] = [i for i in feature_count_df["count"]
                                    /sum(feature_count_df["count"]) * 100]
    feature_count_df = feature_count_df.reset_index()
    feature_count_df.columns = ["feature", "count", "percentage"]
    feature_count_df.replace("five_prime_UTR", "5'-UTR", inplace=True)
    feature_count_df.replace("three_prime_UTR", "3'-UTR", inplace=True)
    return feature_count_df

def summarize_features(lifted_point_csv, unlifted_point_csv):
    """
    Summarize the points with a feature in both lifted and unlifted ones.
    Calculate percentage for both and total.
    Save the result as CSV.
    """
    lifted_df = count_feature_points(lifted_point_csv).sort_values("feature")
    unlifted_df = count_feature_points(unlifted_point_csv).sort_values("feature")
    # how = "outer", in case of any missing featuers in either df.
    # NA to 0 for calculation.
    sum_df = pd.merge(lifted_df, unlifted_df, how = "outer", on = "feature").fillna(0)
    sum_df.columns = ["feature", "lifted", "lifted_pct", "unlifted", "unlifted_pct"]
    sum_df["total"] = sum_df.lifted + sum_df.unlifted
    sum_df["total_pct"] = sum_df.total/sum(sum_df.total) * 100
    # Round the numbers. Note this may result in the percentages being 
    # slightly off 100 in some columns.
    for i in ["lifted_pct", "unlifted_pct", "total_pct"]:
        sum_df[i] = [round(j, 2) for j in sum_df[i]]
    feature_list = ["pseudogene", "CDS", "5'-UTR", "3'-UTR", 
                    "RNA", "intron", "intergenic", "multiple"]
    sum_df["feature"] = pd.Categorical(sum_df["feature"], feature_list)
    sum_df = sum_df.sort_values("feature").reset_index(drop = True)
    sum_df["lifted_feature_pct"] = pd.DataFrame([round(i, 2) for i in 
                                    sum_df["lifted"]/sum_df["total"] * 100])
    sum_df["unlifted_feature_pct"] = pd.DataFrame([round(i, 2) for i in 
                                    sum_df["unlifted"]/sum_df["total"] * 100])
    print("")
    print(sum_df.to_string(index = False))
    print("")
    out_name = os.path.basename(lifted_point_csv).rsplit("_", 1)[0] + "_summary.csv"
    sum_df.to_csv(out_name, index = False)
    print(f"The result was saved as {out_name}\n")


lifted_point_csv = sys.argv[1]      # *_qt.csv
unlifted_point_csv = sys.argv[2]    # *_qu.csv

summarize_features(lifted_point_csv, unlifted_point_csv)
