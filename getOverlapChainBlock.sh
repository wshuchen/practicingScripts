#!/bin/bash

# getOverlapChainBlock.sh
# by Wenshu Chen, 2023

# Examine a chain file for overlapping blocks.

# Example from hg38ToCanFam3.over.chain, the second block is part of the first block,
# which was mapped to a different region of canFam3 genome.
# chain 31418 chr1 248956422 + 101255884 101256672 chr22 61439934 - 52056696 52057169 570061
# chain 9425 chr1 248956422 + 101255991 101256245 chr29 41845238 - 34355552 34355768 720475

# Usage:
# ./getOverLapChainBlock.sh [chain_file]

chain_file=$1

# Get the query chromosome names.
zgrep chain ${chain_file} | cut -d" " -f3 | sort -u > query_chrs

# Could just look at those "regular" chromosomes.
# zgrep chain ${chain_file} | grep -vE "Un|random|alt|fix" |\
#     cut -d" " -f3 | sort -u > query_chrs

# Get the overlap chain blocks for each chromosome.
# The idea is that after sorting, if the first end - next start > 0, then
# there is an overlap.
for i in $(cat query_chrs); do
    echo "Working on ${i}"
    zgrep chain ${chain_file} | awk -v chr=${i} '$3 == chr' > ${i}_block

    # Get the query start and end positions, and
    # delete the first of the start and the last of the end.
    cut -d" " -f6,7 ${i}_block | tr " " "\t" | sort -k1,1n > ${i}_start_end
    paste <(cut -f1 ${i}_start_end | tail -n +2) \
          <(cut -f2 ${i}_start_end | head -n -1) |\
          sort -u > ${i}_start_end_merged

    # Get the overlap block positions.
    awk '$2 - $1 > 0' ${i}_start_end_merged > ${i}_overlap_pos

    # Get the overlap chain blocks.
    while read line; do
        pos=($line) # (start end)
        awk -v s=${pos[0]} -v e=${pos[1]} '$6 == s || $7 == e' ${i}_block \
            >> chr_temp_out
    done < ${i}_overlap_pos
done

# There are duplicate chain blocks in the output because of ovelapping between pairs.
# Can sort them out but keep them for easier to see the overlap.
out_file=$(basename $chain_file | cut -d"." -f1)_overlap_blocks
sort -k3 -k7,7n chr_temp_out > ${out_file}

n=$(cat ${out_file} | wc -l)
pair=$(( $n/2 ))
echo ""
echo "There are ${pair} overlap pairs in the blocks of ${chain_file}. "
echo "The result was saved as ${out_file}"
echo ""

# Clean up.
rm chr* query_chrs
