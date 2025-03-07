#!/bin/bash

# getChainGenomeCoverage.sh
# by Wenshu Chen, 2023

# Calculate the chain coverage for each chromosome and chromosome genome.
# Chain file may not include all the chromosomes, therefore,
# chromosome size files are neccesary for the genome coverage.
# Assume we will be using "genome".chrom.sizes file from UCSC.
# chr1	248956422

# Usage:
# ./getChainGenomeCoverage.sh [chain_file] [query_chr_file] [target_chr_file]

chain_file=$1
query_chr_file=$2
target_chr_file=$3

# chain lines, delimited with single space:
# chain 737886269 chr1 248956422 + 85019147 120505888 chr4 130910915 - 46034 30438425 21
# chain 618644190 chr9 138394717 + 51343 27329259 chr1 274330532 - 52672422 78464505 28

# Get the chromosome chain lines (headers) and chromosome names. 
# Exclude contigs:
# Unplaced("Un"), unlocalized("random"), haplotype("alt"), fix("fix").
# About that nice and messy "_": 
# It wasn't always there for those "extra" ones, and,
# for example, CavPor3 has chromosome names "scaffold_*".
zgrep chain ${chain_file} | grep -vE "Un|random|alt|fix" > chr_chain
cut -d" " -f3,4 chr_chain | tr " " "\t" | sort -u > query_chr_size
cut -d" " -f8,9 chr_chain | tr " " "\t" | sort -u > target_chr_size

# Calculate the coverage.
# Query
while read line; do
    chr_size=($line)
    chr=${chr_size[0]}
    size=${chr_size[1]}

    awk -v c=${chr} '$3 == c' chr_chain |\
        cut -d" " -f3,6,7 | tr " " "\t" | sort -k2,2n \
        > ${chr}_start_end

    # Merge overlapping coordinates before adding up
    bedtools merge -i ${chr}_start_end |\
        awk -v c=${chr} -v z=${size} \
            '{l += $3-$2 } END {printf "%s\t%d\t%d\t%.3f\n", c, z, l, l/z*100}' \
            >> query_chr_coverage
done < query_chr_size

# Target
while read line; do
    chr_size=($line)
    chr=${chr_size[0]}
    size=${chr_size[1]}

    awk -v c=${chr} '$8 == c' chr_chain |\
        cut -d" " -f8,11,12 | tr " " "\t" | sort -k2,2n \
        > ${chr}_start_end

    bedtools merge -i ${chr}_start_end |\
        awk -v c=${chr} -v z=${size} \
            '{l += $3-$2 } END {printf "%s\t%d\t%d\t%.3f\n", c, z, l, l/z*100}' \
            >> target_chr_coverage
done < target_chr_size

echo ""
for i in query target; do
    echo -n "Chain coverage of ${i} genome (%): "
    if [[ ${i} == "query" ]]; then
        chr_file=${query_chr_file}
    else
        chr_file=${target_chr_file}
    fi
    genome_size=$(grep -vE "Un|random|alt|fix" ${chr_file} |\
                  awk '{g += $2} END {print g}')
    chain_size=$(awk '{c += $3} END {print c}' ${i}_chr_coverage)
    awk -v gs=${genome_size} -v cs=${chain_size} 'BEGIN {printf "%.3f\n", cs/gs*100}'
done
echo ""

chain_name=$(basename ${chain_file} | cut -d"." -f1)
mv query_chr_coverage ${chain_name}_query_chr_coverage
mv target_chr_coverage ${chain_name}_target_chr_coverage

# Clean up.
rm chr* scaffold* *size
