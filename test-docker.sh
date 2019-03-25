#!/bin/bash
# Test of PseudogeneSearcher.

set -e

echo "Requirements: docker, samtools, wget, gunzip"

echo "Create test dir"
mkdir -p test && cd test

echo "Download and index the reference sequence"
if [ ! -f hg38.fa ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && \
    gunzip hg38.fa.gz;
fi
if [ ! -f hg38.fa.fai ]; then
    samtools faidx hg38.fa
fi

# Download the data
wget -nc https://zgm.wum.edu.pl/share/ifae7aeZ/chr13.bam && \
wget -nc https://zgm.wum.edu.pl/share/ifae7aeZ/chr13.bam.bai

# Run

docker run --rm \
    -v $PWD/chr13.bam:/input/file.bam \
    -v $PWD/chr13.bam.bai:/input/file.bai \
    -v $PWD/hg38.fa:/reference/reference.fa \
    -v $PWD/hg38.fa.fai:/reference/reference.fa.fai \
    -v $PWD:/output \
    zgmwum/psedogenesearcher \
    --only-soft-clipped \
    --threads 12 \
    --standard-genome hg38 \
    --bam /input/file.bam \
    --output-json /output/output.json \
    --output-vcf /output/output.vcf \
    --output /output/output.log \
    --reference /reference/reference.fa

# Check the positive result
cat output.vcf | grep PASS | grep SKA3 && echo "[  OK  ] SKA3 pseudogene found"
cat output.vcf | grep PASS | grep SKA3 || echo "[ FAIL ] SKA3 pseudogene NOT found"

