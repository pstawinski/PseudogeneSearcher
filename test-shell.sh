#!/bin/bash
# Test of PseudogeneSearcher.

set -e

echo "Requirements: java, maven, samtools, wget, gunzip"

echo "Compile"
mvn clean compile package assembly:single

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
java -jar ../target/PseudogeneSearcher-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
    --reference hg38.fa \
    --bam chr13.bam \
    --only-soft-clipped \
    --threads 12 \
    --output output.log \
    --output-json output.json \
    --output-vcf output.vcf \
    --standard-genome hg38

# Check the positive result
cat output.vcf | grep PASS | grep SKA3 && echo "[ OK ] SKA3 pseudogene found"
cat output.vcf | grep PASS | grep SKA3 || echo "[FAIL] SKA3 pseudogene NOT found"

