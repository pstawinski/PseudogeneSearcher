# PseudogeneSearcher

## Purpose
This application identifies pseudogenes in Whole-Exome and Whole-Genome Next Generation Sequencing data.

## Usage
### Using docker

```
docker run --rm \
    -v file.bam:/input/file.bam \
    -v file.bai:/input/file.bai \
    -v hg38.fa:/reference/reference.fa \
    -v /path:/output \
    zgmwum/psedogenesearcher \
    --only-soft-clipped \
    --threads 12 \
    --standard-genome hg38 \
    --bam /input/file.bam \
    --output-json /output/out.json \
    --output-vcf /output/out.vcf \
    --reference /reference/reference.fa
```

### Using java
Standard usage may look like this
```
java -jar PseudogeneSearcher-0.0.1-SNAPSHOT.jar \
    --reference /reference/hg38.fasta \
    --bam /input/file.bam \
    --only-soft-clipped \
    --threads 12 \
    --output /output/1.output \
    --output-json /tmp/1.json \
    --output-vcf /tmp/1.vcf \
    --standard-genome hg38
```

### Program arguments
```TODO: will be published by the 2019.03.25```

## Output
Program will create both .vcf and .json file. However, I encourage to use the .json file and to make postprocessing of the file to remove
false positives according to adjust the precision and recall due to specific needs.

In our analysis we are using Apache Spark for postprocessing.

```
spark.read.json("input.json")
    .filter($"pass"===true) // take only 'passed' findings
    .filter($"fracIntronesCoveredIS">0.8) // take only these retropseudogenes where fraction of intrones with support of removal by the insert size is greater than 80%
    .filter($"fracIntronesCoveredClip">0.8) // take only these retropseudogenes where fraction of intrones with support of removal by the soft clipped reads aligning to the exon-exon junction is greater than 80%
    .withColumn("introns", size($"intronDetails"))
    .filter($"introns" > 2) // take only these retropseudogenes that have more than two intrones
    .write.json("filtered.json")
```

If you have a population of samples you can do population analysis to find unique retropseudogenes or samples with unusual amount of restropseudogenes by executing the code below. To raise the precision of population analysis the expectation of number of reads being higher that 2 for each intron was added. Please adjust it according to your specific needs.

```
spark.read.json("/aDirectoryWithJsons/")
    .filter($"pass"===true)
    .filter($"fracIntronesCoveredIS">0.8)
    .filter($"fracIntronesCoveredClip">0.8)
    .withColumn("introns", size($"intronDetails"))
    .filter($"introns" > 2)
    .withColumn("introne", explode($"intronDetails"))
    .filter($"introne.readsSupportingRemovalByIS" > 2)
    .filter($"introne.readsSupportingRemovalByAlignment" > 2)
    .groupBy($"sampleName",$"motherTranscript.gene" as "gene")
    .agg(collect_set($"pass"))
    .groupBy("sampleName")
    .agg(collect_set($"gene") as "gene")
    .sort(size($"gene"))
    .coalesce(1)
    .write.json("sampleToPseudogenes.json")
```


```
spark.read.json("/aDirectoryWithJsons/")
    .filter($"pass"===true)
    .filter($"fracIntronesCoveredIS">0.8)
    .filter($"fracIntronesCoveredClip">0.8)
    .withColumn("introns", size($"intronDetails"))
    .filter($"introns" > 2)
    .withColumn("introne", explode($"intronDetails"))
    .filter($"introne.readsSupportingRemovalByIS" > 2)
    .filter($"introne.readsSupportingRemovalByAlignment" > 2)
    .groupBy($"sampleName",$"motherTranscript.gene" as "gene")
    .agg(collect_set($"pass"))
    .groupBy("gene")
    .agg(collect_set($"sampleName") as "sampleName")
    .sort(size($"sampleName"))
    .coalesce(1)
    .write.json("pseudogeneToSamples.json")
```

## Testing
Example data and testing procedure will be published on 2019.03.25.


## Known issues
If your sample is so big, that you get into memory issues:

Recommended solution is to use the ```--position``` argument to process only a part of the huge input bam in a sigle run. For example: ```--position chr2:1-242194529``` to test ```chr2```.


## Development

### Compile and package
To reate an executable ```jar``` run:
```
mvn clean compile package
```
The ```jar``` file will be created in the ```target``` directory.

### Dockerize
To create and push a docker image run:
```
mvn compile jib:build
```
