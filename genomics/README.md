# Genomics Read Mapping Example

Adapted from [Genomics Tutorial](https://genomics.sschmeier.com/introduction.html) 
by Carrie Brown - carrie.brown@psu.edu

This tutorial will walk through mapping sample reads to an assembled reference genome, 
perform post-processing, and generate coverage statistics and plot using BWA, samtools 
and R on the ROAR cluster.

## Download sample data

```
wget -O trimmed.tar.gz https://osf.io/m3wpr/download
tar xvzf trimmed.tar.gz

wget -O assembly.tar.gz https://osf.io/t2zpm/download
tar xvzf assembly.tar.gz
```

## Map sample reads to reference genome using BWA

```
# load RISE software stack
module use /gpfs/group/RISE/sw7/modules

# search for BWA software
module spider bwa

# load BWA version 0.7.15
module load bwa/0.7.15

# view BWA help dialog
bwa

# view bwa index help dialog 
bwa index

# index scaffold
bwa index ./assembly/scaffolds.fasta

# create output directory
mkdir mappings

# view bwa mem help dialog
bwa mem

# map sample reads to scaffolded reference
bwa mem ./assembly/scaffolds.fasta ./trimmed/evol1_R1.fastq.gz ./trimmed/evol1_R2.fastq.gz > mappings/evol1.sam
bwa mem ./assembly/scaffolds.fasta ./trimmed/evol2_R1.fastq.gz ./trimmed/evol2_R2.fastq.gz > mappings/evol2.sam
```

## Post Processing:

```
# load samtools
module load samtools/1.10

# sort sam files to be name-sorted then convert to bam
samtools sort -n -O sam mappings/evol1.sam | samtools fixmate -m -O bam - mappings/evol1.fixmate.bam

# sort again into coordinate-order
samtools sort -O bam -o mappings/evol1.sorted.bam mappings/evol1.fixmate.bam

# remove duplicate reads (not always desired depending on experimental outcomes)
samtools markdup -r -S mappings/evol1.sorted.bam mappings/evol1.sorted.dedup.bam

# generate mapping statistics with samtools
samtools flagstat mappings/evol1.sorted.dedup.bam
samtools depth mappings/evol1.sorted.dedup.bam | gzip > mappings/evol1.depth.txt.gz

# extract results for contig NODE20 for R analysis
zcat mappings/evol1.depth.txt.gz | egrep '^NODE_20_' | gzip >  mappings/NODE_20.depth.txt.gz

# load R
module load r

# generate coverage plot
Rscript coverage-plot.R
```
