#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=1:00:00
#PBS -A wff3_g_g_lc_icds-training

module use /gpfs/group/RISE/sw7/modules
module load bwa

cd ~/scratch/compbio
pwd

# index scaffold
bwa index ./assembly/scaffolds.fasta

# map reads to scaffold
bwa mem ./assembly/scaffolds.fasta ./trimmed/evol1_R1_fastq.gz ./trimmed/evol1_R2.fastq.gz > mappings/evol1.sam
bwa mem ./assembly/scaffolds.fasta ./trimmed/evol2_R1_fastq.gz ./trimmed/evol2_R2.fastq.gz > mappings/evol2.sam

