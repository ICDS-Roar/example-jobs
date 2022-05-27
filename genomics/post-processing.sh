#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=1:00:00
#PBS -A wff3_g_g_lc_icds-training

module use /gpfs/group/RISE/sw7/modules
module load r/3.5.2

cd /gpfs/scratch/cxb5898/compbio/
pwd
Rscript coverage-plot.R
