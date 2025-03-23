#!/bin/sh
#SBATCH -p com
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -J testname
#SBATCH -o testname.o
#SBATCH -e testname.e
#SBATCH --mem=1G
export PATH=~/mambaforge/bin/conda:$PATH
#gffread f153.gff -T -o  f153.gtf
mafft --auto --thread 12 atpH.bro.fasta > aligned_atph.fasta

