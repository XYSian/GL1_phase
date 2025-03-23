#!/bin/sh
#SBATCH -p com
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -J testname
#SBATCH -o testname.o
#SBATCH -e testname.e
#SBATCH --mem=1G
export PATH=~/mambaforge/bin/conda:$PATH
#gffread f153.gff -T -o  f153.gtf
iqtree -s 12.5.renamed.mafft -T 18 -m MFP -bb 1000 -bnni -redo -o  17
#mafft --auto --thread 18 all2.fasta >all2.mafft
