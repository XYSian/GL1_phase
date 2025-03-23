#!/bin/sh
#SBATCH -p com
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J testname
#SBATCH -o testname.o
#SBATCH -e testname.e
#SBATCH --mem=1G
export PATH=~/mambaforge/bin/conda:$PATH

#makeblastdb -in nac_pep.fasta  -dbtype prot
#blastp -query  nac_pep.fasta  -db  nac_pep.fasta  -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out nac.blast
duplicate_gene_classifier ./nac > nac.log
