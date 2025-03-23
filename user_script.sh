
IFS=$'
	'
set -vxeuo pipefail
hostname
pwd
date
# Substitution will be similar to snakemake "shell".
python -m falcon_phase.mains.emit_haplotigs ../phasing/phased.txt ../prep-emit/BC.bed ../prep-emit/p_h_ctg.fasta bedtools pseudohap

rm tmp_phase1.txt
rm tmp_phase0.txt
rm ../prep-emit/p_h_ctg.fasta
rm ../prep-emit/p_h_ctg.fasta.fai

mv f0.fasta phased.0.fasta
mv f1.fasta phased.1.fasta
mv b0.bed phased.0.bed
mv b1.bed phased.1.bed

date
