#!/bin/bash
#SBATCH -A b2011088
#SBATCH -p core -n 4
#SBATCH -t 3:00:00
#SBATCH -J sample_hmmscan
#SBATCH --mail-user christian.pou@ki.se
#SBATCH --mail-type ALL
#SBATCH -e log/sample_hmmscan.err
#$BATCH -o log/sample_hmmscan.out 

set -euo pipefail

module load hmmer
module load emboss

sample_name=$(basename assembly/*_allctgs.fa _allctgs.fa)
threads=16
#First param is number of procs
if [ -z $0 ]; then threads=$0; fi

echo "Processing ${sample_name} with ${threads} procs"

mkdir -p tax_assign/hmmscan/input

if [ ! -e tax_assign/hmmscan/input/${sample_name}_diamond_nohits.fa ];
then
cat <(python scripts/extract_diamond_nohits.py assembly/*_allctgs.fa tax_assign/diamond/*_allctgs_diamond_nr.sam.gz) \
	<(python scripts/extract_diamond_nohits.py assembly/*_se.fa tax_assign/diamond/*_se_diamond_nr.sam.gz)  | seqtk seq -L 400 > tax_assign/hmmscan/input/${sample_name}_diamond_nohits.fa
fi

cd tax_assign/hmmscan
make -rf ../../steps/tax_hmmscan.mak SAMPLE=${sample_name} threads=${threads} UPPMAX=1 in_fasta=input/${sample_name}_diamond_nohits.fa
