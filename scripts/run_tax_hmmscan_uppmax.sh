#!/bin/bash
#SBATCH -A b2011088
#SBATCH -p core -n 4
#SBATCH -t 3:00:00
#SBATCH -J sample_hmmscan
#SBATCH --mail-user christian.pou@ki.se
#SBATCH --mail-type ALL
#SBATCH -e log/sample_hmmscan-%j.err
#SBATCH -o log/sample_hmmscan-%j.out

set -euo pipefail

module load bioinfo-tools
module load hmmer/3.1b2
module load emboss/6.5.7
module load seqtk

sample_name=$(basename 4_assembly/*_contigs.fa _contigs.fa)
threads=16
#First param is number of procs
if [ -z $0 ]; then threads=$0; fi

echo "Processing ${sample_name} with ${threads} procs"

mkdir -p 5_tax_diamond/hmmscan/input

if [ ! -e 5_tax_diamond/hmmscan/input/${sample_name}_diamond_nohits.fa ];
then
cat <(python scripts/extract_diamond_nohits.py 4_assembly/*_contigs.fa 5_tax_diamond/diamond/*_contigs_diamond_nr.sam.gz) \
	<(python scripts/extract_diamond_nohits.py 4_assembly/*_merged.fa 5_tax_diamond/diamond/*_merged_diamond_nr.sam.gz)  | seqtk seq -L 350 - > 5_tax_diamond/hmmscan/input/${sample_name}_diamond_nohits.fa
fi

if [[ -s 5_tax_diamond/hmmscan/input/${sample_name}_diamond_nohits.fa ]];
then
cd 5_tax_diamond/hmmscan
make -rf ../../steps/tax_hmmscan.mak SAMPLE=${sample_name} threads=${threads} UPPMAX=1 in_fasta=input/${sample_name}_diamond_nohits.fa
else
echo "File 5_tax_diamond/hmmscan/input/${sample_name}_diamond_nohits.fa is empty!" >&2
fi
