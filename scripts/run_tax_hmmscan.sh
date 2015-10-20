#!/bin/bash
set -euo pipefail

sample_name=$(basename assembly/*_allctgs.fa _allctgs.fa)
threads=16
#First param is number of procs
if [ -z $0 ]; then threads=$0; fi

echo "Processing ${sample_name} with ${threads} procs"

mkdir -p tax_assign/hmmscan/input

if [ ! -e tax_assign/hmmscan/input/${sample_name}_diamond_nohits.fa ];
then
cat <(python scripts/extract_diamond_nohits.py assembly/*_allctgs.fa tax_assign/diamond/*_allctgs_diamond_nr.sam*) \
	<(python scripts/extract_diamond_nohits.py assembly/*_se.fa tax_assign/diamond/*_se_diamond_nr.sam*)  | seqtk seq -L 400 > tax_assign/hmmscan/input/${sample_name}_diamond_nohits.fa
fi

cd tax_assign/hmmscan
make -f ../../steps/tax_hmmscan.mak SAMPLE=${sample_name} threads=${threads} in_fasta=input/${sample_name}_diamond_nohits.fa
