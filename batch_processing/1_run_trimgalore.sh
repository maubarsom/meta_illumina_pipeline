#!/bin/bash
set -euo verbose -o pipefail

READS_DIR=$1
OUT_DIR=$2

module load trimgalore

mkdir -p 1_trimgalore

trim_galore -q 20 --fastqc --fastqc_args "-k 10 -t 16" --illumina --paired --gzip \
	--stringency 5 --length 60 --output_dir ${OUT_DIR} --trim1 \
	--retain_unpaired -r1 85 -r2 85 $(ls ${READS_DIR}/*.fastq.gz | sort)

cp ${OUT_DIR}/*.txt ${OUT_DIR}/*fastqc.zip ${OUT_DIR}/*.html 1_trimgalore/
