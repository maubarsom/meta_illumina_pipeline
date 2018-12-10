#!/bin/bash
set -euo verbose -o pipefail
IFS=$'\n\t'

BOWTIE_DB=/labcommon/db/bowtie2/grch38

SAMPLE_ID=$1
IN_DIR=$2

OUT_FILE_R1=${SAMPLE_ID}_1.fq.gz
OUT_FILE_R2=${SAMPLE_ID}_2.fq.gz

if [ -s ${OUT_FILE_R1} ] || [ -s ${OUT_FILE_R2} ]; then
	echo "File ${OUT_FILE_R1} and/or ${OUT_FILE_R2} already exist. Delete file(s) to regenerate."
	exit 0;
fi

#Create named PIPEs
mkfifo R1_PIPE R2_PIPE

mkdir -p log

#Requires cutadapt
cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -f 'fastq' \
	-o R1_PIPE ${IN_DIR}/*_val_1.fq.gz 2>&1 > log/${SAMPLE_ID}_overhangs_R1.log &

cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -f 'fastq' \
	-o R2_PIPE ${IN_DIR}/*_val_2.fq.gz 2>&1 > log/${SAMPLE_ID}_overhangs_R2.log &

bowtie2 --local --very-sensitive-local -t -p 12 -x ${BOWTIE_DB} -1 R1_PIPE -2 R2_PIPE | \
	tee >(samtools view -hSb - | samtools flagstat - > log/${SAMPLE_ID}_bowtie_pe.flagstat ) | \
	samtools view -hSb -f12 -F256 - | samtools fastq -1 ${OUT_FILE_R1} -2 ${OUT_FILE_R2} -

#Remove named PIPEs
rm R1_PIPE R2_PIPE
