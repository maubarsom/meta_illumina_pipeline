#!/bin/bash

set -euo verbose -o pipefail

module load bowtie2
module load samtools

BOWTIE_DB=/labcommon/db/bowtie2/grch38

SAMPLE_ID=$1
IN_DIR=$2

#Create named PIPEs
mkfifo R1_PIPE R2_PIPE

#Requires cutadapt
cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 \
	-o R1_PIPE $(IN_DIR)/*_val_1.fq.gz 2>&1 > overhangs_R1.log &

cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 \
	-o R2_PIPE $(IN_DIR)/*_val_2.fq.gz 2>&1 > overhangs_R2.log &

bowtie2 --local --very-sensitive-local -t -p 8 -x ${BOWTIE_DB} -1 R1_PIPE -2 R2_PIPE | tee >(samtools view -hSb - | samtools flgstat - > bowtie_pe.flagstat ) | samtools view -hSb -f12 -F256 - | samtools fastq - | gzip > ${SAMPLE_ID}_pe.fq.gz

#Remove named PIPEs
rm R1_PIPE R2_PIPE
