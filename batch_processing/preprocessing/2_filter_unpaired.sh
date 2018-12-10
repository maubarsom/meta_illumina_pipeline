#!/bin/bash
set -euo verbose -o pipefail
IFS='\n\t'

BOWTIE_DB=/labcommon/db/bowtie2/grch38

SAMPLE_ID=$1
IN_DIR=$2

OUT_FILE=${SAMPLE_ID}_unpaired.fq.gz

#If output file exists, do nothing
if [ -s "${OUT_FILE}" ]; then
	echo "File ${OUT_FILE} already exists. Delete file to regenerate"
	exit 0;
fi

#Create named PIPEs
mkfifo SE_PIPE

mkdir -p log

#Merge unpaired into a single file
cat ${IN_DIR}/*_unpaired_2.fq.gz >> ${IN_DIR}/*_unpaired_1.fq.gz 
rm ${IN_DIR}/*_unpaired_2.fq.gz

#Requires cutadapt
cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 \
	-o SE_PIPE ${IN_DIR}/*_unpaired_1.fq.gz 2>&1 > log/${SAMPLE_ID}_overhangs_unpaired.log &

bowtie2 --local --very-sensitive-local -t -p 8 -x ${BOWTIE_DB} -U SE_PIPE | \
	tee >(samtools view -hSb - | samtools flagstat - > log/${SAMPLE_ID}_bowtie_unpaired.flagstat ) | \
	samtools view -hSb -f4 -F256 - | samtools fastq - | gzip > ${OUT_FILE}

#Remove named PIPEs
rm SE_PIPE
