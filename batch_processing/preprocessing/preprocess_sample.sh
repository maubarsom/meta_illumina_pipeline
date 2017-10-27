#!/bin/bash
set -euo verbose -o pipefail
IFS=$'\n\t'

module load trimgalore/0.4.1
module load bowtie2/2.3.0
module load samtools/1.3

SAMPLE_ID=$1
READS_DIR=$2
SCRIPTS_DIR=/labcommon/viral_discovery/meta_illumina_pipeline/batch_processing/preprocessing

OLD_PATH=${PATH}
export PATH=${SCRIPTS_DIR}:${PATH}

make -j 3 -f ${SCRIPTS_DIR}/Makefile \
			sample_id=${SAMPLE_ID} \
			read_folder=${READS_DIR}

#Calculate checksum of output files
for x in $(ls *.fq.gz);
do
	md5sum ${x} > ${x}.md5
done

#Clean tmp folder if everything worked alright
make -f ${SCRIPTS_DIR}/Makefile sample_id=${SAMPLE_ID} clean

#Restore old path
export PATH=${OLD_PATH}
