#!/bin/bash
set -euo verbose -o pipefail

SAMPLE_ID=$1
READS_DIR=$2
SCRIPTS_DIR=/labcommon/viral_discovery/meta_illumina_pipeline/batch_processing

OLD_PATH=${PATH}
export PATH=${SCRIPTS_DIR}:${PATH}

make -j 3 -f ${SCRIPTS_DIR}/Makefile \
			sample_id=${SAMPLE_ID} \
			read_folder=${READS_DIR}

export PATH=${OLD_PATH}
