#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

preprocess_script=/labcommon/viral_discovery/meta_illumina_pipeline/batch_processing/preprocessing/preprocess_sample.sh

for batch_id in $(cd /media/virushd/samples/illumina; ls -d B.Andersson_*);
do
for sample_id in $(cd /media/virushd/samples/illumina/${batch_id}; ls -d P*);
do
	mkdir -p ${sample_id}
	cd ${sample_id}
	bash ${preprocess_script} ${sample_id} /media/virushd/samples/illumina/${batch_id}/${sample_id}
	cd ../
done
