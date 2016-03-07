#!/bin/bash
if [ -z pipeline_stats.ipynb ];
then
	echo "You should run the pipeline stats before cleaning the tmp files"
	exit 1
fi

#Remove tmp folder
rm -r tmp/

#Remove temporary qf files
find 2_qf -name "*.fq.gz" -delete

#Remove
rm 3_hostfiltering/*.fq

#Remove assembler-specific folders
rm -r 4_assembly/fermi
rm -r 4_assembly/megahit
rm -rf 4_assembly/spades

#Remove tmp mapping files to determine non-incorporated reads
rm 4_assembly/singletons/*.sam
