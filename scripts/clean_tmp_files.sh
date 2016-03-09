#!/bin/bash
if [ -z pipeline_stats.ipynb ];
then
	echo "You should run the pipeline stats before cleaning the tmp files"
	exit 1
fi

#***********************
# Tmp folder
#***********************
rm -r tmp/

#***********************
# QUALITY FILTERING
#***********************
#Remove temporary qf files
find 2_qf -name "*.fq.gz" -delete

#***********************
# HOST FILTERING
#***********************
#Remove files mapping to human reference
rm -r 3_hostfiltering/bowtie2
#Compress effective reads
find 3_hostfiltering -name "*.fq" -exec gzip {} \;

#***********************
# ASSEMBLY
#***********************
#Remove assembler-specific folders
rm -r 4_assembly/fermi
rm -r 4_assembly/megahit
rm -rf 4_assembly/spades

#Remove tmp mapping files to determine non-incorporated reads
#All reads are kept from 3_hostfiltering step
rm -r 4_assembly/singletons
#Remove softlinks to mapped reads to assembly
unlink 4_assembly/*.fa

#***********************
# Tax assign
#***********************
#Remove diamond alignment files
rm 5_tax_diamond/diamond/*.daa
