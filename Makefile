# Viral discovery pipeline for Illumina (MiSeq) data

# Author: Mauricio Barrientos-Somarribas
# Email:  mauricio.barrientos@ki.se

# Copyright 2014 Mauricio Barrientos-Somarribas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

SHELL := /bin/bash

#Required variables

ifndef sample_name
$(error Variable sample_name not set.)
endif

export sample_name

input_files := $(wildcard reads/*.fastq.gz) $(wildcard reads/*.fq.gz) $(wildcard reads/*.fq) $(wildcard reads/*.fastq)

ifneq "$(words $(input_files))" "2"
$(error Invalid number of paired-end read files in reads folder)
endif

#Run params
threads:=16

#Databases
db_path:=/labcommon/db
bwa_hg_ref:= $(db_path)/iGenomes/Homo_Sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
blastdb_folder:=$(db_path)/blastdb
kraken_db:=$(db_path)/krakendb/kraken140311/
swissprot_fasta:=$(db_path)/fasta/uniprot_sprot.fasta

#Paths to binaries
FragGeneScan_bin:=/labcommon/tools/FragGeneScan1.18/run_FragGeneScan.pl

#Tool parameters
blast_params:= -evalue 1 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

#Logging info
export log_name := $(CURDIR)$(sample_name)_$(shell date +%s).log
export log_file := >( tee -a $(log_name) >&2 )

#Prefixes for output files
qf_prefix := $(sample_name)_q20h
bwa_pre := $(nesoni_pre)_grch37
sortsam_pre := $(bwa_pre)_sort
filtersam_pre := $(sortsam_pre)_filter
sam2fq_pre := $(filtersam_pre)_sam2fq
diginorm_pre:=$(sam2fq_pre)_dgnrm

#Outputs
kraken_reports:= kraken_raymeta.report kraken_fermi.report kraken_abyss.report
phmmer_files := raymeta_fgs_phmmer.tbl fermi_fgs_phmmer.tbl abyss_fgs_phmmer.tbl
blastp_out := raymeta_fgs_blastp.xml fermi_fgs_blastp.xml abyss_fgs_blastp.xml reads_pe_fgs_blastp.xml reads_se_fgs_blastp.xml
blastx_out := raymeta_blastx.xml fermi_blastx.xml abyss_blastx.xml

#Avoid the parallel execution of rules in this makefile
.NOTPARALLEL:

.PHONY: all raw_qc qf_qc quality_filtering contamination_rm

all: raw_qc qf_qc contamination_rm

#QC raw reads
raw_qc: $(input_files)
	mkdir -p $@
	if [ ! -r $@/qc.mak ]; then cp scripts/qc.mak $@/; fi
	cd raw_qc && $(MAKE) -f qc.mak read_folder=../reads/ step=raw

#Quality filtering
quality_filtering: $(input_files)
	mkdir -p $@
	if [ ! -r quality_filtering/quality_filtering.mak ]; then cp scripts/quality_filtering.mak $@/; fi
	cd $@ && $(MAKE) -f quality_filtering.mak read_folder=../reads/

#QC Quality filtering
qf_qc: quality_filtering
	mkdir -p $@
	if [ ! -r $@/qc.mak ]; then cp scripts/qc.mak $@; fi
	cd $@ && $(MAKE) -f qc.mak read_folder=../quality_filtering/ step=qf

#Contamination removal (human)
contamination_rm: quality_filtering
	mkdir -p $@
	if [ ! -r $@/contamination_rm.mak ]; then cp scripts/contamination_rm.mak $@; fi
	cd $@ && $(MAKE) -f contamination_rm.mak read_folder=../quality_filtering/ step=rmcont

#Assembly step

#Post-Assembly - Separate singletons

#Taxonomic / Functional Annotation
