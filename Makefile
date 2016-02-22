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

ifndef read_folder
read_folder := reads/
$(warning 'Read folder is assumed to be $(read_folder)')
endif

#Raw read detection parameters
ifndef raw_fq_ext
raw_fq_ext=fastq.gz
endif

ifndef raw_R1_filter
#raw_R1_filter=_R1_
raw_R1_filter=_1.
endif

ifndef raw_R2_filter
#raw_R2_filter=_R2_
raw_R2_filter=_2.
endif

#Run params from
ifndef cfg_file
$(error Config file variable 'cfg_file' not set)
endif
include $(cfg_file)

#Logging
timestamp := $(shell date +%s)
log_name = $(sample_name)_$@_$(timestamp).log
log_file = >(tee -a log/$(log_name) >&2)

.PHONY: all 1_raw_qc 2_qf 2_qf_qc 3_hostfiltering 4_assembly 5_tax_diamond 5_metaphlan

all: 1_raw_qc 2_qf 2_qf_qc 3_hostfiltering 4_assembly 5_tax_diamond 5_metaphlan

#QC raw reads
1_raw_qc: $(read_folder)
	mkdir -p $@/log
	cd $@ && $(MAKE) -rf ../steps/qc.mak sample_name=$(sample_name) read_folder=../reads/ step=raw fq_ext=$(raw_fq_ext) \
		R1_filter=$(raw_R1_filter) R2_filter=$(raw_R2_filter) RAW=1 fastqc &>> $(log_file)
	-cd $@ && $(MAKE) -rf ../steps/qc.mak read_folder=../reads/ step=raw clean-tmp

#Quality filtering
2_qf: $(read_folder)
	mkdir -p $@/log
	cd $@ && $(MAKE) -rf ../steps/quality_filtering.mak sample_name=$(sample_name) read_folder=../reads/ \
		fq_ext=$(raw_fq_ext) R1_filter=$(raw_R1_filter) R2_filter=$(raw_R2_filter) &>> $(log_file)

#QC Quality filtering
2_qf_qc: 2_qf
	mkdir -p $@/log
	cd $@ &&  $(MAKE) -rf ../steps/qc.mak sample_name=$(sample_name) read_folder=../$^/ step=qc fq_ext=fq.gz fastqc &>> $(log_file)
	-cd $@ && $(MAKE) -rf ../steps/qc.mak sample_name=$(sample_name) read_folder=../$^/ step=qc clean-tmp

#Contamination removal (human and phiX174)
3_hostfiltering: 2_qf
	mkdir -p $@/log
	cd $@ && $(MAKE) -rf ../steps/contamination_rm.mak sample_name=$(sample_name)_qf read_folder=../$^/ step=hf &>> $(log_file)

#Assembly step
4_assembly: 3_hostfiltering
	mkdir -p $@/log
	cd $@ && $(MAKE) -rf ../steps/assembly.mak sample_name=$(sample_name)_qf_hf read_folder=../$^/ step=asm &>> $(log_file)

#Taxonomic / Functional Annotation
5_tax_diamond: 4_assembly
	mkdir -p $@/log
	#Diamond blastx against NR
	cd $@ && $(MAKE) -rf ../steps/tax_assign.mak sample_name=$(sample_name)_qf_hf_asm in_folder=../$^/ step=tax megan &>> $(log_file)

5_metaphlan:
	mkdir -p $@/log
	cd $@ && $(MAKE) -rf ../steps/metaphlan.mak read_folder=../reads/ raw &>> $(log_file)

# 5_tax_blast: 4_assembly
# 	echo "This rule has not been implemented"
#
# 5_tax_kraken: 4_assembly
# 	echo "This rule has not been implemented"
