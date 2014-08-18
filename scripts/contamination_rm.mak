# Human contamination removal pipeline

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

#Make parameters
SHELL := /bin/bash

ifndef sample_name
$(error Variable 'sample_name' is not defined)
endif

ifndef read_folder
$(error Variable 'read_folder' is not defined)
endif

ifndef prev_steps
prev_steps := qf
$(info 'prev_steps' is assumed to be $(prev_steps))
endif

ifndef step
step:=rmcont
$(warning Variable 'step' has been defined as '$(step)')
endif

ifndef STRATEGY
$(error Variable 'STRATEGY' is not defined. Select 'bwa' or 'bwastampy')
endif

#Input and Output file prefixes
IN_PREFIX:= $(sample_name)_$(prev_steps)
OUT_PREFIX:= $(IN_PREFIX)_$(step)

#Run parameters
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

#Databases
stampy_grch38 := /labcommon/db/stampy/GRCh38/grch38
stampy_grch37 := /labcommon/db/stampy/GRCh37.p13/grch37_p13
bwa_grch37 := /labcommon/db/bwa/GRCh37.p13/grch37_p13

#Input files
R1 := $(read_folder)/$(IN_PREFIX)_R1.fq.gz
R2 := $(read_folder)/$(IN_PREFIX)_R2.fq.gz
singles := $(read_folder)/$(IN_PREFIX)_single.fq.gz

#Log file
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Intermediate file names
bwa_pre:= $(IN_PREFIX)_bwa
bwastampy_pre := $(IN_PREFIX)_bwastampy

mapping_pre := $(IN_PREFIX)_$(STRATEGY)
sort_pre := $(mapping_pre)_sort
filter_pre := $(sort_pre)_filter

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: $(OUT_PREFIX)_pe.fq $(OUT_PREFIX)_se.fq

#*************************************************************************
#Create output files from the strategy
#*************************************************************************
$(OUT_PREFIX)_%.fq: $(filter_pre)_%.fq
	ln -s $^ $@

#*************************************************************************
#Map to GRCh37 with BWA MEM
#*************************************************************************
$(bwa_pre)_pe.sam: $(R1) $(R2)
$(bwa_pre)_se.sam: $(singles)

$(bwa_pre)_pe.sam $(bwa_pre)_se.sam:
	@echo -e "\nMapping $^ to GRCh37 with BWA MEM @"`date`"\n\n" >> $(log_file)
	bwa mem -t $(threads) -T 20 -M $(bwa_grch37) $^ > $@ 2> $(log_file)

#Prepare bwa output for Stampy
$(bwa_pre)_pe.bam $(bwa_pre)_se.bam: $(bwa_pre)_%.bam: $(bwa_pre)_%.sam
	samtools view -hSb -o $@ $^ 2> $(log_file)

#*************************************************************************
#Improve BWA mappings with Stampy
#*************************************************************************
$(bwastampy_pre)_%.sam: $(bwa_pre)_%.bam
	@echo -e "\nCorrecting $^ with Stampy @"`date`"\n\n" >> $(log_file)
	stampy.py -t $(threads) -g $(stampy_grch37) -h $(stampy_grch37) --bamkeepgoodreads -o $@ -M $^

#*************************************************************************
#Convert to bam removing secondary mappings
#*************************************************************************
$(bwastampy_pre)_pe.bam $(bwastampy_pre)_se.bam: %.bam:%.sam
	@echo -e "\nRemoving secondary mappings from $^ @ `date` \n\n" >> $(log_file)
	samtools view -F 256 -hSb -o $@ $^ 2> $(log_file)

#*************************************************************************
#Extract unmapped reads using Picard Tools
#*************************************************************************
$(sort_pre)_pe.bam $(sort_pre)_se.bam: $(sort_pre)_%.bam : $(mapping_pre)_%.bam
	@echo -e "\nSort bam file by queryname\n\n" >> $(log_file)
	run_picard SortSam.jar INPUT=$^ OUTPUT=$@ SORT_ORDER=queryname 2> $(log_file)

#Keep only reads that did not map confidently (with both pairs)
$(filter_pre)_pe.bam $(filter_pre)_se.bam: $(filter_pre)_%.bam : $(sort_pre)_%.bam
	@echo -e "\nExtract $* reads that did not map to GRCh37\n\n" > $(log_file)
	run_picard FilterSamReads.jar INPUT=$^ OUTPUT=$@ FILTER=excludeAligned SORT_ORDER=queryname WRITE_READS_FILES=False 2> $(log_file)

$(filter_pre)_%.fq : $(filter_pre)_%.bam
	@echo -e "\nConvert bam to interleaved fastq\n\n" > $(log_file)
	run_picard SamToFastq.jar INPUT=$^ FASTQ=$@ INTERLEAVE=True 2> $(log_file)
