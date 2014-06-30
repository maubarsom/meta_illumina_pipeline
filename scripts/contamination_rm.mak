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

ifndef qf_prefix
qf_prefix := $(sample_name)_qf
$(info 'qf_prefix' is assumed to be $(qf_prefix))
endif

threads := 16
stampy_grch38 := /labcommon/db/stampy/GRCh38/grch38
stampy_grch37 := /labcommon/db/stampy/GRCh37.p13/grch37_p13
bwa_grch37 := /labcommon/db/bwa/GRCh37.p13/grch37_p13

#Input files
R1 := $(read_folder)/$(qf_prefix)_R1.fq.gz
R2 := $(read_folder)/$(qf_prefix)_R2.fq.gz
singles := $(read_folder)/$(qf_prefix)_single.fq.gz

#Output files
OUT_PREFIX:= $(sample_name)_qf_rmcont

#Log file
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Intermediate file names
bwa_pre := $(qf_prefix)_bwa
stampy_pre := $(qf_prefix)_stampy
bwastampy_pre := $(qf_prefix)_bwastampy

sortsam_pre := $(bwastampy_pre)_sort
filtersam_pre := $(sortsam_pre)_filter
sam2fq_pre := $(filtersam_pre)_sam2fq

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: $(filtersam_pre)_pe.fq
# all: $(OUT_PREFIX)_pe.fq

#*************************************************************************
#Map to GRCh37 with BWA MEM
#*************************************************************************
$(bwa_pre)_pe.sam: $(R1) $(R2)
$(bwa_pre)_se.sam: $(single_input)

$(bwa_pre)_pe.sam $(bwa_pre)_se.sam:
	@echo -e "\nMapping $^ to GRCh37 with BWA MEM\n\n" > $(log_file)
	bwa mem -t $(threads) -T 20 -M $(bwa_grch37) $^ > $@ 2> $(log_file)

#*************************************************************************
#Improve BWA mappings with Stampy
#*************************************************************************
$(bwastampy_pre)_%.sam: $(bwa_pre)_%.bam
	stampy.py -t $(threads) -g $(stampy_grch37) -h $(stampy_grch37) --bamkeepgoodreads -o $@ -M $^

#*************************************************************************
#Map to GRCh37 with Stampy
#************************************************************************
$(stampy_pre)_pe.sam: $(R1) $(R2)
$(stampy_pre)_se.sam: $(singles)

$(stampy_pre)_pe.sam $(stampy_pre)_se.sam:
	stampy.py -t $(threads) -g $(stampy_grch37) -h $(stampy_grch37) -o $@ -M $^ 2>> $(log_file)

#Convert from sam to bam removing secondary mappings
%.bam: %.sam
	@echo -e "Converting $^ to .bam\n\n" > $(log_file)
	samtools view -F 256 -hSb -o $@ $^ 2> $(log_file)

#*************************************************************************
#Extract unmapped reads using Picard Tools
#*************************************************************************
$(sortsam_pre)_pe.bam $(sortsam_pre)_se.bam: $(sortsam_pre)_%.bam : $(bwastampy_pre)_%.bam
	@echo -e "\nSort bam file by queryname\n\n" > $(log_file)
	run_picard SortSam.jar INPUT=$^ OUTPUT=$@ SORT_ORDER=queryname 2> $(log_file)

#Keep only reads that did not map confidently (with both pairs)
$(filtersam_pre)_pe.bam $(filtersam_pre)_se.bam: $(filtersam_pre)_%.bam : $(sortsam_pre)_%.bam
	@echo -e "\nExtract $* reads that did not map to GRCh37\n\n" > $(log_file)
	run_picard FilterSamReads.jar INPUT=$^ OUTPUT=$@ FILTER=excludeAligned SORT_ORDER=queryname WRITE_READS_FILES=False 2> $(log_file)

#Convert unmapped reads to Fastq for assembly
#%_R1.fq %_R2.fq: $(filtersam_pre).bam
#	run_picard SamToFastq.jar INPUT=$^ FASTQ=$(OUT_PREFIX)_R1.fq SECOND_END_FASTQ=$(OUT_PREFIX)_R2.fq 2>> $(log_file)

#Alternatively use interleaved paired-end fastq format
# $(OUT_PREFIX)_pe.fq $(OUT_PREFIX)_se.fq: $(OUT_PREFIX)_%.fq : $(filtersam_pre)_%.bam
$(filtersam_pre)_%.fq : $(filtersam_pre)_%.bam
	@echo -e "\nConvert bam to interleaved fastq\n\n" > $(log_file)
	run_picard SamToFastq.jar INPUT=$^ FASTQ=$@ INTERLEAVE=True 2> $(log_file)
