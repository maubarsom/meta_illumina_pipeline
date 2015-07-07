# Raw reads quality filtering pipeline

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

ifndef step
$(warning Variable 'step' has been defined as 'qf')
step:=qf
endif

ifndef STRATEGY
$(info Default quality filtering strategy is nesoni+cutadapt+flash)
STRATEGY=3_mergepairs
endif

#Outfile
OUT_PREFIX := $(sample_name)_$(step)

#Input read autodetection settings
R1_filter:=_1.
R2_filter:=_2.
fq_ext:=fastq.gz

#filter_fx( substring, list): Returns all items in list that contains substring
filter_fx = $(foreach file,$(2),$(if $(findstring $(1),$(file)),$(file)))

#Reads
R1 := $(call filter_fx,$(R1_filter),$(wildcard $(read_folder)/*.$(fq_ext)))
R2 := $(call filter_fx,$(R2_filter),$(wildcard $(read_folder)/*.$(fq_ext)))


ifneq ($(words $(R1) $(R2)),2)
$(error More than one R1 or R2 $(words $(R1) $(R2)))
endif

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

#Triming parameters
MIN_READ_LENGTH := 50

#Pair merging parameters
pairmerge_min_ovlp := 10

#Output name generators (notice the = instead of := to set the appropriate directory)
nesoni_out_prefix = $(dir $@)$*

#Delete produced files if step fails
.DELETE_ON_ERROR:
.SECONDARY:

.PHONY: all

all: $(OUT_PREFIX)_R1.fq.gz $(OUT_PREFIX)_R2.fq.gz $(OUT_PREFIX)_single.fq.gz

$(OUT_PREFIX)_%.fq.gz: $(STRATEGY)/$(sample_name)_%.fq.gz
	ln -fs $^ $@

#*************************************************************************
#Calls to trimmers
#*************************************************************************
#Light quality trimming, phred > 5
1_nesoni/%_R1.fq.gz 1_nesoni/%_R2.fq.gz 1_nesoni/%_single.fq.gz: $(R1) $(R2)
	mkdir -p $(dir $@)
	$(NESONI_BIN) clip --adaptor-clip no --homopolymers yes --qoffset 33 --quality 5 --length $(MIN_READ_LENGTH) \
		--out-separate yes $(nesoni_out_prefix) pairs: $^

2_cutadapt/%_R1.fq.gz 2_cutadapt/%_R2.fq.gz 2_cutadapt/%_single.fq.gz: 1_nesoni/%_R1.fq.gz 1_nesoni/%_R2.fq.gz 1_nesoni/%_single.fq.gz
	mkdir -p $(dir $@)
	#Remove Illumina double(or single index) adapters from fwd and rev pairs
	$(CUTADAPT_BIN) -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
					--overlap=5 --error-rate=0.1 --minimum-length $(MIN_READ_LENGTH) \
					-o $(TMP_DIR)/cutadapt_r1.fq.gz -p $(TMP_DIR)/cutadapt_r2.fq.gz $< $(word 2,$^)
	#Remove Illumina double(or single index) adapters from singletons
	$(CUTADAPT_BIN) -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
					--overlap=5 --error-rate=0.1 --minimum-length $(MIN_READ_LENGTH) \
					-o $(TMP_DIR)/cutadapt_single.fq.gz $(word 3,$^)
	#Remove PCR overhang adapter for all reads
	$(CUTADAPT_BIN) --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -o $(dir $@)/$*_R1.fq.gz $(TMP_DIR)/cutadapt_r1.fq.gz
	$(CUTADAPT_BIN) --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -o $(dir $@)/$*_R2.fq.gz $(TMP_DIR)/cutadapt_r2.fq.gz
	$(CUTADAPT_BIN) --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -o $(dir $@)/$*_single.fq.gz $(TMP_DIR)/cutadapt_single.fq.gz
	-rm $(TMP_DIR)/cutadapt_*.fq.gz

3_mergepairs/%_R1.fq.gz 3_mergepairs/%_R2.fq.gz 3_mergepairs/%_single.fq.gz: 2_cutadapt/%_R1.fq.gz 2_cutadapt/%_R2.fq.gz 2_cutadapt/%_single.fq.gz
	mkdir -p $(dir $@)
	$(FLASH_BIN) -m $(pairmerge_min_ovlp) -M $(READ_LEN)  $(word 1,$^) $(word 2,$^) -o out -d $(dir $@)
	cd $(dir $@) && mv out.notCombined_1.fastq $*_R1.fq && gzip $*_R1.fq
	cd $(dir $@) && mv out.notCombined_2.fastq $*_R2.fq && gzip $*_R2.fq
	cd $(dir $@) && mv out.extendedFrags.fastq $*_merged.fq
	cd $(dir $@) && gunzip -c ../$(word 3,$^) | cat $*_merged.fq | gzip > $*_single.fq.gz

.PHONY: clean
clean:
	-rm *.fq.gz
	-rm *.log #Makefile log
	-rm *.log.txt #Nesoni log
