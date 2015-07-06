# Pipeline for read QC inspection

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

#External parameters
# 1) basename
# 2) step
# 3) read_folder

ifndef sample_name
$(error Variable 'sample_name' is not defined)
endif

ifndef step
$(error Variable 'step' is not defined)
endif

ifndef read_folder
$(warning Variable 'read_folder' will be assumed to be "./")
read_folder := ./
endif

ifndef TMP_DIR
TMP_DIR := /tmp
endif

#Outfile
OUT_PREFIX:=$(sample_name)_$(step)

#Input read autodetection settings
fq_ext:=fastq.gz
#Only used for raw mode
R1_filter:=_1.
R2_filter:=_2.

#For raw mode
# filter_fx( subst, list) : Returns items in list that contain subst
filter_fx = $(foreach file,$(2),$(if $(findstring $(1),$(file)),$(file)))

R1_files := $(call filter_fx,$(R1_filter),$(wildcard $(read_folder)/*.$(fq_ext)))
R2_files := $(call filter_fx,$(R2_filter),$(wildcard $(read_folder)/*.$(fq_ext)))

#Definitions to create fastqc file targets depending on the fastq extension
fastqc_ext := $(if $(findstring fastq,$(fq_ext)),_fastqc.zip,.fq_fastqc.zip)
fastqc_txt_ext := $(if $(findstring fastq,$(fq_ext)),_fastqc.txt,.fq_fastqc.txt)
generate_fastqc_targets = $(patsubst %.$(fq_ext),%$(fastqc_txt_ext),$(notdir $(1)))
smart_cat := $(if $(findstring gz,$(fq_ext)),zcat,$(if $(findstring bz2,$(fq_ext)),bzcat,cat))

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

#Targets!

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all fastqc jellyfish2 sga

all: fastqc jellyfish2

ifdef RAW
fastqc: $(sample_name)_raw_1$(fastqc_txt_ext) $(sample_name)_raw_2$(fastqc_txt_ext)
else
fastqc: $(call generate_fastqc_targets, $(wildcard $(read_folder)/*.$(fq_ext)) )
endif

#Computationally intensive
jellyfish2: $(OUT_PREFIX)_k17.hist.pdf
sga: $(OUT_PREFIX)_sga_preqc.pdf

#********************************************************
#*********************RAW MODE **************************
#********************************************************
#Deals automatically with multiple R1 and R2 files ,
# merging them together into a single R1 and R2 file

ifdef RAW
ifneq "$(words $(R1_files))" "$(words $(R2_files))"
$(error Different number of R1 ($(words $(R1_files))) and R2 ($(words $(R2_files))) files)
endif

ifeq "0" "$(words $(R1_files))"
$(error "No R1 or R2 files")
endif

ifeq "1" "$(words $(R1_files))"
# If only one R1 and one R2 file, create links to TMP_DIR
$(TMP_DIR)/$(sample_name)_raw_1.$(fq_ext) : $(R1_files)
	ln -s $(shell pwd)/$^ $@

$(TMP_DIR)/$(sample_name)_raw_2.$(fq_ext) : $(R2_files)
	ln -s $(shell pwd)/$^ $@
else
# If more than one R1 and R2
# Merge all Forward and Reverse reads into a single file
$(TMP_DIR)/$(sample_name)_raw_1.$(fq_ext) : $(R1_files)
	cat $^ > $@

$(TMP_DIR)/$(sample_name)_raw_2.$(fq_ext) : $(R2_files)
	cat $^ > $@
#*********
endif

endif

#*************************************************************************
# FASTQC
#*************************************************************************
FASTQC_RECIPE = $(FASTQC_BIN) --noextract -k 10 -t $(threads) -o ./ $^

#QF mode
%$(fastqc_ext): $(read_folder)/%.$(fq_ext)
	$(FASTQC_RECIPE)

#Raw mode
$(sample_name)_raw_%$(fastqc_ext): $(TMP_DIR)/$(sample_name)_raw_%.$(fq_ext)
	$(FASTQC_RECIPE)

%$(fastqc_txt_ext): %$(fastqc_ext)
	unzip -p $<  $(basename $<)/fastqc_data.txt > $@

#*************************************************************************
#SGA PREQC - only in RAW mode
#*************************************************************************
# First, preprocess the data to remove ambiguous basecalls
$(OUT_PREFIX)_sga.fq: $(TMP_DIR)/$(sample_name)_raw_1.$(fq_ext) $(TMP_DIR)/$(sample_name)_raw_2.$(fq_ext)
	$(SGA_BIN) preprocess --pe-mode 1 -o $@ $^ >&2

# Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	$(SGA_BIN) index -a ropebwt -t $(threads) --no-reverse $(notdir $^) >&2

#Run SGA preqc
%_sga.preqc: %_sga.fq %_sga.sai
	$(SGA_BIN) preqc -t $(threads) --force-EM $< > $@

#Create SGA preqc report
%_sga_preqc.pdf: %_sga.preqc
	$(SGA_PREQC_REPORT_BIN) -o $(basename $@) $^

#*************************************************************************
#JELLYFISH 2 -
#*************************************************************************
plot_kmer_histogram.R:
	ln -s ../scripts/plot_kmer_histogram.R

%_k17.jf: $(wildcard $(read_folder)/*.$(fq_ext))
	$(smart_cat) $^ | $(JELLYFISH2_BIN) count -s 8G -C -m 17 -t $(threads) -o $@

%_k17.hist: %_k17.jf
	$(JELLYFISH2_BIN) histo -t $(threads) $^ -o $@

%_k17.hist.pdf: %_k17.hist | plot_kmer_histogram.R
	Rscript plot_kmer_histogram.R $^ $@

#*************************************************************************
#PRINSEQ
#*************************************************************************
%_pe_stats.txt: $(R1) $(R2)
	gunzip -c $< > $(TMP_DIR)/tmp_R1.fq
	gunzip -c $(word 2,$^) > $(TMP_DIR)/tmp_R2.fq
	$(PRINSEQ_BIN) -fastq $(TMP_DIR)/tmp_R1.fq -fastq2 $(TMP_DIR)/tmp_R2.fq -stats_all > $@
	rm $(TMP_DIR)/tmp_R1.fq $(TMP_DIR)/tmp_R2.fq

%_se_stats.txt: $(SINGLE)
	gunzip -c $< > $(TMP_DIR)/tmp_single.fq
	$(PRINSEQ_BIN) -fastq $(TMP_DIR)/tmp_single.fq -stats_all > $@
	rm $(TMP_DIR)/tmp_single.fq

#*************************************************************************
#CLEANING RULES
#*************************************************************************
.PHONY: clean-tmp clean-out

clean-tmp:
	-rm $(OUT_PREFIX)_sga.{fq,sai,bwt}
	-rm *.jf
	-rm plot_kmer_histogram.R

clean-out:
	-rm *_fastqc.zip #Fastqc
	-rm *.hist *.hist.pdf #Jellyfish
	-rm *.preqc *.pdf #SGA preqc
