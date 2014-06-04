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
$(error Variable 'basename' is not defined)
endif

ifndef step
$(error Variable 'step' is not defined)
endif

read_folder := ./

ifeq $(read_folder) "./"
	$(warning Variable 'read_folder' has been assumed to be "./")
endif

#Outfile
OUT_PREFIX:=$(sample_name)_$(step)

fastq_suffix:= .fastq .fastq.gz .fq .fq.gz
fastqc_suffix = .fastqc.zip fastqc.zip fq_fastqc.zip fq_fastqc.zip

#Reads
R1:= $(wildcard $(read_folder)/*R1*.  $(read_folder)/*_1.)

R2:= $(wildcard $(read_folder)/*R2*.fastq.gz $(read_folder)/*R2*.fq.gz )
R2+= $(wildcard $(read_folder)/*_2.fastq.gz $(read_folder)/*_2.fq.gz)

ifneq ($(words $(R1)) $(words $(R2)),1 1)
	$(error More than one R1 or R2)
endif

#Logging
log_name := $(CURDIR)/qc_$(step)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Run params
threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: $(OUT_PREFIX)_sga_preqc.pdf $(OUT_PREFIX)_k17.hist.pdf $(patsubst %.fastq.gz,%_fastqc.zip,$(R1) $(R2))

#*************************************************************************
#Import helper scripts - Check paths are ok
#*************************************************************************
plot_kmer_histogram.R:
	ln -s ../scripts/plot_kmer_histogram.R

#*************************************************************************
#FASTQC
#*************************************************************************
#@TODO Check if this thing actually works
%_fastqc.zip: %.fastq.gz
%_fastqc.zip: %.fastq
%.fq_fastqc.zip: %.fq.gz
%.fq_fastqc.zip: %.fq

%_fastqc.zip %.fq_fastqc.zip:
	fastqc --noextract -k 10 $^ 2>> $(log_file)

#*************************************************************************
#SGA PREQC
#*************************************************************************
# First, preprocess the data to remove ambiguous basecalls
$(OUT_PREFIX)_sga.fq: $(R1) $(R2)
	sga preprocess --pe-mode 1 -o $@ $^ >&2 2>> $(log_file)

# Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	sga index -a ropebwt -t $(threads) --no-reverse $(notdir $^) >&2 2>> $(log_file)

#Run SGA preqc
%_sga.preqc: %_sga.fq %_sga.sai
	sga preqc -t $(threads) --force-EM $< > $@ 2>> $(log_file)

#Create SGA preqc report
%_sga_preqc.pdf: %_sga.preqc
	sga-preqc-report.py -o $(basename $@) $^

#*************************************************************************
#JELLYFISH 2
#*************************************************************************
%_k17.jf: %_sga.fq
	jellyfish2 count -s 8G -C -m 17 -t $(threads) -o $@ $^ 2>> $(log_file)

%_k17.hist: %_k17.jf
	jellyfish2 histo -t $(threads) $^ -o $@ 2>> $(log_file)

%_k17.hist.pdf: %_k17.hist | plot_kmer_histogram.R
	Rscript plot_kmer_histogram.R $^ $@

#*************************************************************************
#CLEANING RULES
#*************************************************************************
.PHONY: clean-tmp clean-out

clean-tmp:
	-rm $(OUT_PREFIX)_sga.{fq,sai,bwt}
	-rm *.log #Makefile log
	-rm *.jf
	-rm plot_kmer_histogram.R

clean-out:
	-rm *_fastqc.zip #Fastqc
	-rm *.hist *.hist.pdf #Jellyfish
	-rm *.preqc *.pdf #SGA preqc
