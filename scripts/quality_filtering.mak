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
$(warning Variable 'read_folder' will be assumed to be "./")
read_folder := ./
endif

ifndef STRATEGY
$(info Default quality filtering strategy is nesoni)
STRATEGY=nesoni_qf
endif

#Outfile
OUT_PREFIX := $(sample_name)_$(step)

#Reads
R1 := $(wildcard $(read_folder)/*R1*.f*q.gz  $(read_folder)/*_1.f*q.gz)
R2 := $(wildcard $(read_folder)/*R2*.f*q.gz  $(read_folder)/*_2.f*q.gz)

ifneq ($(words $(R1) $(R2)),2)
$(error More than one R1 or R2 $(words $(R1) $(R2)))
endif

OUT_PREFIX := $(sample_name)_qf

#Logging
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Run params
threads:=16

#SGA parameters
sga_ec_kmer := 41
sga_cov_filter := 2

#Binary paths
TRIMMOMATIC_BIN := java -jar /labcommon/tools/Trimmomatic-0.32/trimmomatic-0.32.jar
NESONI_BIN := nesoni
SGA_BIN := sga

#Output basenames
nesoni_1 := nesoni_q20hL75
nesoni_2 := nesoni_alt
trimmy_1 := trimmy
sga_1:= sga

nesoni_trimmy := nesoni_def_trimmy
nesoni_prinseq := nesoni_def_prinseq

#Output name generators (notice the = instead of := to set the appropriate directory)
nesoni_out_prefix = $(dir $@)$(sample_name)
trimmomatic_out = $(addprefix, $(dir $@)$(sample_name)_ , $(addsuffix .fq.gz R1 R1_single R2 R2_single ))

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: $(OUT_PREFIX)_R1.fq.gz $(OUT_PREFIX)_R2.fq.gz $(OUT_PREFIX)_single.fq.gz

$(OUT_PREFIX)_%.fq.gz: $(STRATEGY)/$(sample_name)_%.fq.gz
	ln -s $^ $@

#*************************************************************************
#Calls to trimmers
#*************************************************************************
nesoni_qf/%_R1.fq.gz nesoni_qf/%_R2.fq.gz nesoni_qf/%_single.fq.gz: $(R1) $(R2)
	mkdir -p $(dir $@)
	$(NESONI_BIN) clip --homopolymers yes --quality 20 --length 75 --out-separate yes \
		$(nesoni_out_prefix) pairs: $^ 2>> $(log_file)

nesoni_alt/%_R1.fq.gz nesoni_alt/%_R2.fq.gz nesoni_alt/%_single.fq.gz: $(R1) $(R2)
	mkdir -p $(dir $@)
	$(NESONI_BIN) clip --match 15 --max-errors 2 --homopolymers yes --quality 20 --length 75 \
		--out-separate yes $(nesoni_out_prefix) pairs: $^ 2>> $(log_file)

trimmy/$(sample_name)_R%.fq.gz: $(R1) $(R2)
	mkdir -p $(dir $@)
	$(TRIMMOMATIC_BIN) PE -threads 16 -phred33 -trimlog trimmomatic.log $^ $(trimmomatic_out) \
		ILLUMINACLIP:$(adaptor_file):2:30:10 MINLEN:75 2>> $(log_file)

#*************************************************************************
#Mixed quality filtering
#*************************************************************************
$(nesoni_trimmy)/R%.fq.gz: $(nesoni_1)/R1.fq.gz $(nesoni_1)/R2.fq.gz | $(nesoni_trimmy)
	$(TRIMMOMATIC_BIN) PE -threads 16 -phred33 -trimlog trimmomatic.log $^ $(trimmomatic_out) \
		ILLUMINACLIP:$(adaptor_file):2:30:10 MINLEN:75 2>> $(log_file)

$(nesoni_prinseq)/cosillo.fq: $(R1) $(R2) | $(nesoni_prinseq)
	prinseq-lite.pl -fastq $(R1) -fastq2 $(R2) -params prinseq_params.txt 2>> $(log_file)

#*************************************************************************
#SGA quality filtering steps
#*************************************************************************
#1) Preprocess the data to remove ambiguous basecalls
sga_qc/$(sample_name)_sga.fq: $(R1) $(R2)
	$(SGA_BIN) preprocess --pe-mode 1 -o $@ $^ 2>> $(log_file)

#2) Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) --no-reverse $(notdir $^) 2>> $(log_file)

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.fq: %.k$(sga_ec_kmer).ec.fq : %.fq %.sai
	$(SGA_BIN) correct -k $(sga_ec_kmer) --discard --learn -t $(threads) -o $@ $< 2>> $(log_file)

# Index the corrected data
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.sai: %.k$(sga_ec_kmer).ec.sai : %.k$(sga_ec_kmer).ec.fq
	$(SGA_BIN) index -a ropebwt -t $(threads) $^ 2>> $(log_file)

# Remove exact-match duplicates and reads with low-frequency k-mers
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.filter.pass.fa: $(sga_1).k$(sga_ec_kmer).ec.fq $(sga_1).k$(sga_ec_kmer).ec.sai
	$(SGA_BIN) filter -x $(sga_cov_filter) -t $(threads) --homopolymer-check \
		--low-complexity-check $< 2>> $(log_file)

.PHONY: clean
clean:
	-rm *.fq.gz
	-rm *.fq_fastqc.zip #Fastqc
	-rm *.log #Makefile log
	-rm *.log.txt #Nesoni log
