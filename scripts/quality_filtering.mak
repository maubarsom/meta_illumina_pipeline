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

ifndef step
$(warning Variable 'step' has been defined as 'qf')
step:=qf
endif

ifndef STRATEGY
$(info Default quality filtering strategy is cutadapt+nesoni+prinseq-lite)
STRATEGY=3_prinseq
endif

#Outfile
OUT_PREFIX := $(sample_name)_$(step)

#Reads
R1 := $(wildcard $(read_folder)/*R1*.f*q.gz  $(read_folder)/*_1.f*q.gz)
R2 := $(wildcard $(read_folder)/*R2*.f*q.gz  $(read_folder)/*_2.f*q.gz)

ifneq ($(words $(R1) $(R2)),2)
$(error More than one R1 or R2 $(words $(R1) $(R2)))
endif

#Logging
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

#Prinseq params
#out_format: fasta 1, fastq 3
#Remove duplicates (derep) : 1=exact dup, 2= 5' dup 3= 3' dup 4= rev comp exact dup
#Low complexity filters(lc_method): dust, entropy
prinseq_params:= -verbose -out_format 3 -log prinseq.log -min_len 75 -derep 1 -lc_method dust -lc_threshold 39

#SGA parameters
sga_ec_kmer := 41
sga_cov_filter := 2

#Output name generators (notice the = instead of := to set the appropriate directory)
nesoni_out_prefix = $(dir $@)$*
prinseq_out_prefix = $(dir $@)$*

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: $(OUT_PREFIX)_R1.fq.gz $(OUT_PREFIX)_R2.fq.gz $(OUT_PREFIX)_single.fq.gz

$(OUT_PREFIX)_%.fq.gz: $(STRATEGY)/$(sample_name)_%.fq.gz
	ln -fs $^ $@

#*************************************************************************
#Calls to trimmers
#*************************************************************************
1_cutadapt/$(sample_name)_R1.fq.gz: $(R1)
	mkdir -p $(dir $@)
	$(CUTADAPT_BIN) -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 -o $@ $^ >> $(log_file)

1_cutadapt/$(sample_name)_R2.fq.gz: $(R2)
	mkdir -p $(dir $@)
	$(CUTADAPT_BIN) -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 -o $@ $^ >>$(log_file)

#You have to specify quality is phred33 because with cutadapt clipped fragments nesoni fails to detect encoding
2_nesoni/%_R1.fq.gz 2_nesoni/%_R2.fq.gz 2_nesoni/%_single.fq.gz: 1_cutadapt/%_R1.fq.gz 1_cutadapt/%_R2.fq.gz
	mkdir -p $(dir $@)
	$(NESONI_BIN) clip --adaptor-clip no --homopolymers yes --qoffset 33 --quality 20 --length 75 \
		--out-separate yes $(nesoni_out_prefix) pairs: $^ 2>> $(log_file)

#Rule to plug to prinseq as it does not accept .gz files
2_nesoni/%_R1.fq 2_nesoni/%_R2.fq 2_nesoni/%_single.fq: 1_cutadapt/%_R1.fq.gz 1_cutadapt/%_R2.fq.gz
	mkdir -p $(dir $@)
	$(NESONI_BIN) clip --adaptor-clip no --homopolymers yes --qoffset 33 --quality 20 --length 75 \
		--out-separate yes --gzip no $(nesoni_out_prefix) pairs: $^ 2>> $(log_file)

3_prinseq/%_R1.fq.gz 3_prinseq/%_R2.fq.gz 3_prinseq/%_1_singletons.fastq 3_prinseq/%_2_singletons.fastq: 2_nesoni/%_R1.fq 2_nesoni/%_R2.fq
	mkdir -p $(dir $@)
	$(PRINSEQ_BIN) -fastq $< -fastq2 $(word 2,$^) $(prinseq_params) -out_good $(prinseq_out_prefix) -out_bad $(prinseq_out_prefix)_BAD 2>> $(log_file)
	mv $(prinseq_out_prefix)_1.fastq $(prinseq_out_prefix)_R1.fq && gzip $(prinseq_out_prefix)_R1.fq
	mv $(prinseq_out_prefix)_2.fastq $(prinseq_out_prefix)_R2.fq && gzip $(prinseq_out_prefix)_R2.fq
	if [ ! -e $(prinseq_out_prefix)_1_singletons.fastq ]; then touch $(prinseq_out_prefix)_1_singletons.fastq; fi
	if [ ! -e $(prinseq_out_prefix)_2_singletons.fastq ]; then touch $(prinseq_out_prefix)_2_singletons.fastq; fi

3_prinseq/%_single.fastq: 2_nesoni/%_single.fq
	mkdir -p $(dir $@)
	$(PRINSEQ_BIN) -fastq $^ $(prinseq_params) -out_good $(prinseq_out_prefix)_single -out_bad $(prinseq_out_prefix)_single_BAD 2>> $(log_file)

3_prinseq/%_single.fq.gz: 3_prinseq/%_single.fastq 3_prinseq/%_1_singletons.fastq 3_prinseq/%_2_singletons.fastq
	cat $^ | gzip > $@

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
	-rm *.log #Makefile log
	-rm *.log.txt #Nesoni log
