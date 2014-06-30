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

# ifndef read_folder
# $(warning Variable 'read_folder' will be assumed to be "./")
# read_folder := ./
# endif

# ifndef STRATEGY
# $(info Default quality filtering strategy is cutadapt+nesoni+prinseq-lite)
# STRATEGY=3_prinseq
# endif

#Outfile
OUT_PREFIX := $(sample_name)_$(step)

#Reads
INPUT_PAIRED_END := $(read_folder)/$(sample_name)_qf_rmcont_pe.fq
INPUT_SINGLE_END := $(read_folder)/$(sample_name)_qf_rmcont_se.fq

OUT_PREFIX := $(sample_name)_qf

#Logging
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Run params
threads:=16

#Prinseq params
#out_format: fasta 1, fastq 3
#Remove duplicates (derep) : 1=exact dup, 2= 5' dup 3= 3' dup 4= rev comp exact dup
#Low complexity filters(lc_method): dust, entropy
prinseq_params:= -verbose -out_format 3 -log prinseq.log -min_len 75 -derep 1 -lc_method dust -lc_threshold 35

#SGA parameters
sga_ec_kmer := 41
sga_cov_filter := 2
sga_min_overlap := 200
sga_assemble_overlap := 111
TRIM_LENGTH := 400

#Binary paths
SGA_BIN := sga

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: raymeta_ctgs_filt.fa fermi_ctgs_filt.fa abyss_ctgs_filt.fa sga_ctgs_filt.fa

$(OUT_PREFIX)_%.fq.gz: $(STRATEGY)/$(sample_name)_%.fq.gz
	ln -s $^ $@

#*************************************************************************
#MetaRay
#*************************************************************************
raymeta/Contigs.fasta: $(INPUT_PAIRED_END)
	@echo -e "\nAssembling reads with Ray Meta\n\n" > $(log_file)
	mpiexec -n 16 Ray Meta -i $^ -o $(dir $@) 2> $(log_file)

raymeta_contigs.fa: raymeta/Contigs.fasta
	ln -s $^ $@

#*************************************************************************
#Fermi
#*************************************************************************
#Runs Fermi assembler until the 4th step.
#Qualities are not neccessary for Kraken/Blast classification
fermi/fmdef.p4.fa.gz: $(INPUT_PAIRED_END)
	mkdir -p fermi
	cd fermi && run-fermi.pl -t $(threads) -c ../$^ > assembly.mak 2> $(log_file)
	cd fermi && $(MAKE) -f assembly.mak -j $(threads) $(notdir $@) 2> $(log_file)

fermi_contigs.fa: fermi/fmdef.p4.fa.gz
	gunzip -c $^ > $@

#*************************************************************************
#Abyss
#*************************************************************************
abyss/abyssk47-contigs.fa: $(INPUT_PAIRED_END)
	mkdir -p abyss
	cd abyss && abyss-pe k=47 name='abyssk47' in='../$^' np=16 2> $(log_file)

abyss_contigs.fa: abyss/abyssk47-contigs.fa
	ln -s $^ $@

#*************************************************************************
#SGA quality filtering steps
#*************************************************************************
#1) Preprocess the data to remove ambiguous basecalls
#pe-mode=2 : interleaved
sga/$(sample_name)_sga.fq: $(INPUT_PAIRED_END)
	$(SGA_BIN) preprocess --pe-mode 2 -o $@ $^ 2>> $(log_file)

#2) Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) --no-reverse $(notdir $^) 2>> $(log_file)

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
%.k$(sga_ec_kmer).ec.fq : %_sga.fq %_sga.sai
	$(SGA_BIN) correct -k $(sga_ec_kmer) --discard --learn -t $(threads) -o $@ $< 2>> $(log_file)

# Index the corrected data
%.ec.sai: %.ec.fq
	$(SGA_BIN) index -a ropebwt -t $(threads) $^ 2>> $(log_file)

# Remove exact-match duplicates and reads with low-frequency k-mers
%.ec.filter.pass.fa: %.ec.fq %.ec.sai
	$(SGA_BIN) filter -x $(sga_cov_filter) -t $(threads) --homopolymer-check \
		--low-complexity-check $< 2>> $(log_file)

# Compute the structure of the string graph
%.ec.filter.pass.asqg.gz :%.ec.filter.pass.fa
	$(SGA_BIN) overlap -m $(sga_min_overlap) -t $(threads) $^ 2>> $(log_file)

# Perform the contig assembly
%-contigs.fa: %.ec.filter.pass.asqg.gz
	$(SGA_BIN) assemble -m $(sga_assemble_overlap) --min-branch-length $(TRIM_LENGTH) -o $* $^ 2>> $(log_file)

sga_contigs.fa: sga/$(sample_name)-contigs.fa
	ln $^ $@

#*************************************************************************
#Extract contigs > 550 bp
#*************************************************************************
#Keep only contigs greater than 750bp
%_ctgs_filt.fa : %_contigs.fa
	seqtk seq -L 500 $^ > $@ 2> $(log_file)

.PHONY: clean
clean:
	-rm *.fq.gz
	-rm *.fq_fastqc.zip #Fastqc
	-rm *.log #Makefile log
	-rm *.log.txt #Nesoni log
