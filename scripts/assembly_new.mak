# Assembly strategies pipeline

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

ifndef prev_steps
prev_steps := qf_rmcont
$(warning 'prev_steps' is assumed to be $(prev_steps))
endif

ifndef step
step:= asm
$(warning Variable 'step' has been defined as '$(step)')
endif

#Input and Output file prefixes
IN_PREFIX := $(sample_name)_$(prev_steps)
OUT_PREFIX:= $(IN_PREFIX)_$(step)

#Reads
INPUT_PAIRED_END := $(read_folder)/$(IN_PREFIX)_pe.fq
INPUT_SINGLE_END := $(read_folder)/$(IN_PREFIX)_se.fq


#Logging
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Run params
threads:=16

#Abyss parameters
abyss_kmer:= 30

#SGA parameters
sga_ec_kmer := 41
sga_cov_filter := 2
sga_min_overlap := 200
sga_assemble_overlap := 111
TRIM_LENGTH := 400

#Masurca parameters
masurca_kmer := 21
masurca_pe_stats := 400 80
masurca_se_stats := 250 50
#Binary paths
SGA_BIN := sga
MASURCA_BIN:=/labcommon/tools/MaSuRCA-2.2.1/bin

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: $(OUT_PREFIX)_raymeta_ctgs_filt.fa $(OUT_PREFIX)_fermi_ctgs_filt.fa $(OUT_PREFIX)_sga_ctgs_filt.fa $(OUT_PREFIX)_masurca_ctgs_filt.fa

$(OUT_PREFIX)_%.fq.gz: $(STRATEGY)/$(sample_name)_%.fq.gz
	ln -s $^ $@

#*************************************************************************
#MetaRay
#*************************************************************************
raymeta/Contigs.fasta: $(INPUT_PAIRED_END)
	@echo -e "\nAssembling reads with Ray Meta\n\n" > $(log_file)
	mpiexec -n 16 Ray Meta -i $^ -o $(dir $@) 2> $(log_file)

$(OUT_PREFIX)_raymeta_contigs.fa: raymeta/Contigs.fasta
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

$(OUT_PREFIX)_fermi_contigs.fa: fermi/fmdef.p4.fa.gz
	gunzip -c $^ > $@

#*************************************************************************
#Abyss
#*************************************************************************
abyss/abyssk$(abyss_kmer)-contigs.fa: $(INPUT_PAIRED_END)
	mkdir -p abyss
	cd abyss && abyss-pe k=$(abyss_kmer) name='abyssk$(abyss_kmer)' in='../$^' np=16 2> $(log_file)

$(OUT_PREFIX)_abyss_contigs.fa: abyss/abyssk$(abyss_kmer)-contigs.fa
	ln -s $^ $@

#*************************************************************************
#MaSuRCA
#*************************************************************************
paired_ends/%_R1.fq : $(read_folder)/%_pe.fq
	mkdir -p $(dir $@)
	seqtk seq -1 $^ > $(firstword $@)

paired_ends/%_R2.fq: $(read_folder)/%_pe.fq
	mkdir -p $(dir $@)
	seqtk seq -2 $^ > $(lastword $@)

masurca/masurca.cfg: paired_ends/$(IN_PREFIX)_R1.fq paired_ends/$(IN_PREFIX)_R2.fq $(INPUT_SINGLE_END)
	mkdir -p $(dir $@)
	masurca -g masurca.cfg.1
	sed -e "/^PE=/cPE=pe $(masurca_pe_stats)  $(addprefix ../,$(wordlist 1,2,$^))" \
		-e "/^JUMP/c\PE=se $(masurca_se_stats)  ../$(lastword $^)" \
		-e "/^OTHER/d" \
		-e "/^USE_LINKING_MATES/s/0/1/" \
		-e  "/^GRAPH_KMER_SIZE/s/auto/$(masurca_kmer)/" \
		-e  "/^NUM_THREADS/s/\$$NUM_THREADS/$(threads)/" \
		masurca.cfg.1 > $@
		rm masurca.cfg.1

masurca/CA/10-gapclose/genome.ctg.fasta: masurca/masurca.cfg
	cd $(dir $@) && masurca $(notdir $<)
	cd $(dir $@) && bash assemble.sh

$(OUT_PREFIX)_masurca_contigs.fa: masurca/CA/10-gapclose/genome.ctg.fasta
	ln -s $^ $@
#*************************************************************************
#SGA quality filtering steps
#*************************************************************************
#1) Preprocess the data to remove ambiguous basecalls
#pe-mode=2 : interleaved
sga/$(sample_name)_sga.fq: $(INPUT_PAIRED_END)
	mkdir -p sga/
	$(SGA_BIN) preprocess --pe-mode 2 -o $@ $^ 2>> $(log_file)

#2) Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) --no-reverse $(notdir $^) 2>> $(log_file)

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
%.ec.fq : %_sga.fq %_sga.sai
	cd $(dir $@) && $(SGA_BIN) correct -k $(sga_ec_kmer) --discard --learn -t $(threads) -o $(notdir $@) $(notdir $<) 2>> $(log_file)

# Index the corrected data
%.ec.sai: %.ec.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) $(notdir $^) 2>> $(log_file)

# Remove exact-match duplicates and reads with low-frequency k-mers
%.ec.filter.pass.fa: %.ec.fq %.ec.sai
	cd $(dir $@) && $(SGA_BIN) filter -x $(sga_cov_filter) -t $(threads) --homopolymer-check \
		--low-complexity-check $(notdir $<) 2>> $(log_file)

# Compute the structure of the string graph
%.ec.filter.pass.asqg.gz: %.ec.filter.pass.fa
	cd $(dir $@) && $(SGA_BIN) overlap -m $(sga_min_overlap) -t $(threads) $(notdir $^) 2>> $(log_file)

# Perform the contig assembly
sga/%-contigs.fa: sga/%.ec.filter.pass.asqg.gz
	cd $(dir $@) && $(SGA_BIN) assemble -m $(sga_assemble_overlap) --min-branch-length $(TRIM_LENGTH) -o $* $(notdir $^) 2>> $(log_file)

$(OUT_PREFIX)_sga_contigs.fa: sga/$(sample_name)-contigs.fa
	ln $^ $@

#*************************************************************************
#Extract contigs > 500 bp
#*************************************************************************
%_ctgs_filt.fa : %_contigs.fa
	seqtk seq -L 500 $^ > $@ 2> $(log_file)

.PHONY: clean
clean:
	-rm -r sga/
	-rm -r abyss/
	-rm -r raymeta/
	-rm -r fermi
	-rm *.contigs.fa *.ctgs_filt.fa
