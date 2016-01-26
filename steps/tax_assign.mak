# Homology search pipeline

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

ifndef in_folder
$(error Variable 'in_folder' is not defined)
endif

ifndef in_steps
in_steps := qf_rmcont_asm
$(info 'in_steps' is assumed to be $(in_steps))
endif

ifndef step
step:= tax
$(warning Variable 'step' has been defined as '$(step)')
endif

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

#Input and Output file prefixes
IN_CTG_PREFIX := $(sample_name)_$(ctg_steps)
IN_READ_PREFIX := $(sample_name)_$(read_steps)
OUT_PREFIX:= $(IN_CTG_PREFIX)_$(step)

#Blast parameters
blast_params:= -evalue 1 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

ctg_outfile = $(1)/$(IN_CTG_PREFIX)_allctgs_$(2)
read_outfiles = $(addsuffix _$(2),$(addprefix $(1)/$(IN_READ_PREFIX)_,$(3)))

#Delete produced files if step fails
.DELETE_ON_ERROR:
#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.INTERMEDIATE: $(TMP_DIR)/%.fa

.PHONY: all
.PHONY: kraken_reports blastn_vir blastn_nt
.PHONY: blastp_vir blastp_nr blastp_sprot
.PHONY: blastx_vir blastx_nr blastx_sprot
.PHONY: hmmscan_pfam hmmscan_vfam phmmer_vir phmmer_sprot

all: diamond_nr
# all: hmmscan_pfam hmmscan_vfam
# all: phmmer_vir

#Outputs
diamond_nr : $(call ctg_outfile,diamond,diamond_nr.sam.gz.md5)
diamond_nr : $(call read_outfiles,diamond,diamond_nr.sam.gz.md5,pe se)

phmmer_vir : $(call ctg_outfile,phmmer,fgs_phmmer_refseqvir.tbl)
phmmer_sprot : $(call ctg_outfile,phmmer,fgs_phmmer_sprot.tbl)

hmmscan_pfam : $(call ctg_outfile,hmmscan,fgs_hmmscan_pfam.tbl)
hmmscan_pfam : $(call read_outfiles,hmmscan,fgs_hmmscan_pfam.tbl,pe se)

hmmscan_vfam : $(call ctg_outfile,hmmscan,fgs_hmmscan_vfam.tbl)
hmmscan_vfam : $(call read_outfiles,hmmscan,fgs_hmmscan_vfam.tbl,pe se)

#*************************************************************************
#Diamond (Tubingen) - Proteins
#*************************************************************************
#Contigs to NR
#Can add --sensitive flag for slower but more accurate results
#--seg yes/no for low complexity masking
diamond/%_diamond_nr.daa : $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(DIAMOND_BIN) blastx --sensitive -p $(threads) --db $(diamond_nr) --query $< --daa $@ --tmpdir $(TMP_DIR) --seg yes

diamond/%_diamond_nr.daa : $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(DIAMOND_BIN) blastx --sensitive -p $(threads) --db $(diamond_nr) --query $< --daa $@ --tmpdir $(TMP_DIR) --seg yes

diamond/%.sam.gz : diamond/%.daa
	$(DIAMOND_BIN) view --daa $^ --out $(dir $@)/$*.sam --outfmt sam
	gzip $(dir $@)/$*.sam

$(TMP_DIR)/%_diamond_nohits.fa: diamond/%_diamond_nr.sam.gz ../assembly/%.fa
	python ../scripts/extract_diamond_nohits.py $< $(word 2)

#*************************************************************************
#EMBOSS ORF prediction
#*************************************************************************
#Contig analysis
orf/%_emboss_find0.faa: $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(FGS_BIN) -genome=$^ -out=$(basename $@) -complete=0 -train=illumina_10

#Read analysis
orf/%_emboss_find1.faa : $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(FGS_BIN) -genome=$< -out=$(basename $@) -complete=0 -train=illumina_10

#*************************************************************************
#FragGeneScan
#*************************************************************************
#Contig analysis
fgs/%_fgs.faa: $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(FGS_BIN) -genome=$^ -out=$(basename $@) -complete=0 -train=illumina_10

#Read analysis
fgs/%_fgs.faa : $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(FGS_BIN) -genome=$< -out=$(basename $@) -complete=0 -train=illumina_10

#*************************************************************************
#HMMSCAN
#*************************************************************************
#Optional --domtblout $(basename $@).dom

#Contigs against pfam
hmmscan/%_fgs_hmmscan_pfam.tbl : $(pfam_hmm_db) fgs/%_fgs.faa
	mkdir -p $(dir $@)
	$(HMMSCAN_BIN) --cpu $(threads) --noali --tblout $@ --domtblout $(basename $@)_dom.tbl $^ > /dev/null

#Contigs against vFam (Skewes-Cox,2014)
hmmscan/%_fgs_hmmscan_vfam.tbl : $(vfam_hmm_db) fgs/%_fgs.faa
	mkdir -p $(dir $@)
	$(HMMSCAN_BIN) --cpu $(threads) --noali --tblout $@ --domtblout $(basename $@)_dom.tbl $^ > /dev/null

#*************************************************************************
# Calculate checksums
#*************************************************************************
%.md5: %
	md5sum $< > $@

#*************************************************************************
#CLEANING RULES
#*************************************************************************
.PHONY: clean-tmp

clean-tmp:
	-rm $(TMP_DIR)/*.fa
