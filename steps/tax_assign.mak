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

ifndef step
step:= tax
$(warning Variable 'step' has been defined as '$(step)')
endif

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

#Input and Output file prefixes
INPUT_CONTIG := $(sample_name)_contigs
INPUT_PE := $(sample_name)_pe
INPUT_SINGLE := $(sample_name)_single
INPUT_MERGED := $(sample_name)_merged

OUT_PREFIX:= $(IN_CTG_PREFIX)_$(step)

#Blast parameters
blast_params:= -evalue 1 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

generate_outfiles = $(addsuffix _$(2),$(addprefix $(1)/$(sample_name)_,$(3)))

#Delete produced files if step fails
.DELETE_ON_ERROR:
#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.INTERMEDIATE: $(TMP_DIR)/%.fa

.PHONY: all
.PHONY: hmmscan_pfam hmmscan_vfam
.PHONY: phmmer_vir phmmer_sprot

all: megan #diamond_nr
# all: hmmscan_pfam hmmscan_vfam
# all: phmmer_vir phmmer_sprot

#Outputs
diamond_nr :   $(call generate_outfiles,diamond,diamond_nr.sam.gz.md5,contigs pe single merged)

megan: megan/$(sample_name)_diamond_nr.rma

phmmer_vir :   $(call generate_outfiles,phmmer,fgs_phmmer_refseqvir.tbl,contigs)

phmmer_sprot : $(call generate_outfiles,phmmer,fgs_phmmer_sprot.tbl,contigs )

hmmscan_pfam : $(call generate_outfiles,hmmscan,fgs_hmmscan_pfam.tbl,contigs pe single merged)

hmmscan_vfam : $(call generate_outfiles,hmmscan,fgs_hmmscan_vfam.tbl,contigs pe single merged)

#*************************************************************************
#Diamond (Tubingen) - Proteins
#*************************************************************************
#Contigs to NR
#Can add --sensitive flag for slower but more accurate results
#--seg yes/no for low complexity masking
diamond/%_diamond_nr.daa : $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(DIAMOND_BIN) blastx --sensitive -p $(threads) --db $(diamond_nr) --query $< --daa $@ --tmpdir $(TMP_DIR) --seg yes

diamond/%.sam.gz : diamond/%.daa
	$(DIAMOND_BIN) view --daa $^ --out $(dir $@)/$*.sam --outfmt sam
	gzip $(dir $@)/$*.sam

$(TMP_DIR)/%_diamond_nohits.fa: diamond/%_diamond_nr.sam.gz ../assembly/%.fa
	python ../scripts/extract_diamond_nohits.py $< $(word 2)

#*************************************************************************
# MEGAN 5 parsing of diamond results
#*************************************************************************
#The script assumes files with contigs pe single and merged.sam.gz
megan/$(sample_name)_diamond_nr.rma: $(call generate_outfiles,diamond,diamond_nr.sam.gz,contigs pe single merged) | ../scripts/megan_diamond_nr.m4
	mkdir -p $(dir $@)
	m4 -DTAX2GI=$(megan_gi2tax) -DPREFIX=$(sample_name) -DOUT_FILE=$@ $| > $(TMP_DIR)/megan_script.txt
	$(MEGAN_BIN) -g -c $(TMP_DIR)/megan_script.txt

#*************************************************************************
#FragGeneScan
#*************************************************************************
#Contig analysis
fgs/%_fgs.faa: $(in_folder)/%.fa
	mkdir -p $(dir $@)
	$(FGS_BIN) -genome=$^ -out=$(basename $@) -complete=0 -train=illumina_10

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
