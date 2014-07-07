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

ifndef read_folder
$(error Variable 'read_folder' is not defined)
endif

ifndef ctg_folder
$(error Variable 'ctg_folder' is not defined)
endif

ifndef prev_steps
prev_steps := qf_rmcont_asm
$(info 'prev_steps' is assumed to be $(prev_steps))
endif

ifndef step
step:= tax
$(warning Variable 'step' has been defined as '$(step)')
endif

#Run params
threads:=16

# ASSEMBLERS:= masurca raymeta fermi sga #abyss
ASSEMBLERS:= masurca

#Input and Output file prefixes
IN_PREFIX := $(sample_name)_$(prev_steps)
OUT_PREFIX:= $(IN_PREFIX)_$(step)

#Reads
READS_PAIRED_END := $(read_folder)/$(IN_PREFIX)_pe.fq
READS_SINGLE_END := $(read_folder)/$(IN_PREFIX)_se.fq

#Logging
log_name := $(CURDIR)/$(OUT_PREFIX)_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Databases
kraken_db := /labcommon/db/krakendb/kraken140311
refseq_virus_fna := /labcommon/db/fasta/refseq_viral/viral.1.1.genomic.fna
refseq_virus_faa := /labcommon/db/fasta/refseq_viral/viral.1.protein.faa
swissprot_faa := /labcommon/db/fasta/uniprot_sprot.fasta
tax_dmp_nucl := /labcommon/db/taxdb/gi_taxid_nucl.dmp
tax_dmp_prot := /labcommon/db/taxdb/gi_taxid_prot.dmp

#Blast parameters
blast_params:= -evalue 1 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

produce_outfiles = $(addsuffix _$(2),$(addprefix $(1)/$(IN_PREFIX)_,$(3)))

#Binary paths
FGS_PATH:= /labcommon/tools/FragGeneScan1.18


#Delete produced files if step fails
.DELETE_ON_ERROR:
#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all kraken_reports phmmer_vir

# all: kraken_reports phmmer_out blastp_out blastx_out
all: kraken_reports phmmer_vir phmmer_sprot blastn_vir blastp_vir blastp_sprot

#Outputs

kraken_reports: $(call produce_outfiles,kraken/,kraken.report,$(ASSEMBLERS))

phmmer_vir : $(call produce_outfiles,phmmer,fgs_phmmer_refseqvir.tbl,$(ASSEMBLERS))
phmmer_sprot : $(call produce_outfiles,phmmer,fgs_phmmer_sprot.tbl,$(ASSEMBLERS))

blastn_vir : $(call produce_outfiles,blastn,blastn_refseqvir.xml,$(ASSEMBLERS))

blastp_vir : $(call produce_outfiles,blastp,fgs_blastp_refseqvir.xml,$(ASSEMBLERS))
blastp_sprot : $(call produce_outfiles,blastp,fgs_blastp_sprot.xml,$(ASSEMBLERS))

blastx_vir : $(call produce_outfiles,blastx,blastx_refseqvir.xml,$(ASSEMBLERS))
blastx_sprot : $(call produce_outfiles,blastx,blastx_sprot.xml,$(ASSEMBLERS))

#*************************************************************************
#Call to Kraken - Salzberg
#*************************************************************************
#Other flags: --fastq-input
kraken/$(IN_PREFIX)_%_kraken.out: $(ctg_folder)/$(IN_PREFIX)_%_ctgs_filt.fa
	mkdir -p kraken
	@echo -e "\nClassifying $* contigs with Kraken\n\n" > $(log_file)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2>> $(log_file)

%_kraken.report: %_kraken.out
	@echo -e "\nCreating Kraken report for $* \n\n" > $(log_file)
	kraken-report --db $(kraken_db) $^ > $@ 2>> $(log_file)

#*************************************************************************
#FragGeneScan
#*************************************************************************
#Contig analysis
fgs/%_fgs.faa: $(ctg_folder)/%_ctgs_filt.fa
	mkdir -p $(dir $@)
	$(FGS_PATH)/run_FragGeneScan.pl -genome=$^ -out=$(basename $@) -complete=0 -train=illumina_10 2>> $(log_file)

# #Read analysis
# fgs/reads_%_fgs.faa : $(sam2fq_pre)_%.fa
# 	mkdir -p $(dir $@)
# 	$(FGS_PATH)/run_FragGeneScan.pl -genome=$< -out=$(basename $@) -complete=0 -train=illumina_10 2>> $(log_file)

#All filtered reads analysis
#Convert fq to fa
$(sam2fq_pre)_%.fa: $(sam2fq_pre)_%.fq
	seqtk seq -A $^ > $@ 2> $(log_file)
	#sed -n '1~4s/^@/>/p;2~4p' $^

#*************************************************************************
#Phmmer
#*************************************************************************
#Optional --domtblout $(basename $@).dom

#Contigs against swissprot
phmmer/%_fgs_phmmer_sprot.tbl : fgs/%_fgs.faa $(swissprot_faa)
	mkdir -p phmmer
	phmmer --cpu $(threads) --noali --tblout $@ $^ > /dev/null 2>> $(log_file)

#Contigs against refseq virus proteins
phmmer/%_fgs_phmmer_refseqvir.tbl : fgs/%_fgs.faa $(refseq_virus_faa)
	mkdir -p phmmer
	phmmer --cpu $(threads) --noali --tblout $@ $^ > /dev/null 2>> $(log_file)

#*************************************************************************
#BlastN - Nucleotides
#*************************************************************************
#Create blastdb from RefSeq virus both nucleotide and protein
refseq_virus_fna_blastdb: $(refseq_virus_fna)
	mkdir -p $@
	cd $@ && makeblastdb -dbtype nucl -out $@ -title $@ -parse_seqids -taxid_map $(tax_dmp_nucl) -in $^

#Blastn against refseq virus
blastn/%_blastn_refseqvir.xml: $(ctg_folder)/%_ctgs_filt.fa | refseq_virus_fna_blastdb
	mkdir -p $(dir $@)
	blastn -task blastn $(blast_params) $(blastn_params) -db $|/$| -query $^ -out $@ 2>> $(log_file)

#*************************************************************************
#BlastX - Proteins
#*************************************************************************
refseq_virus_faa_blastdb: $(refseq_virus_faa)
	mkdir -p $@
	cd $@ && makeblastdb -dbtype prot -out $@ -title $@ -parse_seqids -taxid_map $(tax_dmp_prot) -in $^

blastx/%_blastx_sprot.xml : $(ctg_folder)/%_ctgs_filt.fa
	mkdir -p $(dir $@)
	blastx $(blast_params) -db $(blastdb_folder)/nr/swissprot -query $< -out $@ 2>> $(log_file)

blastx/%_blastx_refseqvir.xml : $(ctg_folder)/%_ctgs_filt.fa
	mkdir -p $(dir $@)
	blastx $(blast_params) -db $|/$| -query $< -out $@ 2>> $(log_file)

#BlastP the predicted ORF to swissprot
blastp/%_fgs_blastp.xml: fgs/%_fgs.faa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_folder)/nr/swissprot -query $< -out $@ 2>> $(log_file)
