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

ifndef read_steps
read_steps := qf_rmcont
$(info 'read_steps' is assumed to be $(read_steps))
endif

ifndef ctg_steps
ctg_steps := qf_rmcont_asm
$(info 'ctg_steps' is assumed to be $(ctg_steps))
endif

ifndef step
step:= tax
$(warning Variable 'step' has been defined as '$(step)')
endif

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

ifndef ASSEMBLERS
	$(error Define ASSEMBLERS variable in make.cfg file)
endif

#Input and Output file prefixes
IN_CTG_PREFIX := $(sample_name)_$(ctg_steps)
IN_READ_PREFIX := $(sample_name)_$(read_steps)
OUT_PREFIX:= $(IN_CTG_PREFIX)_$(step)

#Reads
READS_PAIRED_END := $(read_folder)/$(IN_READ_PREFIX)_pe.fq
READS_SINGLE_END := $(read_folder)/$(IN_READ_PREFIX)_se.fq

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
blastdb_folder:=/labcommon/db/blastdb/
pfam_hmm_db := /labcommon/db/hmmerdb/Pfam-A.hmm
vfam_hmm_db := /labcommon/db/hmmerdb/vFam-A_2014.hmm

#Blast parameters
blast_params:= -evalue 1 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

ctg_outfiles = $(addsuffix _$(2),$(addprefix $(1)/$(IN_CTG_PREFIX)_,$(3)))
read_outfiles = $(addsuffix _$(2),$(addprefix $(1)/$(IN_READ_PREFIX)_,$(3)))

#Binary paths
FGS_PATH:= /labcommon/tools/FragGeneScan1.18


#Delete produced files if step fails
.DELETE_ON_ERROR:
#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all
.PHONY: kraken_reports blastn_vir blastn_nt
.PHONY: blastp_vir blastp_nr blastp_sprot
.PHONY: blastx_vir blastx_nr blastx_sprot
.PHONY: hmmscan_pfam hmmscan_vfam phmmer_vir phmmer_sprot

all: kraken_reports blastn_vir blastn_nt
#all: blastp_vir blastp_nr
all: blastx_vir blastx_nr
all: hmmscan_pfam hmmscan_vfam phmmer_vir

#Outputs

kraken_reports: $(call ctg_outfiles,kraken,kraken.report,$(ASSEMBLERS))

phmmer_vir : $(call ctg_outfiles,phmmer,fgs_phmmer_refseqvir.tbl,$(ASSEMBLERS))
phmmer_sprot : $(call ctg_outfiles,phmmer,fgs_phmmer_sprot.tbl,$(ASSEMBLERS))

hmmscan_pfam : $(call ctg_outfiles,hmmscan,fgs_hmmscan_pfam.tbl,$(ASSEMBLERS))
hmmscan_vfam : $(call ctg_outfiles,hmmscan,fgs_hmmscan_vfam.tbl,$(ASSEMBLERS))

blastn_nt : $(call ctg_outfiles,blastn,blastn_nt.xml,$(ASSEMBLERS))
blastn_vir : $(call ctg_outfiles,blastn,blastn_refseqvir.xml,$(ASSEMBLERS))
blastn_vir : $(call read_outfiles,blastn,blastn_refseqvir.xml,pe se)

blastp_vir : $(call ctg_outfiles,blastp,fgs_blastp_refseqvir.xml,$(ASSEMBLERS))
blastp_vir : $(call read_outfiles,blastp,fgs_blastp_refseqvir.xml,pe se)
blastp_nr : $(call ctg_outfiles,blastp,fgs_blastp_nr.xml,$(ASSEMBLERS))
blastp_sprot : $(call ctg_outfiles,blastp,fgs_blastp_sprot.xml,$(ASSEMBLERS))

blastx_nr : $(call ctg_outfiles,blastx,blastx_nr.xml,$(ASSEMBLERS))
blastx_sprot : $(call ctg_outfiles,blastx,blastx_sprot.xml,$(ASSEMBLERS))
blastx_vir : $(call ctg_outfiles,blastx,blastx_refseqvir.xml,$(ASSEMBLERS))
blastx_vir : $(call read_outfiles,blastx,blastx_refseqvir.xml,pe se)

#*************************************************************************
#Convert Reads from Fastq to Fasta
#*************************************************************************
#Convert fq to fa
reads_fa/%.fa: $(read_folder)/%.fq
	mkdir -p $(dir $@)
	seqtk seq -A $^ > $@ 2> $(log_file)

#*************************************************************************
#Call to Kraken - Salzberg
#*************************************************************************
#Other flags: --fastq-input
kraken/$(IN_CTG_PREFIX)_%_kraken.out: $(ctg_folder)/$(IN_CTG_PREFIX)_%_ctgs_filt.fa
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

#Read analysis
fgs/%_fgs.faa : reads_fa/%.fa
	mkdir -p $(dir $@)
	$(FGS_PATH)/run_FragGeneScan.pl -genome=$< -out=$(basename $@) -complete=0 -train=illumina_10 2>> $(log_file)

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
#HMMSCAN
#*************************************************************************
#Optional --domtblout $(basename $@).dom

#Contigs against pfam
hmmscan/%_fgs_hmmscan_pfam.tbl : $(pfam_hmm_db) fgs/%_fgs.faa
	mkdir -p $(dir $@)
	hmmscan --cpu $(threads) --noali --tblout $@ $^ > /dev/null 2>> $(log_file)

#Contigs against vFam (Skewes-Cox,2014)
hmmscan/%_fgs_hmmscan_vfam.tbl : $(vfam_hmm_db) fgs/%_fgs.faa
	mkdir -p $(dir $@)
	hmmscan --cpu $(threads) --noali --tblout $@ $^ > /dev/null 2>> $(log_file)

#*************************************************************************
#BlastN - Nucleotides
#*************************************************************************
#Create blastdb from RefSeq virus both nucleotide and protein
refseq_virus_fna_blastdb: $(refseq_virus_fna)
	mkdir -p $@
	cd $@ && makeblastdb -dbtype nucl -out $@ -title $@ -parse_seqids -taxid_map $(tax_dmp_nucl) -in $^

#Reads to Refseq Viral nucleotides
blastn/%_blastn_refseqvir.xml: reads_fa/%.fa | refseq_virus_fna_blastdb
	mkdir -p $(dir $@)
	blastn -task blastn $(blast_params) $(blastn_params) -db $|/$| -query $^ -out $@ 2>> $(log_file)

#Contigs to Refseq Viral nucleotides
blastn/%_blastn_refseqvir.xml: $(ctg_folder)/%_ctgs_filt.fa | refseq_virus_fna_blastdb
	mkdir -p $(dir $@)
	blastn -task blastn $(blast_params) $(blastn_params) -db $|/$| -query $^ -out $@ 2>> $(log_file)

#Contigs to nt
blastn/%_blastn_nt.xml: $(ctg_folder)/%_ctgs_filt.fa
	mkdir -p $(dir $@)
	blastn -task blastn $(blast_params) $(blastn_params) -db $(blastdb_folder)/nt/nt -query $^ -out $@ 2>> $(log_file)

#*************************************************************************
#BlastX - Proteins
#*************************************************************************
#Build blastdb from refseq viral proteins
refseq_virus_faa_blastdb: $(refseq_virus_faa)
	mkdir -p $@
	cd $@ && makeblastdb -dbtype prot -out $@ -title $@ -parse_seqids -taxid_map $(tax_dmp_prot) -in $^

#Contigs to NR
blastx/%_blastx_nr.xml : $(ctg_folder)/%_ctgs_filt.fa
	mkdir -p $(dir $@)
	blastx $(blast_params) -db $(blastdb_folder)/nr/nr -query $< -out $@ 2>> $(log_file)

#Contigs to Swissprot
blastx/%_blastx_sprot.xml : $(ctg_folder)/%_ctgs_filt.fa
	mkdir -p $(dir $@)
	blastx $(blast_params) -db $(blastdb_folder)/nr/swissprot -query $< -out $@ 2>> $(log_file)

#Contigs to Refseq Virus Proteins
blastx/%_blastx_refseqvir.xml : $(ctg_folder)/%_ctgs_filt.fa | refseq_virus_faa_blastdb
	mkdir -p $(dir $@)
	blastx $(blast_params) -db $|/$| -query $< -out $@ 2>> $(log_file)

#Reads to Refseq Virus Proteins
blastx/%_blastx_refseqvir.xml : reads_fa/%.fa | refseq_virus_faa_blastdb
	mkdir -p $(dir $@)
	blastx $(blast_params) -db $|/$| -query $< -out $@ 2>> $(log_file)

#*************************************************************************
#BlastP - Predicted ORF to Proteins
#*************************************************************************
blastp/%_fgs_blastp_nr.xml: fgs/%_fgs.faa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_folder)/nr/nr -query $< -out $@ 2>> $(log_file)

blastp/%_fgs_blastp_sprot.xml: fgs/%_fgs.faa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_folder)/nr/swissprot -query $< -out $@ 2>> $(log_file)

blastp/%_fgs_blastp_refseqvir.xml: fgs/%_fgs.faa | refseq_virus_faa_blastdb
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $|/$| -query $< -out $@ 2>> $(log_file)
