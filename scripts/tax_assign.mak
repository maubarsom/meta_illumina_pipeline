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
$(warning Variable 'read_folder' will be assumed to be "./")
read_folder := ./
endif

ifndef prev_steps
prev_steps := qf_rmcont_asm
$(info 'prev_steps' is assumed to be $(prev_steps))
endif

ifndef step
step:= tax
$(warning Variable 'step' has been defined as '$(step)')
endif

ASSEMBLERS:=raymeta fermi sga abyss

#Input and Output file prefixes
IN_PREFIX := $(sample_name)_$(prev_steps)
OUT_PREFIX:= $(IN_PREFIX)_asm

#Reads
INPUT_PAIRED_END := $(read_folder)/$(IN_PREFIX)_pe.fq
INPUT_SINGLE_END := $(read_folder)/$(IN_PREFIX)_se.fq


#*************************************************************************
#Call to Kraken - Salzberg
#*************************************************************************

#Other flags: --fastq-input
kraken_raymeta.out kraken_fermi.out kraken_abyss.out: kraken_%.out : %_contigs_1k.fa
	@echo -e "\nClassifying $* contigs with Kraken\n\n" > $(log_file)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2> $(log_file)

kraken_raymeta.report kraken_fermi.report kraken_abyss.report: kraken_%.report : kraken_%.out
	@echo -e "\nCreating Kraken report for $* \n\n" > $(log_file)
	kraken-report --db $(kraken_db) $^ > $@ 2> $(log_file)

#*************************************************************************
#FragGeneScan
#*************************************************************************
#Contig analysis
raymeta_fgs.faa fermi_fgs.faa abyss_fgs.faa: %_fgs.faa : %_contigs_1k.fa
	$(FragGeneScan_folder)/run_FragGeneScan.pl -genome=$^ -out=$*_fgs -complete=0 -train=illumina_10 2> $(log_file)

#All filtered reads analysis
#Convert fq to fa
$(sam2fq_pre)_%.fa: $(sam2fq_pre)_%.fq
	seqtk seq -A $^ > $@ 2> $(log_file)
	#sed -n '1~4s/^@/>/p;2~4p' $^

reads_pe_fgs.faa reads_se_fgs.faa: reads_%_fgs.faa : $(sam2fq_pre)_%.fa
	$(FragGeneScan_bin) -genome=$< -out=$(basename $@) -complete=0 -train=illumina_10 2> $(log_file)


#*************************************************************************
#Phmmer
#*************************************************************************
#Optional --domtblout $(basename $@).dom
#Contigs
raymeta_fgs_phmmer.tbl fermi_fgs_phmmer.tbl abyss_fgs_phmmer.tbl reads_pe_fgs_phmmer.tbl reads_se_fgs_phmmer.tbl: %_fgs_phmmer.tbl : %_fgs.faa $(swissprot_fasta)
	phmmer --cpu $(threads) --noali --tblout $@ $^ >/dev/null 2> $(log_file)

#*************************************************************************
#Blast+
#*************************************************************************
#Create index from RefSeq virus
refseq_virus: $(refseq_virus_fasta)
	mkdir -p $@
	cd $@ && makeblastdb -dbtype nucl -out $@ -title $@ -parse_seqids -taxid_map /labcommon/db/taxdb/gi_taxid_nucl.dmp -in $^

#Blastn against refseq virus
raymeta_blastn_vir.xml fermi_blastn_vir.xml abyss_blastn_vir.xml: %_blastn_vir.xml: %_contigs_1k.fa | refseq_virus
	blastn -task blastn $(blast_params) $(blastn_params) -query $^ -db refseq_virus/refseq_virus > $@

reads_pe_blastn_vir.xml reads_se_blastn_vir.xml: reads_%_blastn_vir.xml: $(sam2fq_pre)_%.fa | refseq_virus
	blastn -task blastn $(blast_params) $(blastn_params) -query $^ -db refseq_virus/refseq_virus > $@

#Blastx filtered contigs
#Uses differnt genetic code
raymeta_blastx.xml fermi_blastx.xml abyss_blastx.xml: %_blastx.xml : %_contigs_1k.fa
	blastx -query $< -db $(blastdb_folder)/nr/swissprot -query_genetic_code 11 -out $@ $(blast_params) 2> $(log_file)

#Blastp fgs proteins
#Reads and contigs , singletons only pending
raymeta_fgs_blastp.xml fermi_fgs_blastp.xml abyss_fgs_blastp.xml reads_pe_fgs_blastp.xml reads_se_fgs_blastp.xml: %_fgs_blastp.xml: %_fgs.faa
	blastp $(blast_params) -db $(blastdb_folder)/nr/swissprot -query $< -out $@  2> $(log_file)
