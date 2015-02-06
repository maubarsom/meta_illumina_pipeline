SHELL := /bin/bash

#External parameters
# 1) basename
# 2) read_folder

ifndef sample_name
$(error Variable 'sample_name' is not defined)
endif

ifndef read_folder
$(error Variable 'read_folder' is not defined)
endif

ifndef threads
threads := $(shell nproc)
endif

ifndef TMP_DIR
$(shell mkdir -p tmp)
TMP_DIR := ./tmp
endif

.PHONY: raw filtered asm

raw: $(sample_name)_raw_metaphlan.txt $(sample_name)_raw_metaphlan.biom
filtered: $(sample_name)_nohuman_metaphlan.txt $(sample_name)_nohuman_metaphlan.biom
asm: $(sample_name)_asm_metaphlan.txt $(sample_name)_asm_metaphlan.biom

#******************************************************************
# Raw reads
#******************************************************************
$(TMP_DIR)/reads_raw.fq : $(wildcard $(read_folder)/*.fastq.gz)
	zcat $^ >> $@

%_raw_metaphlan.txt %_raw_metaphlan.biom: $(TMP_DIR)/reads_raw.fq
	mkdir -p metaphlan
	$(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
		--bowtie2out $(TMP_DIR)/raw.bowtie2.bz2 --nproc $(threads) --input_type multifastq --biom $(word 2,$@) \
		$< $(word 1,$@)

#******************************************************************
# Process reads after contamination removal
#******************************************************************
$(TMP_DIR)/reads_filt.fq : $(wildcard $(read_folder)/*_?e.fq)
	cat $^ >> $@

%_filt_metaphlan.txt %_filt_metaphlan.biom: $(TMP_DIR)/reads_filt.fq
	mkdir -p metaphlan
	$(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
		--bowtie2out $(TMP_DIR)/filt.bowtie2.bz2 --nproc $(threads) --input_type multifastq --biom $(word 2,$@) \
		$< $(word 1,$@)

#******************************************************************
# Classification of contigs
#******************************************************************
%_asm_metaphlan.txt %_asm_metaphlan.biom: $(wildcard $(read_folder)/*_allctgs.fa)
	mkdir -p metaphlan
	$(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
		--bowtie2out $(TMP_DIR)/asm.bowtie2.bz2 --nproc $(threads) --input_type multifasta --biom $(word 2,$@) \
		$< $(word 1,$@)
