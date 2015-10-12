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

%_raw_metaphlan.txt: $(TMP_DIR)/reads_raw.fq
	$(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
		--bowtie2out $(TMP_DIR)/raw.bowtie2.bz2 --nproc $(threads) --input_type multifastq \
		--sample_id_key $(sample_name) --biom $*_raw_metaphlan.biom $< $@
	-rm $(TMP_DIR)/raw.bowtie2.bz2

#******************************************************************
# Process reads after contamination removal
#******************************************************************
$(TMP_DIR)/reads_nohuman.fq : $(wildcard $(read_folder)/*_?e.fq)
	cat $^ >> $@

%_nohuman_metaphlan.txt: $(TMP_DIR)/reads_nohuman.fq
	$(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
		--bowtie2out $(TMP_DIR )/nohuman.bowtie2.bz2 --nproc $(threads) --input_type multifastq \
		--sample_id_key $(sample_name) --biom $*_nohuman_metaphlan.biom $< $(word 1,$@)
	-rm $(TMP_DIR)/nohuman.bowtie2.bz2

#******************************************************************
# Classification of contigs
#******************************************************************
%_asm_metaphlan.txt: $(wildcard $(read_folder)/*_allctgs.fa)
	$(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
		--bowtie2out $(TMP_DIR)/asm.bowtie2.bz2 --nproc $(threads) --input_type multifasta \
		--sample_id_key $(sample_name) --biom $*_asm_metaphlan.biom $< $@
	-rm $(TMP_DIR)/asm.bowtie2.bz2
