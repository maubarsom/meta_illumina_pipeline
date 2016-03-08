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

ifndef step
step:= asm
$(warning Variable 'step' has been defined as '$(step)')
endif

#Input and Output file prefixes
OUT_PREFIX:= $(sample_name)_$(step)

#Reads
INPUT_PE := $(read_folder)/$(sample_name)_pe.fq
INPUT_SINGLE := $(read_folder)/$(sample_name)_single.fq
INPUT_MERGED := $(read_folder)/$(sample_name)_merged.fq

#Run params
ifndef threads
	$(error Define threads variable in make.cfg file)
endif

ifndef ASSEMBLERS
	$(error Define ASSEMBLERS variable in make.cfg file)
endif

#Creates a OUT_PREFIX_assembler_ctgs_filt.fa for each assembler
CTG_FILES:= $(addsuffix _ctgs_filt.fa,$(addprefix contigs_filt/$(OUT_PREFIX)_,$(ASSEMBLERS)))

#Ray
ray_kmer := 31

#Fermi parameters
fermi_overlap := 40

#Fermi pe parameters
fermi_pe_overlap := 40
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
masurca_pe_stats := 300 80
masurca_se_stats := 250 50

#Samtools parameters
#If paired-end exclude if read mapped in proper pair. If single-end include if unmapped
samtools_filter_flag = $(if $(filter pe,$*),-F2,-f4)

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.INTERMEDIATE: $(TMP_DIR)/singletons_pe.fq $(TMP_DIR)/singletons_single.fq $(TMP_DIR)/singletons_merged.fq $(TMP_DIR)/$(OUT_PREFIX)_contigs.fa.bwt

.PHONY: all

all: $(OUT_PREFIX)_contigs.fa $(OUT_PREFIX)_pe.fa $(OUT_PREFIX)_single.fa $(OUT_PREFIX)_merged.fa

#Concatenate the contigs from all assembly strategies into a single file
$(OUT_PREFIX)_contigs.fa: $(CTG_FILES)
	cat $^ > $@

$(OUT_PREFIX)_%.fa: singletons/singletons_%.fa
	ln -s $^ $@

#*************************************************************************
#MetaRay
#*************************************************************************
raymeta/Contigs.fasta: $(INPUT_PE) $(INPUT_SINGLE)
	@echo -e "\nAssembling reads with Ray Meta\n\n"
	mpiexec -n 16 Ray Meta -k $(ray_kmer) -i $< -s $(word 2,$^) -o $(dir $@)

#Ray Assembler duplicates contigs for some reason, probably due to MPI
contigs/$(OUT_PREFIX)_raymeta_contigs.fa: raymeta/Contigs.fasta
	../scripts/deduplicate_raymeta_ctgs.py -o $@ $^

#*************************************************************************
#Fermi - all
#*************************************************************************
#Runs Fermi assembler with both single and paired-ends assuming they are all single-ends
#Fermi outputs only until the 2nd step
#Qualities are not necessary for Kraken/Blast classification
fermi/fmdef.p2.mag.gz: $(INPUT_PE) $(INPUT_MERGED) $(INPUT_SINGLE)
	mkdir -p $(dir $@)
	cd $(dir $@) && run-fermi.pl -t $(threads) -k $(fermi_overlap) $(addprefix ../,$^) > assembly.mak
	cd $(dir $@) && $(MAKE) -f assembly.mak -j $(threads) $(notdir $@)

contigs/$(OUT_PREFIX)_fermi_contigs.fa: fermi/fmdef.p2.mag.gz
	mkdir -p $(dir $@)
	gunzip $^ |	$(SEQTK_BIN) seq -A - > $@

#*************************************************************************
#Megahit
#*************************************************************************
megahit/final.contigs.fa: $(INPUT_PE) $(INPUT_MERGED) $(INPUT_SINGLE)
	$(MEGAHIT_BIN) -m 5e10 -l $$(( 2*$(READ_LEN) )) --k-step 4 --k-max 81 --12 $< -r $(word 2,$^),$(word 3,$^) --cpu-only -t $(threads) -o megahit

contigs/$(OUT_PREFIX)_megahit_contigs.fa: megahit/final.contigs.fa
	mkdir -p $(dir $@)
	cp $^ $@

#*************************************************************************
#Spades
#*************************************************************************
spades/contigs.fasta: $(INPUT_PE) $(INPUT_SINGLE) $(INPUT_MERGED)
	mkdir -p $(dir $@)
	$(SPADES_BIN) -t $(threads) --pe1-12 $< --pe1-s $(word 2,$^) --s1 $(word 3,$^) -o $(dir $@) --tmp-dir $(TMP_DIR) -m 64

contigs/$(OUT_PREFIX)_spades_contigs.fa: spades/contigs.fasta
	mkdir -p $(dir $@)
	cp $^ $@

#*************************************************************************
#MaSuRCA -- OBSOLETE : Rules must be updated !
#*************************************************************************
paired_ends/%_R1.fq : $(read_folder)/%_pe.fq
	mkdir -p $(dir $@)
	$(SEQTK_BIN) seq -1 $^ > $(firstword $@)

paired_ends/%_R2.fq: $(read_folder)/%_pe.fq
	mkdir -p $(dir $@)
	$(SEQTK_BIN) seq -2 $^ > $(lastword $@)

masurca/masurca.cfg: paired_ends/$(IN_PREFIX)_R1.fq paired_ends/$(IN_PREFIX)_R2.fq $(INPUT_SINGLE)
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
	cd masurca && masurca $(notdir $<)
	cd masurca && bash assemble.sh

contigs/$(OUT_PREFIX)_masurca_contigs.fa: masurca/CA/10-gapclose/genome.ctg.fasta
	mkdir -p $(dir $@)
	cp $^ $@

#*************************************************************************
#SGA -- OBSOLETE : Rules must be updated !
#*************************************************************************
#1) Preprocess the data to remove ambiguous basecalls
#pe-mode=2 : interleaved
sga/$(sample_name)_sga.fq: $(INPUT_PE)
	mkdir -p sga/
	$(SGA_BIN) preprocess --pe-mode 2 -o $@ $^

#2) Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) --no-reverse $(notdir $^)

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
%.ec.fq : %_sga.fq %_sga.sai
	cd $(dir $@) && $(SGA_BIN) correct -k $(sga_ec_kmer) --discard --learn -t $(threads) -o $(notdir $@) $(notdir $<)

# Index the corrected data
%.ec.sai: %.ec.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) $(notdir $^)

# Remove exact-match duplicates and reads with low-frequency k-mers
%.ec.filter.pass.fa: %.ec.fq %.ec.sai
	cd $(dir $@) && $(SGA_BIN) filter -x $(sga_cov_filter) -t $(threads) --homopolymer-check \
		--low-complexity-check $(notdir $<)

# Compute the structure of the string graph
%.ec.filter.pass.asqg.gz: %.ec.filter.pass.fa
	cd $(dir $@) && $(SGA_BIN) overlap -m $(sga_min_overlap) -t $(threads) $(notdir $^)

# Perform the contig assembly
sga/%-contigs.fa: sga/%.ec.filter.pass.asqg.gz
	cd $(dir $@) && $(SGA_BIN) assemble -m $(sga_assemble_overlap) --min-branch-length $(TRIM_LENGTH) -o $* $(notdir $^)

contigs/$(OUT_PREFIX)_sga_contigs.fa: sga/$(sample_name)-contigs.fa
	mkdir -p $(dir $@)
	cp $^ $@

#*************************************************************************
#Extract contigs > 500 bp
#*************************************************************************
#Seqtk is used to sample sequences > 500bp
#Awk renames contigs to make them friendly for RAPSEARCH and other tools that do not support spaces in the names
contigs_filt/%_ctgs_filt.fa: contigs/%_contigs.fa
	mkdir -p $(dir $@)
	$(SEQTK_BIN) seq -L 500 $< | awk -vPREFIX=$*  'BEGIN {counter=1; split(PREFIX,fields,"_asm_"); ASM=fields[2];} /^>/ {print ">" ASM "_" counter ;counter+=1;} ! /^>/{print $0;}' > $@

#*************************************************************************
#Extract singletons
#*************************************************************************
$(TMP_DIR)/$(OUT_PREFIX)_contigs.fa.bwt: $(OUT_PREFIX)_contigs.fa
	$(BWA_BIN) index -p $(basename $@) $<

singletons/pe_to_contigs.bam: $(INPUT_PE) $(TMP_DIR)/$(OUT_PREFIX)_contigs.fa.bwt
	mkdir -p singletons
	$(BWA_BIN) mem -t $(threads) -T 30 -M -p $(basename $(word 2,$^)) $< | $(SAMTOOLS_BIN) view -F 256 -hSb -o $@ -

singletons/single_to_contigs.bam: $(INPUT_SINGLE) $(TMP_DIR)/$(OUT_PREFIX)_contigs.fa.bwt
	mkdir -p singletons
	$(BWA_BIN) mem -t $(threads) -T 30 -M    $(basename $(word 2,$^)) $< | $(SAMTOOLS_BIN) view -F 256 -hSb -o $@ -

singletons/merged_to_contigs.bam: $(INPUT_MERGED) $(TMP_DIR)/$(OUT_PREFIX)_contigs.fa.bwt
	mkdir -p singletons
	$(BWA_BIN) mem -t $(threads) -T 30 -M    $(basename $(word 2,$^)) $< | $(SAMTOOLS_BIN) view -F 256 -hSb -o $@ -

singletons/singletons_pe.fa singletons/singletons_single.fa singletons/singletons_merged.fa: $(TMP_DIR)/singletons_%.fa: singletons/%_to_contigs.bam
	$(SAMTOOLS_BIN) view $(samtools_filter_flag) -hb $< | $(PICARD_BIN) SamToFastq INPUT=/dev/stdin FASTQ=/dev/stdout INTERLEAVE=True | $(SEQTK_BIN) seq -A - > $@

.PHONY: clean
clean:
	-rm -r sga/
	-rm -r abyss/
	-rm -r raymeta/
	-rm -r fermi/
	-rm -r masurca/
	-rm *.contigs.fa *.ctgs_filt.fa
