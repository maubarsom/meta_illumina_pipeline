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

.PHONY: all kraken_reports

all: kraken_reports

kraken_reports: $(call ctg_outfile,kraken,kraken.report)
#kraken_reports: $(call read_outfiles,kraken,kraken.report,pe se)all

#*************************************************************************
#Call to Kraken - Salzberg
#*************************************************************************
#Other flags: --fastq-input
kraken/%_kraken.out: $(ctg_folder)/%.fa
	mkdir -p kraken
	@echo -e "\nClassifying $* with Kraken\n\n"
	kraken --preload --db $(kraken_db) --threads $(threads) $^ > $@

%_kraken.report: %_kraken.out
	@echo -e "\nCreating Kraken report for $* \n\n"
	kraken-report --db $(kraken_db) $^ > $@
