#Requires:
# wget
# git
# GNU sed as sed
# Python 2.7 and pip
# tar , unzip
# make and cmake(for bamtools)

SHELL := /bin/bash

ROOT_FOLDER := $(shell pwd)

.PHONY: fastqc jellyfish2

.PHONY: cutadapt nesoni
.PHONY: abyss raymeta masurca

all: seqtk samtools prinseq sga
all: fastqc jellyfish2
all: cutadapt nesoni fqtrim
all: bwa stampy picard-tools
all: fermi abyss raymeta masurca
all: ncbi-blast hmmer

#Create bin folder
$(shell mkdir -p bin/)

seqtk/:
	git clone https://github.com/lh3/seqtk.git
	cd seqtk && make
	cp seqtk/seqtk bin/

samtools/:
	wget http://sourceforge.net/projects/samtools/files/latest/download?source=files -O samtools.tar.bz2
	tar -xjf samtools.tar.bz2
	mv samtools-*/ samtools/
	cd samtools/ && make
	cp samtools/samtools bin/

prinseq/:
	wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download -O prinseq.tar.gz
	tar -xzf prinseq.tar.gz
	mv prinseq-* prinseq/
	cp prinseq/prinseq-lite.pl bin/
	chmod ug+x bin/prinseq-lite.pl

#Quality Control
fastqc:
	wget -N http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
	unzip fastqc_*.zip
	rm -rf fastqc
	mv FastQC fastqc
	chmod ug+x fastqc/fastqc
	@echo "This installation assumes you have at least 16 gigs of ram"
	sed -i "s/Xmx250m/Xmx16g/" fastqc/fastqc
	cd bin && ln -sf ../fastqc/fastqc

jellyfish2:
	echo "No rule available yet"

#SGA and sga-dependencies
sga: sparsehash bamtools jemalloc
	git clone https://github.com/jts/sga.git
	cd sga/src && ./autogen.sh
	cd sga/src && ./configure --with-sparsehash=$(ROOT_FOLDER)/sparsehash --with-bamtools=$(ROOT_FOLDER)/bamtools --with-jemalloc=$(ROOT_FOLDER)/jemalloc/lib --prefix=$(ROOT_FOLDER)
	cd sga/src && make && make install
	#Copy preqc report script to bin
	cp sga/src/bin/sga-preqc-report.py bin/

sparsehash:
	wget -N https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz
	tar -xzf sparsehash*.tar.gz
	cd sparsehash-*/ && ./configure --prefix=$(ROOT_FOLDER)/sparsehash && make && make install

bamtools:
	git clone https://github.com/pezmaster31/bamtools.git
	mkdir -p bamtools/build
	cd bamtools/build && cmake ..
	cd bamtools/build && make

jemalloc:
	wget -N http://www.canonware.com/download/jemalloc/jemalloc-3.6.0.tar.bz2
	tar -xjf jemalloc-*.tar.bz2
	cd jemalloc-*/ && ./configure --prefix=$(ROOT_FOLDER)/jemalloc && make && make install

#QF
cutadapt:
	pip install cutadapt

nesoni:
	pip install nesoni

fqtrim:
	wget -N http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.94.tar.gz
	tar -xzf fqtrim-*.tar.gz
	mv fqtrim-*/ fqtrim
	cd fqtrim && make release
	cp fqtrim/fqtrim bin/

#Cont rm
bwa/:
	wget -N http://sourceforge.net/projects/bio-bwa/files/latest/download?source=files -O bwa.tar.bz2
	tar -xjf bwa.tar.bz2
	cp bwa.kit/bwa bin/

#Requires python2.x and python2.x-config, where x = 6 or 7
#This is breaking, maybe requires python-dev package?
stampy/:
	wget -N http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
	tar -xzf Stampy-latest.tgz
	mv stampy-*/ stampy
	cd stampy && make
	cd bin && ln -s ../stampy/stampy.py

picard-tools/:
	wget -N https://github.com/broadinstitute/picard/releases/download/1.124/picard-tools-1.124.zip
	unzip picard-tools-*.zip
	mv picard-tools-*/ picard-tools/

#Assemblers
fermi:
	wget -N https://github.com/downloads/lh3/fermi/fermi-1.1.tar.bz2
	tar -xjf fermi-*
	mv fermi-*/ fermi
	cd fermi && make
	cp fermi/fermi bin/
	cp fermi/run-fermi.pl bin/

#Taxonomic searches
ncbi-blast:
	wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.??+-x64-linux.tar.gz
	tar -xzf ncbi-blast*.tar.gz
	mv ncbi-blast-*/ ncbi-blast
	cp ncbi-blast/bin/* bin

hmmer:
	wget -N ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-*-linux-intel-x86_64.tar.gz
	tar -xzf hmmer-*-linux-intel-*.tar.gz
	mv hmmer-*/ hmmer
	cp hmmer/binaries/?hmm* bin/
	cp hmmer/binaries/hmm* bin/

FragGeneScan:
	wget http://sourceforge.net/projects/fraggenescan/files/latest/download -O fgs.tar.gz
	tar -xzf fgs.tar.gz
	mv FragGeneScan*/ FragGeneScan
	cd FragGeneScan && make clean && make fgs

diamond:
	wget http://www-ab.informatik.uni-tuebingen.de/data/software/diamond/download/public/diamond-intel64-linux.tar.gz
	tar -xzf diamond*.tar.gz
	mv diamond bin/
