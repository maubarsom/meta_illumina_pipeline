#Requires:
# wget
# git
# GNU sed as sed
# Python 2.7 and 3, with pip installed for both
# tar , unzip
# make and cmake(for bamtools)

SHELL := /bin/bash

ROOT_FOLDER := $(shell pwd)

#Python configuration.
# Commands like pip or python might require sudo or a virtualenv
pip_python2 := pip
pip_python3 := source activate py3k && pip
python2_bin := python
python3_bin := source activate py3k && python3

.PHONY: maars viral-discovery
.PHONY: fastqc jellyfish2
.PHONY: cutadapt nesoni
.PHONY: raymeta masurca

#******************************************************
#		MAARS pipeline tools
#*****************************************************
maars: seqtk samtools prinseq
maars: fastqc
maars: cutadapt nesoni
maars: bwa picard-tools
maars: megahit
maars: diamond kraken

#******************************************************
#		Viral discovery tools
#*****************************************************
viral-discovery: seqtk samtools prinseq sga
viral-discovery: fastqc
viral-discovery: cutadapt nesoni
viral-discovery: bwa picard-tools
viral-discovery: spades iva fermi raymeta masurca
viral-discovery: ncbi-blast hmmer

#Create bin folder
$(shell mkdir -p bin/)

#**************************************************************************************
#******************        General purpose tools      *********************************
#**************************************************************************************
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

picard-tools/:
	wget -N https://github.com/broadinstitute/picard/releases/download/1.124/picard-tools-1.124.zip
	unzip picard-tools-*.zip
	mv picard-tools-*/ picard-tools/

#**************************************************************************************
#******************        QC and QF       ********************************************
#**************************************************************************************

fastqc:
	wget -N http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
	unzip fastqc_*.zip
	rm -rf fastqc
	mv FastQC fastqc
	chmod ug+x fastqc/fastqc
	@echo "This installation assumes you have at least 16 gigs of ram"
	sed -i "s/Xmx250m/Xmx16g/" fastqc/fastqc
	cd bin && ln -sf ../fastqc/fastqc

prinseq/:
	wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download -O prinseq.tar.gz
	tar -xzf prinseq.tar.gz
	mv prinseq-* prinseq/
	cp prinseq/prinseq-lite.pl bin/
	chmod ug+x bin/prinseq-lite.pl

jellyfish2:
	echo "No rule available yet" && exit 1

cutadapt:
	$(pip_python2) install cutadapt

nesoni:
	$(pip_python2) install nesoni

fqtrim:
	wget -N http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.94.tar.gz
	tar -xzf fqtrim-*.tar.gz
	mv fqtrim-*/ fqtrim
	cd fqtrim && make release
	cp fqtrim/fqtrim bin/

pandaseq:
	echo "Requires zlib, bzip2 and libtool"
	#sudo apt-get install zlib1g-dev libbz2-dev libltdl-dev libtool
	git clone http://github.com/neufeld/pandaseq.git/
	cd pandaseq
	./autogen.sh && ./configure --prefix=`pwd`/dist && make && make install && make clean
	cp pandaseq/dist/bin/* bin/*

#**************************************************************************************
#******************        READ MAPPERS/ALIGNERS     **********************************
#**************************************************************************************
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

#**************************************************************************************
#******************        ASSEMBLERS     *********************************************
#**************************************************************************************

fermi:
	wget -N https://github.com/downloads/lh3/fermi/fermi-1.1.tar.bz2
	tar -xjf fermi-*
	mv fermi-*/ fermi
	cd fermi && make
	cp fermi/fermi bin/
	cp fermi/run-fermi.pl bin/

megahit:
	git clone https://github.com/voutcn/megahit.git
	cd megahit && make
	cp megahit/megahit bin/

SPAdes:
	wget -N http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
	tar -xzf SPAdes*.tar.gz && mv SPAdes-*/ SPAdes
	cp SPAdes/bin/* bin/

#Assumes python3 installed from anaconda in a virtualenv called py3k
Fastaq/:
	git clone https://github.com/sanger-pathogens/Fastaq.git
	cd Fastaq && $(python3_bin) setup.py install

kmc:
	mkdir -p kmc/
	cd kmc && wget -N http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc
	cd kmc && wget -N http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc_dump
	cp kmc/* bin/

#Assumes python 3 install from anaconda in a virtual env called py3k
iva: Fastaq MUMmer smalt samtools
	$(pip_python3) install networkx
	$(pip_python3) install pysam
	wget -N https://github.com/sanger-pathogens/iva/archive/v0.11.0.tar.gz
	tar -xzf v0.11.0.tar.gz
	mv iva*/ iva
	cd iva && $(python3_bin) setup.py install

MUMmer:
	wget -N http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
	tar -xzf MUMmer*.tar.gz
	mv MUMmer*/ MUMmer
	cd MUMmer && make install
	for x in `find MUMmer/ -maxdepth 1 -executable -not -type d -exec basename {} \; `; do ln `pwd`/MUMmer/$x bin/$x ; done

smalt:
	wget -N http://downloads.sourceforge.net/project/smalt/smalt-0.7.6-static.tar.gz
	tar -xzf smalt*.tar.gz
	mv smalt*/ smalt/
	cd smalt && mkdir -p dist && ./configure --prefix=`pwd`/dist && make install
	cp smalt/dist/bin/* bin/*

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

#************************************************************************
#******************    GENE PREDICTION       ****************************
#************************************************************************
FragGeneScan:
	wget http://sourceforge.net/projects/fraggenescan/files/latest/download -O fgs.tar.gz
	tar -xzf fgs.tar.gz
	mv FragGeneScan*/ FragGeneScan
	cd FragGeneScan && make clean && make fgs

#************************************************************************
#******************        DB SEARCHES       ****************************
#************************************************************************
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

diamond:
	wget http://www-ab.informatik.uni-tuebingen.de/data/software/diamond/download/public/diamond-intel64-linux.tar.gz
	tar -xzf diamond*.tar.gz
	mv diamond bin/

kraken:
	wget -N http://ccb.jhu.edu/software/kraken/dl/kraken-0.10.4-beta.tgz
	tar -xzf kraken*.tgz
	mv kraken*/ kraken
	mkdir -p kraken/dist
	cd kraken && bash install_kraken.sh dist/
	cp kraken/dist/* bin/
