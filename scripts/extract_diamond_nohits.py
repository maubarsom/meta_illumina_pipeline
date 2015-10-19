#!/usr/bin/env python
# Script to extract a fasta with all

# Author: Mauricio Barrientos-Somarribas
# Email:  mauricio.barrientos@ki.se

# Copyright 2015 Mauricio Barrientos-Somarribas
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

import sys
import argparse
import os.path
import time
import gzip
import re

def main(args):
	#Open diamond file
	if args.diamond_sam.endswith(".gz"):
		diamond_fh = gzip.open(args.diamond_sam,"r")
	else:
		diamond_fh = open(args.diamond_sam,"r")

	seqs_with_hits = []
	reading_sam_header = True
	for line in diamond_fh:
		if reading_sam_header and re.match(r"@[A-Za-z]{2}\t",line):
			continue
		else:
			is_header = False
			seqs_with_hits.append( line.split("\t")[0] )

	diamond_fh.close()
	diamond_hits = set(seqs_with_hits)

	#Subsample the fasta with reads that did not appear in the sam
	seqs_in_fasta = 0
	seqs_with_no_hits = 0
	with open(args.fasta_input,"r") as input_f:
		for line in input_f:
			if is_fasta_header(line):
				seqs_in_fasta += 1
				seq_id = extractReadName(line)
				has_diamond_hit = seq_id in diamond_hits
				if not has_diamond_hit:
					seqs_with_no_hits += 1
					args.output_file.write(line)
			else:
				if not has_diamond_hit:
					args.output_file.write(line)

		args.output_file.close()

		sys.stderr.write( "****Stats****\n")
		sys.stderr.write( "# of seqs in fasta:\t{}\n".format(seqs_in_fasta) )
		sys.stderr.write( "# of unique seqs in diamond sam:\t{}\n".format( len(diamond_hits)) )
		sys.stderr.write( "# of seqs written to output:\t{}\n".format(seqs_with_no_hits) )


#*****************End of main**********************
def is_fasta_header(line):
	return line.startswith(">")

def extractReadName(line):
	return line[1:].rstrip("\n")

def validate_args(args):
	if args.fasta_input and not os.path.isfile(args.fasta_input):
		sys.stderr.write("Error! "+args.fasta_input+" does not exist!\n")
		sys.exit(1)

	if not args.output_file:
		args.output_file = sys.stdout
	else:
		args.output_file = open( args.output_file, "w")

	return True

if __name__ == '__main__':
	start_time = time.time()
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Takes a diamond sam output file and the original \
		fasta and lists the sequences that have no hits")
	parser.add_argument("fasta_input",help="Fasta file to process")
	parser.add_argument("diamond_sam",help="Diamond sam output")
	parser.add_argument("-o","--output-file",help="Name for the output file. Default: stdout")

	args = parser.parse_args()

	if validate_args(args):
		main( args )
		sys.stderr.write("Time elapsed: {} seconds\n".format(time.time() - start_time))

	else:
		sys.stderr.write("Invalid arguments. Exiting script")
		sys.exit(1)
