#!/usr/bin/env python
"""
Script that uses output from cutadapt to quickly detect fully overlapping pairs.
It is based on the fact that if sequencing adapters are trimmed from both paired-ends, the resulting
fragment needs to be shorter than the pair-end length

Depends of cutadapt seqio and xopen modules from version 1.6

Author: Mauricio Barrientos-Somarribas
Email:  mauricio.barrientos@ki.se

Copyright 2014 Mauricio Barrientos-Somarribas

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import sys
import argparse
import os.path

#Time of script execution and logging module
import time
import logging

import re
import math

import itertools
from collections import *

from cutadapt import seqio,xopen
from distance import hamming
import ctypes

#Data Analysis libs
import numpy as np

#****************Begin of Main ***************
def main(args):
	#Global vars
	global raw_read_len
	global min_trim

	raw_read_len = args.raw_read_length
	min_trim = args.min_trim

	paired_reader = seqio.PairedSequenceReader(args.R1, args.R2 , fileformat="fastq")
	stats = FragmentStats(raw_read_len)

	try:
		out_r1 = xopen.xopen(args.output_prefix+"_R1.fq.gz", "w")
		out_r2 = xopen.xopen(args.output_prefix+"_R2.fq.gz", "w")
		small_fragments = xopen.xopen(args.output_prefix+"_single.fq.gz","w")
	except Exception:
		logging.error("An error has occured creating one of the output files")
		sys.exit(1)

	for read1, read2 in paired_reader:
		is_short_fragment = False
		#If both paired-ends were trimmed "confidently"
		r1_len, r2_len = len(read1.sequence), len(read2.sequence)
		if max(r1_len,r2_len) < (raw_read_len - min_trim):
			aligned_r1, aligned_r2 = align_sequences(read1,read2)
			is_short_fragment = is_fragment(aligned_r1,aligned_r2)

		if is_short_fragment:
			stats.add_small_fragment( len(aligned_r1.sequence) )
			consensus_fragment = get_consensus(aligned_r1,aligned_r2)
			consensus_fragment.write(small_fragments)
		else:
			stats.add_non_fragment()
			read1.write(out_r1)
			read2.write(out_r2)

	out_r1.close()
	out_r2.close()
	small_fragments.close()
	logging.info(str(stats)+"\n")

#*****************End of Main**********************
class FragmentStats:
	def __init__(self,read_len):
		self.total_reads = 0
		self.fragment_count = 0
		self.fragment_size_dist = np.zeros(read_len,dtype=np.uint32)

	def add_small_fragment(self,fragment_len):
		self.total_reads += 1
		self.fragment_count += 1
		self.fragment_size_dist[fragment_len] += 1

	def add_non_fragment(self):
		self.total_reads += 1

	def __str__(self):
		out_str = "Extract fragment stats\n"
		out_str += "Reads processed: "+str(self.total_reads)+"\n"
		out_str += "Small fragments detected: {0}({1:.2%})\n".format(self.fragment_count,
																	float(self.fragment_count)/self.total_reads)
		out_str += "Small fragment size distribution\n"
		out_str += "Fragment_len\tcount\n"
		out_str += "\n".join([ "{0}\t{1}".format(frag_size,count)
								for frag_size, count in enumerate(self.fragment_size_dist) ] )
		return out_str


def rev_complement(seq):
	return seq[::-1].upper().replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper()

def is_fragment(aligned_r1,aligned_r2):
	r1_len, r2_len = len(aligned_r1.sequence), len(aligned_r2.sequence)
	is_frag = r1_len == r2_len and hamming(aligned_r1.sequence , aligned_r2.sequence ) < (r1_len * 0.05)
	#is_frag = is_frag and abs( r1_len - r2_len ) <= max_fragment_len_dif:
	return is_frag

#Assumes sequences are the same length, with a hamming distance of less than 0.05% of the length
def align_sequences(r1,r2):
	new_r2 = seqio.Sequence(r2.name, rev_complement(r2.sequence), r2.qualities[::-1])
	return r1,new_r2

#Assumes sequences are aligned and the same length
def get_consensus(aligned_r1,aligned_r2):
	consensus_seq = ctypes.create_string_buffer(aligned_r1.sequence)
	consensus_qual = ctypes.create_string_buffer(aligned_r1.qualities)

	for pos in range(len(aligned_r1)):
		if aligned_r2.qualities[pos] > consensus_qual[pos]:
			consensus_seq[pos] = aligned_r2.sequence[pos]
			consensus_qual[pos] = aligned_r2.qualities[pos]

	#Removes illumina's pair info from the consensus sequence name
	consensus_obj = seqio.Sequence(aligned_r1.name.split(" ")[0],
									consensus_seq.value,
									consensus_qual.value)

	return consensus_obj

def validate_args(args):
	valid_args= True
	if not os.path.isfile( args.R1 ):
		logging.error( args.R1+" is not a valid file")
		valid_args = False
	if not os.path.isfile( args.R2 ):
		logging.error( args.R2+" is not a valid file")
		valid_args = False
	if valid_args and not args.output_prefix:
		stem_re = re.match(r"(.+)_R?[12]\.f(ast)?q(\.gz)?$",args.R1 )
		args.output_prefix = stem_re.group(1)
	return valid_args


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Script to process fastq files after adapter removal and extracts sequences fragments smaller than read length")
	parser.add_argument("R1",help="Fastq with forward paired-end")
	parser.add_argument("R2",help="Fastq with reverse paired-end")

	parser.add_argument("-o","--output-prefix", default=None, help="Prefix of the output files" )
	parser.add_argument("--raw_read_length", default=301,type=int, help="Length of raw reads (before adapter trimming). Default: 301" )
	parser.add_argument("--min-trim", default=10,type=int, help="Minimum number of bases trimmed to consider the adapter removed was not spurious. Default: 10" )

	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")

	args = parser.parse_args()

	if validate_args(args):
		#Initialize log
		log_level = logging.INFO
		if args.log_file:
			logging.basicConfig(filename=args.log_file,level=log_level)
		else:
			logging.basicConfig(stream=sys.stderr,level=log_level)

		time_start = time.time()
		main( args )
		logging.info("Time elapsed: "+str(time.time() - time_start)+"\n")
	else:
		logging.error("Invalid arguments. Exiting script\n")
		sys.exit(1)
