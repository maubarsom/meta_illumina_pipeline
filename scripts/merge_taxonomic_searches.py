#!/usr/bin/env python
"""
merge_taxonomic_searches.py is a script that takes the different sequence database searches applied
to a same set of contigs

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

import itertools
from collections import *

from Bio.Blast import NCBIXML
from pymongo import MongoClient

#****************Begin of Main ***************
def main(args):

	#Open Mongodb connection
	mongo_client = MongoClient()
	mongo_db = mongo_client["db_name"]
	mongo_collection = mongo_db["collection_name"]

	assemblers = ["sga","raymeta","masurca","fermi"]
	databases = ["nt","Refseq Viral nucl","nr","Refseq Viral prot","Pfam-A","vFam-A"]


	# for asm in assemblers:
	# 	for db in databases:
	# 		pass

#upsert = insert + update
def upsert_into_mongo(mongo_collection,json_dict):
	mongo_collection.update({"seq_id":json_dict["seq_id"]} ,{"$set": json_dict}, upsert=True )

def hmmr_parser():
	pass

def kraken_parser():

def blast_parser(blast_xml_file,blast_type,asm,db,mongo_collection, evalue_threshold=0.001,max_hits=3, ):
	with open(blast_xml_file,"r") as file_handle:
		iterations = NCBIXML.parse(file_handle)
		logging.info("Parsing XML file...")
		for iteration_num,it in enumerate(iterations):
			hit_num = 0
			#Create dictionary json object for Mongo
			it_json = {"seq_id":it.query,"seq_len":it.query_length,"asm":asm,"classification":None}
			it_json[db] = []
			for hit_num, hit in enumerate(it.alignments):
				hsps = []
				for hsp in hit.hsps:
					if hsp.expect < evalue_threshold:
						hsps.append( { "e-value": hsp.expect, "alignment_length":hsp.align_length, "identities": hsp.identities,"positives": hsp.positives,"gaps": hsp.gaps,"q_start": hsp.query_start, "q_end":hsp.query_end } )
				#Add least one significant hsp
				if len(hsps) > 0:
					it_json[db].append({"hit_id":hit.hit_id,"hit_def":hit.hit_def,"hsps":hsps})

				if hit_num == max_hits:
					break
			#Insert into mongo
			upsert_into_mongo(mongo_collection,it_json)
			#Keep track of number of iterations processed

			if iteration_num % 2000 == 0 :
				logging.info("\tIterations processed: "+str(iteration_num)+"\r")

	logging.info("XML parsed! " + str(iteration_num) + " iterations\n")

#*****************End of Main**********************
def validate_args(args):
	return True

if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="The program agregates hit results from different programs into a mongodb ")

	parser.add_argument("annotation_files",help="Blast xml or hmmer output files to parse",nargs="+")
	parser.add_argument("-o","--output-file", type=argparse.FileType('w'), default=sys.stdout, help="Name of the output file" )
	#Logging
	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")
	#Other args
	parser.add_argument("--max_hits",type=int,default="10",help="Max number of hits to report for each query sequence")

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
