#!/usr/bin/env python
"""
Script to generate sbatch files to run the pipeline in the Uppmax Cluster

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
import os
import os.path

#Time of script execution and logging module
import time
import logging

# import re
#
# import itertools
# from collections import *

#****************Begin of Main ***************
def main(args):
	steps = ["qc","qf","rmcont","asm","tax_diamond","metaphlan2" ]
	modules = {
		"qc": ["+java/sun_jdk1.7.0_25" ],
		"qf" :["+bioinfo-tools","+FLASH/1.2.11","+java/sun_jdk1.7.0_25" ],
		"rmcont":["+bioinfo-tools","+bowtie2/2.2.6","+picard/2.0.1"],
		"asm":["-gcc","+bioinfo-tools","+bwa/0.7.13","+picard/2.0.1","+spades/3.7.0","+fermi/1.1-r751-beta"],
		"tax_blast":["+bioinfo-tools","+blast/2.2.31+"],
		"tax_diamond":["+java/sun_jdk1.8.0_40","+bioinfo-tools","+diamond/0.7.9"],
		"metaphlan2" : ["+bioinfo-tools","+bowtie2/2.2.6"],
	}

	#Folders to add to PATH
	ext_tool_path = {
		"qc": 		["/proj/b2011088/tools/anaconda/bin" ],
		"qf":		["/proj/b2011088/tools/anaconda/bin" ],
		"rmcont":	["/proj/b2011088/tools/anaconda/bin" ],
		"asm":		["/proj/b2011088/tools/anaconda/bin",
				 	"/proj/b2011088/tools/megahit/1.0.2"],
		"tax_blast":	["/proj/b2011088/tools/anaconda/bin" ],
		"tax_diamond":	["/proj/b2011088/tools/anaconda/bin" ],
		"metaphlan2":	["/proj/b2011088/tools/anaconda/bin" ]
	}

	step_partition = {
		"qc":		"core -n 2",
		"qf":		"core -n 3",
		"rmcont":	"core -n 8",
		"asm": 	 	"core -n 8",
		"tax_blast":	"core -n 8",
		"tax_diamond":	"node -n 16",
		"metaphlan2":	"core -n 8"
	}

	step_duration = {
		"qc": 		"06:00:00",
		"qf" : 		"06:00:00",
		"rmcont": 	"06:00:00",
		"asm": 		"06:00:00",
		"tax_blast":	"4-00:00:00",
		"tax_diamond":	"4-00:00:00",
		"metaphlan2":	"3:00:00"
	}

	step_makefile_rule = {
		"qc":		"1_raw_qc",
		"qf":		"2_qf_qc",
		"rmcont": 	"3_hostfiltering",
		"asm": 		"4_assembly",
		"tax_blast": "5_tax_blast",
		"tax_diamond": "5_tax_diamond",
		"metaphlan2":	"5_metaphlan"
	}

	#Create log file for sbatch logs
	log_folder = "log/"
	if not os.path.exists(log_folder):
		os.makedirs(log_folder)

	for step in steps:
		step_filename = os.path.join(args.output_folder,"run_"+step+".sh")
		with open(step_filename, "w") as fh:
			fh.write("#!/bin/bash\n")
			fh.write("#SBATCH -A {}\n".format(args.project_id))
			fh.write("#SBATCH -p {}\n".format(step_partition[step]))
			fh.write("#SBATCH -t {}\n".format(step_duration[step]))
			fh.write("#SBATCH -J {}_{}\n".format(args.sample_name,step))
			fh.write("#SBATCH -o log/{}_{}-%j.out\n".format(args.sample_name,step))
			fh.write("#SBATCH -e log/{}_{}-%j.err\n".format(args.sample_name,step))

			fh.write("#SBATCH --mail-user {}\n".format(args.email) )
			fh.write("#SBATCH --mail-type=ALL\n\n")

			fh.write( modules_to_load(modules[step]) )
			#Unload default python and samtools (using anaconda)
			fh.write( modules_to_load(["-python","-samtools"]))

			fh.write("set -euo pipefail\n\n") #Fast failing for bash script

			#Add external tools to PATH
			fh.write( append_to_path(ext_tool_path[step]) )

			#Assumes cfg/uppmax.cfg is the configuration to use
			#TODO: choose base configuration file from arguments
			fh.write("make -r sample_name={} cfg_file=cfg/uppmax.cfg {}\n".format(args.sample_name,step_makefile_rule[step]))

#*****************End of Main**********************

def modules_to_load(modules):
	module_string = ""
	for mod in modules:
		if mod[0] == "+":
			module_string += "module load "+mod[1:]+"\n"
		elif mod[0] == "-":
			module_string += "module unload "+mod[1:]+"\n"
		else:
			logging.error("Module definition malformed")
			sys.exit(1)

	return module_string

def append_to_path(tool_paths):
	new_path_str = "export PATH=\"" + ":".join(tool_paths) + ":$PATH\"\n"
	return new_path_str

def validate_args(args):
	return True


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Script that generates sbatch scripts to run pipeline in a SLURM-based cluster like Uppmax")

	parser.add_argument("sample_name",help="Sample name that will be used on the pipeline")
	parser.add_argument("-o","--output-folder", default="./", help="Folder where to output the scripts" )
	parser.add_argument("-A","--project_id", default="b2011088",help="Slurm project id")
	parser.add_argument("--email", default="mauricio.barrientos@ki.se",help="Email to send run reports to")
	args = parser.parse_args()

	#Initialize log
	log_level = logging.INFO
	logging.basicConfig(stream=sys.stderr,level=log_level)

	if validate_args(args):
		time_start = time.time()
		main( args )
		logging.info("Time elapsed: "+str(time.time() - time_start)+"\n")
	else:
		logging.error("Invalid arguments. Exiting script\n")
		sys.exit(1)
