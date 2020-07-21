#!/usr/bin/env python3
"""
# Name: Ecsv2merge
# Author: Marie Lefebvre
# Aims: Merge blastx results (.csv) of contigs and singletons
#
"""
import csv, os
import argparse
import logging as log


def main ():
	options = _set_options()
	_set_log_level(options.verbosity)
	_print_csv(options)


def _print_csv(options):
	f = open(options.output, "w+")
	if os.path.exists(options.c_file):
		# contigs file
		with open(options.c_file) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter="\t", quotechar='"')
			for row in csv_reader:
				f.write('\t'.join(row))
				f.write("\n")
	if os.path.exists(options.s_file):
		# singletons file
		with open(options.s_file) as csv_file:
			reader = csv.reader(csv_file, delimiter="\t", quotechar='"')
			headers = True
			for row in reader:
				# skip headers
				if headers is False:
					f.write('\t'.join(row))
					f.write("\n")
				headers = False
	f.close()


def _set_options():
	parser = argparse.ArgumentParser()
	parser.add_argument('-c','--contigs',help='Csv files with blastx results of contigs.', action='store', required=True,dest='c_file')
	parser.add_argument('-s','--singletons',help='Csv files with blastx results of singletons.', action='store',required=True,dest='s_file')
	parser.add_argument('-o','--out',help='The output file.',action='store', type=str,default='./', dest='output')
	parser.add_argument('-v','--verbosity',help='Verbose level', action='store',type=int,choices=[1,2,3,4],default=1)
	args = parser.parse_args()
	return args

def _set_log_level(verbosity):
	if verbosity == 1:
		log_format = '%(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.INFO,format=log_format)
	elif verbosity == 3:
		log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.DEBUG,format=log_format)

if __name__ == "__main__":
	main()