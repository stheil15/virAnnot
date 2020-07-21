#!/usr/bin/env python3
"""
# Name: Filter
# Author: Eden Darnige
# Aims: Filter assembly fasta output by sequence length
#
"""
import argparse
import logging as log
from Bio import SeqIO


def main ():
	options = _set_options()
	_set_log_level(options.verbosity)
	_filter_contigs(options.contigs_in, options.l, options.output)


def _filter_contigs(contigs_in,l,output):
        contigs_file = contigs_in
        filtered_contigs_file = output
        length = l

        with open(contigs_file,'r') as in_handle:
            with open(filtered_contigs_file, 'w') as out_handle:
                for entry in SeqIO.parse(in_handle, "fasta"):
                    if len(entry.seq) >= length:
                        out_handle.write(">%s\n%s\n" % (entry.id, entry.seq))


def _set_options():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input',help='Fasta contig files from assembly output.', type=str, action='store', required=True,dest='contigs_in')
	parser.add_argument('-l','--length',help='Threshold length of contigs to keep.', type=int, action='store',required=True, dest='l')
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