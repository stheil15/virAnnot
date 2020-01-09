#!/usr/bin/env python3
"""
# Name: Rps2merge
# Author: Marie Lefebvre
# Objective: merge OTU results with blastx and rpsblast results
#
"""
import csv, re, os
import argparse
import logging as log
from collections import defaultdict
from glob import glob


def main ():
	options = _set_options()
	_set_log_level(options.verbosity)
	cdd_dict = _read_pfam(options)
	blast_dict = _read_blast(options)
	rps_dict = _read_rps(options)
	merge_dict = _merge_data(cdd_dict, blast_dict, rps_dict)
	_print_csv(merge_dict, options.output)


# Read otu files and retrieve information:
# OTU number, taxonomy, contigs list
# Render dictionnary with cdd as keys
# dict : {'int': [cdd, otu, taxo, [contigs]]}
def _read_pfam(options):
	pfam_folders = glob(options.rps_folder + "/*/")
	otu_file = "cluster_nb_reads.csv"
	pfam_otu_dict = defaultdict(list)
	i = 0
	# for each cdd folder
	for p_f in pfam_folders:
		otu_file_path = p_f + otu_file
		cdd = os.path.basename(os.path.normpath(p_f))
		if os.path.exists(otu_file_path):
			with open(otu_file_path, 'rt') as csvfile:
				spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
				# next(spamreader, None)  # skip the headers
				headers = True
				for row in spamreader:
					row_length = len(row)
					# skip headers
					if headers == False:
						# append datat to dict
						pfam_otu_dict[i].append(cdd) # cdd
						pfam_otu_dict[i].append(row[0]) # otu
						pfam_otu_dict[i].append(str(row[row_length-2])) # taxonomy
						pfam_otu_dict[i].append(row[row_length-1]) # contigs
						i = i + 1
					headers = False
	return(pfam_otu_dict)



# Read input blastx file
# row[1] -> contig name
# row[3] -> query_length
# row[4] -> accession
# row[5] -> description
# row[6] -> organism (can be unknown)
# row[7] -> identity
# row[11] -> evalue
# row[14] -> taxonomy
# Return dictionnary
def _read_blast(options):
	blast_dict = defaultdict(list) # dict result
	# for each blastx files
	for i in range(0, len(options.bx_files)):
		blast_file = options.bx_files[i]
		reader = csv.reader(blast_file, delimiter="\t", quotechar='"')
		headers = True
		for row in reader:
			# skip headers
			if headers == False:
				contig_name = row[1]
				try:
					blast_dict[contig_name].append(row[6]) # organism
				except IndexError:
					blast_dict[contig_name].append("unknown") # organism
				try:
					blast_dict[contig_name].append(row[3]) # query length
				except IndexError:
					blast_dict[contig_name].append("")
				try:
					blast_dict[contig_name].append(row[4]) # accession
				except IndexError:
					blast_dict[contig_name].append("")
				try:
					blast_dict[contig_name].append(row[5]) # description
				except IndexError:
					blast_dict[contig_name].append("")
				try:
					blast_dict[contig_name].append(row[7]) # identity
				except IndexError:
					blast_dict[contig_name].append("")
				try:
					blast_dict[contig_name].append(row[11]) # evalue
					blast_dict[contig_name].append(row[14]) # taxonomy
				except IndexError:
					blast_dict[contig_name].append("") # evalue
					blast_dict[contig_name].append("") # taxonomy
			headers = False
	return(blast_dict)

# Read rps files results
# row[0] -> contig name
# row[4] -> evalue
# row[9] -> taxonomy
# Return dictionnary
def _read_rps(options):
	rps_dict = defaultdict(list)
	# For each rps blast files
	for i in range(0, len(options.pfam_files)):
		rps_file = options.pfam_files[i]
		reader = csv.reader(rps_file, delimiter="\t")
		headers = True
		# for each line
		for row in reader:
			# Skip headers
			if headers is False:
				contig_name = row[0]
				if row[1] != "no_hit":
					rps_dict[contig_name].append(row[4])
					rps_dict[contig_name].append(row[9])
			headers = False
	return(rps_dict)


# Merge pfam, blastx and rps blast informations
# cdd, otu, taxo, contigs, blastx_organism, blastx_identity, blastx_evalue, blastx_taxo, 
# rps_evalue, rps_taxo
def _merge_data(cdd_dict, blast_dict, rps_dict):
	merge_dict = defaultdict(list)
	i = 0
	for key in cdd_dict:
		# Get contig list
		contig_list = cdd_dict[key][3].split(",")
		# # Randomly select a contig in the list
		# random_contig = contig_list[randint(0, len(contig_list)-1)]
		# Select le longest contig in the list
		random_contig = max(contig_list, key=len)
		merge_dict[i].append(cdd_dict[key][0]) # cdd
		merge_dict[i].append(cdd_dict[key][1]) # otu
		merge_dict[i].append(cdd_dict[key][2]) # taxo
		merge_dict[i].append(cdd_dict[key][3]) # contigs
		merge_dict[i].append(random_contig) # selected contig
		# FROM BLASTX
		# Check if contig exists
		if len(blast_dict[random_contig]) != 0:
			merge_dict[i].append(blast_dict[random_contig][0]) # organism
			merge_dict[i].append(blast_dict[random_contig][1]) # query length
			merge_dict[i].append(blast_dict[random_contig][2]) # accession
			merge_dict[i].append(blast_dict[random_contig][3]) # description
			merge_dict[i].append(blast_dict[random_contig][4]) # identity
			merge_dict[i].append(blast_dict[random_contig][5]) # evalue
			merge_dict[i].append(blast_dict[random_contig][6]) # taxonomy
		else:
			merge_dict[i].append("unknown") # organism
			merge_dict[i].append("unknown") # query length
			merge_dict[i].append("unknown") # accession
			merge_dict[i].append("unknown") # description
			merge_dict[i].append("unkonwn") # identity
			merge_dict[i].append("unkonwn") # evalue
			merge_dict[i].append("unkonwn") # taxo
		# FROM RPS
		if len(rps_dict[random_contig]) != 0:
			merge_dict[i].append(rps_dict[random_contig][0]) # evalue
			merge_dict[i].append(rps_dict[random_contig][1]) # taxo
		else:
			merge_dict[i].append("unknown") # evalue
			merge_dict[i].append("unknown") # taxo
		i = i + 1
	return(merge_dict)


# Print merged datat to csv file
def _print_csv(merge_dict,out):
	cwd = os.getcwd()
	f = open(cwd + '/' + out, "w+")
	# set headers
	f.write("cdd\totu\ttaxo\tcontigs\tselected_contig\tblastx_organism\tblastx_query_length\tblastx_accession\tblastx_description\tblastx_identity\tblastx_evalue\tblastx_taxo\trps_evalue\trps_taxo\n")
	for k in merge_dict:
		print(merge_dict[k])
		f.write(merge_dict[k][0] + "\t" + merge_dict[k][1] + "\t") # cdd + otu
		f.write(merge_dict[k][2] + "\t" + merge_dict[k][3] + "\t") # taxo + contigs
		f.write(merge_dict[k][4] + "\t" + merge_dict[k][5] + "\t") # selected contig + organism
		f.write(merge_dict[k][6] + "\t" + merge_dict[k][7] + "\t") # query length + accession
		f.write(merge_dict[k][8] + "\t" + merge_dict[k][9] + "\t") # description + identity
		f.write(merge_dict[k][10] + "\t" + merge_dict[k][11] + "\t") # evalue + taxo
		f.write(merge_dict[k][12] + "\t" + merge_dict[k][13]) # rps evalue + rps taxo
		f.write("\n")
	f.close()


def _set_options():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--blastx',help='Csv files with blastx results.', action='append',required=True, type=argparse.FileType('r'),dest='bx_files')
	parser.add_argument('-p','--pfam',help='Csv files with rpsblast results.', action='append',required=True, type=argparse.FileType('r'),dest='pfam_files')
	parser.add_argument('-r','--rps',help='RPS results folder', action='store', required=True, type=str,dest='rps_folder')
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