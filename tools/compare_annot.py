import sys
import csv, re, os
import argparse
import logging as log
from matplotlib import rcParams
rcParams.update({'figure.autolayout': False})


log.basicConfig(level=log.INFO)
log = log.getLogger(__name__)



def main ():
	options = _set_options()
	raw_data=[]
	raw_headers=[]
	for i in range(0, len(options.blast_list)):
		m,h = _read_blast_file(options.blast_list[i])
		raw_data.append(m)
		raw_headers.append(h)
	if not os.path.exists(options.out):
		os.makedirs(options.out)
	all_data = {}
	for i in range(0, len(options.blast_list)):
		_get_seq_data(all_data,raw_data[i],raw_headers[i],i,options.id_list,options.regex_list)
	_print_csv_file(all_data, options.id_list, options.out, raw_headers)



def _print_csv_file (all_data, id_list, out_dir, raw_headers):
	log.info('creating CSV')
	filename = out_dir + '/' + 'blast_compare.csv'
	csvfile = open(filename, 'w', newline='')
	csv_writer = csv.writer(csvfile, delimiter="\t",quotechar='|', quoting=csv.QUOTE_MINIMAL)
	csv_writer.writerow(['databank'] + raw_headers[0])
	for seq_id in sorted(all_data.keys(), key = lambda x: int(x.split('_')[1])):
		csv_writer.writerow([seq_id])
		for i in range(0,len(id_list)):
			line = []
			line.append(id_list[i])
			if id_list[i] in all_data[seq_id]:
				for h in raw_headers[i]:
					line.append(all_data[seq_id][id_list[i]][h])
			else:
				line.append('None')
			csv_writer.writerow(line)



def _get_seq_data (all_data,data,headers,i,id_list,regex_list):
	for el in data:
		if el['accession'] != '':
			if regex_list:
				found=False
				if len(regex_list) >=1:
					for word in regex_list:
						if(re.search(word,el['taxonomy'])):
							found=True
				else:
					found=True
			else:
				found=True
			if found:
				if el['taxonomy'] != 'unknown':
					tax = el['taxonomy'].split(';')
				else:
					tax = '_undef_'
				if el['query_id'] not in all_data:
					all_data[el['query_id']] = {}
				if id_list[i] not in all_data[el['query_id']]:
					all_data[el['query_id']][id_list[i]] = []
				all_data[el['query_id']][id_list[i]] = el


def _read_blast_file (file):
	log.info('reading csv file ' + file.name)
	reader = csv.reader(file,delimiter="\t")
	data = list(reader)
	h = data[0]
	h[0] = h[0][1:]
	map_obj = []
	for i in range(1,len(data)):
		dict={}
		if len(data[i]) != len(h):
			sys.exit('line and h not the same length.')
		for j in range(0,len(h)):
			dict[h[j]] = data[i][j]
		map_obj.append(dict)

	return map_obj, h


def _set_options ():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--blast',help='A csv Blast file.',action='append',required=True,type=argparse.FileType('r'),dest='blast_list')
	parser.add_argument('-i','--id',help='The ID corresponding to the Blast file..',action='append',required=True,dest='id_list')
	parser.add_argument('-o','--out',help='The output the HTML page.',action='store',type=str,default='./')
	parser.add_argument('--regex',help='Limit matches the provided regex.',action='append',dest='regex_list',default=None)
	parser.add_argument('--levels',help='Maximum number of taxa level',action='store',type=int,dest='max_level',default=4)
	args = parser.parse_args()
	return args



if __name__ == "__main__":
	main()
