import sys
import csv, re, os
import argparse
import logging as log
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
rcParams.update({'figure.autolayout': False})

log.basicConfig(level=log.INFO)
log = log.getLogger(__name__)


def main ():
	options = _set_options()
	keys = []
	levels = range(0,options.max_level)
	# for lvl in levels:
	raw_data=[]
	raw_headers=[]
	for i in range(0, len(options.blast_list)):
		m,h = _read_blast_file(options.blast_list[i])
		raw_data.append(m)
		raw_headers.append(h)
	if not os.path.exists(options.out):
		os.makedirs(options.out)
	for lvl in levels:
		all_data = {}
		for i in range(0, len(options.blast_list)):
			_get_rank_data(all_data,raw_data[i],raw_headers[i],i,len(options.id_list),keys,lvl,options.regex_list)
		for k in sorted(keys):
			_print_csv_files(all_data,options.id_list, k, lvl, options.out)
			_plot_heatmap(all_data, options.id_list, k, lvl, options.out)


def _plot_heatmap (all_data, id_list, key, lvl, out_dir):
	log.info('creating heatmap for ' + key + ' lvl ' + str(lvl))
	filename = out_dir + '/' + 'blast_heatmap_' + key + '_l' + str(lvl) + '.svg'
	data = []
	row_labels = all_data.keys()
	for tax in sorted(row_labels):
		if key == 'identity':
			mean_list = []
			for i in range(0,len(all_data[tax][key])):
				if len(all_data[tax][key][i]) != 0:
						mean_list.append(np.mean(all_data[tax][key][i]))
				else:
					mean_list.append(0)
			data.append(mean_list)
		else:
			data.append(all_data[tax][key])
	np_data = np.array(data)
	axis = plt.subplot(111)
	heatmap = axis.pcolor(np_data, cmap=plt.cm.Greens)
	axis.set_yticks(np.arange(np_data.shape[0])+0.5, minor=False)
	axis.set_xticks(np.arange(np_data.shape[1])+0.5, minor=False)
	axis.set_yticklabels(sorted(row_labels), minor=False)
	axis.set_xticklabels(id_list, minor=False, rotation='vertical')

	fig = plt.gcf()# plt.subplots_adjust(left=0.2,bottom=0.3,top=0.7)
	fig.set_size_inches(20,20)
	plt.tight_layout(w_pad=1.0, h_pad=1.0)
	plt.savefig(filename,dpi=200)
	plt.close('all')


def _print_csv_files (all_data, id_list, key, lvl, out_dir):
	log.info('creating CSV for ' + key + ' lvl ' + str(lvl))
	filename = out_dir + '/' + 'blast_matrix_' + key + '_l' + str(lvl) + '.csv'
	csvfile = open(filename, 'w', newline='')
	csv_writer = csv.writer(csvfile, delimiter="\t",quotechar='|', quoting=csv.QUOTE_MINIMAL)
	csv_writer.writerow(id_list + ['taxonomy'])
	row_labels = all_data.keys()
	for tax in sorted(row_labels):
		if key == 'identity':
			mean_list = []
			for i in range(0,len(all_data[tax][key])):
				if len(all_data[tax][key][i]) != 0:
					mean_list.append(np.mean(all_data[tax][key][i]))
				else:
					mean_list.append(0)
			csv_writer.writerow(mean_list + [tax])
		else:
			csv_writer.writerow(all_data[tax][key] + [tax])



def _get_rank_data (all_data,data,headers,i,id_list_len,keys,lvl,regex_list):
	for el in data:
		if el['accession'] != '':
			found=False
			if len(regex_list) >=1:
				for word in regex_list:
					if(re.search(word,el['taxonomy'])):
						found=True
			else:
				found=True
			if found:
				if el['taxonomy'] != 'unknown':
					tax = el['taxonomy'].split(';')[lvl]
				else:
					tax = '_undef_'
				if tax not in all_data:
					all_data[tax] = {}
					if 'contigs' not in all_data[tax]:
						if 'contigs' not in keys:
							keys.append('contigs')
						all_data[tax]['contigs'] = [0] * id_list_len
					all_data[tax]['contigs'][i] += 1
					if 'nb_reads' in headers:
						if 'reads' not in all_data[tax]:
							all_data[tax]['reads'] = [0] * id_list_len
						if 'reads' not in keys:
							keys.append('reads')
						if el['nb_reads'] != "":
							all_data[tax]['reads'][i] += int(el['nb_reads'])
					if 'query_length' in headers:
						if 'cumul_length' not in all_data[tax]:
							all_data[tax]['cumul_length'] = [0] * id_list_len
						if 'cumul_length' not in keys:
							keys.append('cumul_length')
						if el['query_length'] != "":
							all_data[tax]['cumul_length'][i] += int(el['query_length'])
					if 'percentIdentity' in headers:
						if 'identity' not in all_data[tax]:
							all_data[tax]['identity'] = [[] for j in range(id_list_len)]
						if 'identity' not in keys:
							keys.append('identity')
						all_data[tax]['identity'][i].append(float(el['percentIdentity']))
				else:
					if 'nb_reads' in headers:
						if el['nb_reads'] != "":
							all_data[tax]['reads'][i] += int(el['nb_reads'])
					if 'query_length' in headers:
						if el['query_length'] != "":
							all_data[tax]['cumul_length'][i] += int(el['query_length'])
					if 'percentIdentity' in headers:
						all_data[tax]['identity'][i].append(float(el['percentIdentity']))
					all_data[tax]['contigs'][i] += 1


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
