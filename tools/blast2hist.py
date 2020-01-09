import sys
import csv, re, os
import argparse
import logging as log
import matplotlib.pyplot
import numpy as np
from matplotlib import colors as mcolors
import random


def main ():
	options = _set_options()
	_set_log_level(options.verbosity)
	all_data = {}
	ranks = ('superkingdom','virus_type','family','genus','specie')
	for i in range(0, len(options.blast_list)):
		d,h = _read_blast_file(options.blast_list[i])
		for j in range(0,len(ranks)):
			log.info('sorting by rank ' + ranks[j])
			_get_rank_data(all_data,j,d,h,i,len(options.id_list))
	if not os.path.exists(options.out):
		os.makedirs(options.out)
	plot_files={}
	for r in range(0,len(ranks)):
		if ranks[r] not in plot_files:
			plot_files[ranks[r]]={}
		_get_rank_hist(all_data[r],ranks[r],options.id_list,options.out,plot_files)
	_create_html(plot_files,options.out,ranks)

def _create_html(plot_files,out_dir,ranks):
	log.info('creating html page.')
	fh = open(out_dir + '/index.html', mode='w')
	html = '<html>' + "\n"
	html += '<head>' + "\n"
	html += '<title>' + 'Blast histograms' + '</title>'
	html += '</head>' + "\n"
	html += '<div style="text-align:center">' + "\n"
	html += '<h1 align=center>Blast histograms</h1>' + "\n"
	html += '<body>' + "\n"
	for rank in ranks:
		html += '<h2>' + rank + '</h2>' + "\n"
		html += '</br>' + "\n"
		for key in plot_files[rank]:
			html += '<h3>' + key + '</h3>' + "\n"
			html += '</br>' + "\n"
			html += '<img width="100%" align="center" src="' + plot_files[rank][key]['plot'] + '">' + "\n"
	html += '</body>' + "\n"
	html += '</div>' + "\n"
	html += '</html>' + "\n"
	fh.write(html)
	fh.close()


def _get_rank_hist (data,rank,id_list,out_dir,plot_files):
	log.info('getting histogram for rank ' + rank)
	if not os.path.exists(out_dir + '/data'):
		os.makedirs(out_dir + '/data')
	width=0.4
	ind = np.arange(len(id_list))
	colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
	colors_names = list(colors.keys())
	random.shuffle(colors_names)
	keys = ['contigs','reads','cumul_length']
	for key in keys:
		log.info('getting histogram for key ' + key)
		if key not in plot_files:
			plot_files[rank][key]={}
		plt = matplotlib.pyplot
		plt.cla()
		plt.clf()
		plt.figure(figsize=[30,21], dpi=400)

		matplotlib.rcParams.update({'font.size': 16})
		plt.gca().yaxis.grid(True)
		plt.gcf().subplots_adjust(bottom=0.15)
		x_names=[]
		handles=[]
		tax_id = sorted(list(data.keys()))
		log.debug('treating ' + str(len(tax_id)) + ' tax_id')
		matrix=[]
		for j in range(0,len(tax_id)):
			x_names.append(tax_id[j])
			matrix.append(data[tax_id[j]][key])
		cumul = _get_cumul_matrix(matrix)
		ax1 = plt.subplot2grid((7,4), (0,0), colspan=3, rowspan=3)
		ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
		picked_colors={}
		i=0
		log.debug('creating raw numbers histograms')
		for j in range(0,len(tax_id)):
			if i >= len(colors_names):
				i=0
				picked_colors[tax_id[j]] = colors[colors_names[i]]
			else:
				picked_colors[tax_id[j]] = colors[colors_names[i]]
			if j == 0:
				p = plt.bar(ind,data[tax_id[j]][key],width,color=colors[colors_names[i]],label=tax_id[j])
			else:
				p = plt.bar(ind,data[tax_id[j]][key],width,color=colors[colors_names[i]],label=tax_id[j],bottom=cumul[j-1])
			handles.append(p[0])
			i+=1
		plt.ylabel('Number of ' + key + ' per ' + rank)
		# plt.xticks(ind, id_list, rotation='vertical')
		plt.legend(handles,x_names,loc=(1.05,-1.75),ncol=1)
		# plt.subplot(2,1,2)
		percent = _get_percent_matrix(cumul,matrix)
		cumul_percent = _get_cumul_matrix(percent)
		log.debug('creating percent histograms')
		for j in range(0,len(tax_id)):
			if j == 0:
				p = plt.bar(ind,percent[j],width,color=picked_colors[tax_id[j]],label=tax_id[j])
			else:
				p = plt.bar(ind,percent[j],width,color=picked_colors[tax_id[j]],label=tax_id[j],bottom=cumul_percent[j-1])
		plt.ylabel('Percentage of ' + key + ' per ' + rank)
		plt.xticks(ind, id_list, rotation='vertical')
		plt.savefig(out_dir + '/data/' + key + '_' + rank + ".svg")
		plt.close('all')
		plot_files[rank][key]['plot'] = 'data/' + key + '_' + rank + ".svg"
		plot_files[rank][key]['matrix'] = matrix
		plot_files[rank][key]['percent'] = percent
		plot_files[rank][key]['headers'] = x_names


def _get_percent_matrix (cumul,matrix):
	percent=[]
	for i in range(0,len(matrix)):
		percent.append([])
		for j in range(0, len(matrix[i])):
			if cumul[len(cumul)-1][j] == 0:
				percent[i].append(0)
			else:
				p = matrix[i][j] / cumul[len(cumul)-1][j]
				percent[i].append(p)
	return percent


def _get_cumul_matrix (matrix):
	cumul = []
	for i in range(0,len(matrix)):
		if i == 0:
			cumul.append(matrix[i])
		else:
			cumul.append([])
			for j in range(0,len(matrix[i])):
				cumul[i].append(matrix[i][j] + cumul[i-1][j])
	return cumul


def _get_rank_data (all_data,rank,data,headers,i,id_list_len):
	if rank not in all_data:
		all_data[rank]= {}
	for el in data:
		if(re.search('Viruses',el['taxonomy']) or re.search('Viroids',el['taxonomy'])):
			tax = el['taxonomy'].split(';')[rank]
			if tax not in all_data[rank]:
				all_data[rank][tax] = {}
				if 'contigs' not in all_data[rank][tax]:
					all_data[rank][tax]['contigs'] = [0] * id_list_len
				all_data[rank][tax]['contigs'][i] += 1
				if 'reads' not in all_data[rank][tax]:
					all_data[rank][tax]['reads'] = [0] * id_list_len
				if 'nb_reads' in headers:
					if el['nb_reads'] != "":
						all_data[rank][tax]['reads'][i] += int(el['nb_reads'])
				if 'cumul_length' not in all_data[rank][tax]:
					all_data[rank][tax]['cumul_length'] = [0] * id_list_len
				if 'query_length' in headers:
					all_data[rank][tax]['cumul_length'][i] += int(el['query_length'])
			else:
				if 'nb_reads' in headers:
					if el['nb_reads'] != "":
						all_data[rank][tax]['reads'][i] += int(el['nb_reads'])
				if 'query_length' in headers:
					all_data[rank][tax]['cumul_length'][i] += int(el['query_length'])
				all_data[rank][tax]['contigs'][i] += 1
	return all_data


def _read_blast_file (file):
	log.info('reading ' + file.name)
	reader = csv.reader(file,delimiter="\t")
	data = list(reader)
	headers = data[0]
	headers[0] = headers[0][1:]
	map_obj = []
	for i in range(1,len(data)):
		dict={}
		if len(data[i]) != len(headers):
			print(data[i])
			sys.exit('line and headers not the same length.')
		for j in range(0,len(headers)):
			dict[headers[j]] = data[i][j]
		map_obj.append(dict)
	return map_obj,headers


def _set_options ():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--blast',help='A csv Blast file.',action='append',required=True,type=argparse.FileType('r'),dest='blast_list')
	parser.add_argument('-i','--id',help='The ID corresponding to the Blast file..',action='append',required=True,dest='id_list')
	parser.add_argument('-o','--out',help='The output the HTML page.',action='store',type=str,default='./')
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
