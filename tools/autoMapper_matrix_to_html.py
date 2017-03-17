#!/usr/bin/python
import sys, getopt
import gviz_api
import csv
from optparse import OptionParser


def _print_html (json_obj):
	html = '<html>' + "\n" \
	+ '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>' + "\n" \
	+ '<script>' + "\n" \
	+ 'google.charts.load("current", {packages: ["corechart","table"]});' + "\n" \
	+ 'google.charts.setOnLoadCallback(drawChart);' + "\n" \
	+ "\n" \
	+ 'function drawChart() {' + "\n"
	js_table, html_div_table = _print_json_data_table(json_obj)
	for i in range(0,len(js_table)):
		html += js_table[i]
	html += '}' + "\n" \
	+ '</script>' + "\n" \
	+ '<body>' + "\n"
	for i in range(0,len(html_div_table)):
		if i%2==0:
			html += '<div style="width: 100%; height:600px; display: table;">' + "\n"
		html += html_div_table[i]
		if i%2!=0:
			html += '</div>' + "\n"
	html += '</body>' + "\n" \
	+ '</html>' + "\n" \
	+ '' + "\n"
	return html

def _print_json_data_table (json_obj):
	js_table = []
	html_div_table = []
	aln_link_table = []
	nb_graph=0
	for i in range(0,len(json_obj)):
		if json_obj[i].get('json_scaff'):
			nb_graph +=1
			js_str = 'var json_scaff_data_' + str(i) + ' = new google.visualization.DataTable(' + json_obj[i]['json_scaff'] + ',0.6);' + "\n"
			js_str += 'var chart_' + str(i) + ' = new google.visualization.LineChart(document.getElementById(\'curve_chart_' + str(i) + '\'));' + "\n"
			js_str += 'var options_scaff_' + str(i) + ' = {title: "Scaffolds identity plot", width: "100%", height: 600, chartArea: { width: "70%", left: 50, top: 30}, legend: {position: "right"}, vAxis: {minValue: 0}, hAxis: {minValue: 0}};' + "\n"
			js_str += 'chart_' + str(i) + '.draw(json_scaff_data_' + str(i) + ', options_scaff_' + str(i) + ');' + "\n"
			js_table.append(js_str)
			html_div_str = '<h3>' + json_obj[i]['organism'] + '</h3>' + "\n"
			html_div_str += '<p>' + json_obj[i]['taxonomy'] + '<br>' + '</p>' + "\n"
			html_div_str += '<p>' + json_obj[i]['description'] + '<br>' + '</p>' + "\n"
			html_div_str += '<a href="./' + str(json_obj[i]['fasta_file']) + '"target="_blank">Alignment file</a><br>' + "\n"
			html_div_str += '<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + str(json_obj[i]['accession']) + '"target="_blank">NCBI accession link</a><br>' + "\n"
			html_div_str += '<a href="http://www.ncbi.nlm.nih.gov/taxonomy/' + str(json_obj[i]['tax_id']) + '"target="_blank">NCBI Taxonomy link</a><br>' + "\n"
			html_div_str += '<div id="curve_chart_' + str(i) + '" style="width: "100%"; "></div>' + "\n"

			html_div_table.append(html_div_str)
		if json_obj[i].get('json_read'):
			nb_graph +=1
			js_str = 'var json_read_data_' + str(i) + ' = new google.visualization.DataTable(' + json_obj[i]['json_read'] + ',0.6);' + "\n"
			js_str += 'var scatter_' + str(i) + ' = new google.visualization.AreaChart(document.getElementById(\'scatter_plot_' + str(i) + '\'));' + "\n"
			js_str += 'var options_read_' + str(i) + ' = {title: "Read coverage", width: "100%", height: 600, chartArea: {top: 30, width: "70%", left: 50}, legend: {position: "none"}, vAxis: {minValue: 0}, hAxis: {minValue: 0}};' + "\n"
			js_str += 'scatter_' + str(i) + '.draw(json_read_data_' + str(i) + ', options_read_' + str(i) + ');' + "\n"
			js_table.append(js_str)
			html_div_str = '<div id="scatter_plot_' + str(i) + '" style="width: "100%"; "></div>' + "\n"
			html_div_table.append(html_div_str)
		if len(html_div_table)%2 !=0:
			html_div_table.append('')
		html_div_table.append('<hr>')
	return js_table, html_div_table

def main ():
	options = _set_options()
	map_obj, map_headers = _read_map_file(options.map)

	for i in range(0,len(map_obj)):
		map_obj[i]['json_scaff']=[]
		map_obj[i]['json_read']=[]
		if(map_obj[i]['mat_file'] != '.'):
			desc, data, headers = _parse_csv(map_obj[i]['mat_file'])
			data_table = gviz_api.DataTable(desc)
			data_table.LoadData(data)
			# map_obj[i]['json_scaff']=[]
			map_obj[i]['json_scaff'] = data_table.ToJSon(columns_order=headers)
		if(map_obj[i]['read_map_file'] != '.'):
			desc, data, headers = _parse_reads_matrix(map_obj[i]['read_map_file'])
			data_table = gviz_api.DataTable(desc)
			data_table.LoadData(data)
			# map_obj[i]['json_read']=[]
			map_obj[i]['json_read'] = data_table.ToJSon(columns_order=headers)
	html = _print_html(map_obj)
	print html

def _parse_reads_matrix (file):
	fh = open(file)
	reader = csv.reader(fh,delimiter=';')
	data = list(reader)
	headers = data[0]
	headers[0] = headers[0][1:]
	description = {}
	for i in range(len(headers)):
		description[headers[i]] = ("number", headers[i])
	matrix = []
	for i in range(1,len(data)):
		dict={}
		for j in range(0,len(data[i])):
			if(data[i][j] == 'null'):
				dict[headers[j]] = None
			elif(data[i][j] == ''):
				dict[headers[j]] = None
			else:
				dict[headers[j]] = float(data[i][j])
		matrix.append(dict)
	return description, matrix, headers


def _read_map_file (file):
	fh = open(file)
	reader = csv.reader(fh,delimiter="\t")
	data = list(reader)
	headers = data[0]
	headers[0] = headers[0][1:]
	map_obj = []
	for i in range(1,len(data)):
		dict={}
		if len(data[i]) != len(headers):
			sys.exit('line and headers not the same length.')
		for j in range(0,len(headers)):
			dict[headers[j]] = data[i][j]
		map_obj.append(dict)
	return map_obj,headers


def _set_options ():
	usage = "usage: %prog -m map.txt"
	parser = OptionParser(usage)
	parser.add_option('-m','--map',dest='map',help='A csv file mapping the autoMapper.pl informations. #tax_id,accession,organism,description,fasta_file,aln_file,read_file,ctg_file',action='store',metavar='FILE')
	(options, args) = parser.parse_args()
	if not options.map:
		parser.error('You must provide a map file.')
	return options


def _parse_csv (file):
	fh = open(file)
	reader = csv.reader(fh,delimiter=";")
	data = list(reader)
	headers = data[0]
	description = {}
	for i in range(len(headers)):
		description[headers[i]] = ("number", headers[i])
	matrix = []
	for i in range(1,len(data)):
		dict={}
		for j in range(0,len(data[i])):
			if(data[i][j] == ''):
				dict[headers[j]] = None
			elif(data[i][j] == 'null'):
				dict[headers[j]] = None
			else:
				dict[headers[j]] = float(data[i][j])
		matrix.append(dict)

	# print description
	# print matrix
	return description,matrix,headers

def _help ():
	print 'autoMapper_matrix_to_html.py -m matrix.txt'

if __name__ == "__main__":
	main()
