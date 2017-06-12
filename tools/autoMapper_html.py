#!/usr/bin/python
import sys, getopt
import csv
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main ():
    options = _set_options()
    data,headers = _read_map_file(options.map)
    html = _print_html(options.title,data,headers)
    index_file = options.out + '/index.html'
    fh = open(index_file, mode='w')
    fh.write(html)
    fh.close


def _print_html (title,data,headers):
    html = '<html>' + "\n"
    html += '<head>' + "\n"
    html += '<title>' + title + '</title>' + "\n"
    html += '</head>'
    html += '<div style="text-align:center"' + "\n"
    html += '<h1 align=center>' + title + '</h1>' + "\n"
    html += '<body>' + "\n"
    html += _print_data(data,headers)
    html += '</body>' + "\n"
    html += '</div>'
    html += '</html>' + "\n"
    return html


def _print_data (data,headers):
    html=''
    for ref in data:
        html += '<h2 align=center>' + ref['organism'].replace('_', ' ') + '</h2>' + "\n"
        html += '<p>tax_id: ' + ref['tax_id'] + '</p>' + "\n"
        html += '<p>taxonomy: ' + ref['taxonomy'] + '</p>' + "\n"
        for s in ref['coords_file'].split(','):
            html += '<a href="' + s + '">' + s + '</a>' + '</br>' + "\n"
        for s in ref['fasta_file'].split(','):
            html += '<a href="' + s + '">' + s + '</a>' + '</br>' + "\n"
        for s in ref['coverage'].split(','):
            html += '<p>' + s + '</p>' + "\n"
        html += '<img src=' + ref['svg_file'] + ' style="width:100%;height:80%">' + "\n"
        html += "\n"
    return html


def _read_map_file (file):
	reader = csv.reader(file,delimiter="\t")
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
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--map',help='The map file produced by autoMapper.',action='store',type=argparse.FileType('r'),required=True)
    parser.add_argument('-t','--title',help='The title for the HTML page.',action='store',type=str,required=True)
    parser.add_argument('-o','--out',help='The title for the HTML page.',action='store',type=str,required=True)
    args = parser.parse_args()
    return args



if __name__ == "__main__":
	main()
