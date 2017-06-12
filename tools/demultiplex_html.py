#!/usr/bin/python
import os
import csv
import argparse
import logging


def main ():
    options = _set_options()
    _create_html(options)



def _create_html(options):
    if not os.path.exists(options.out):
        os.mkdir(options.out)
    fh = open(options.out + '/index.html', mode='w')
    html = '<html>' + "\n"
    html += '<head>' + "\n"
    html += '<title>' + 'Demultiplex statistics' + '</title>' + "\n"
    html += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>' + "\n"
    html += '<script type="text/javascript">' + "\n"
    html += 'google.charts.load(\'current\', {\'packages\':[\'table\']});' + "\n"
    array_names, tmp = _get_google_scripts(options)
    html += tmp
    html += '</head>' + "\n"
    html += '<div style="text-align:center">' + "\n"
    html += '<h1 align=center>Demultiplex statistics</h1>' + "\n"
    html += '<body>' + "\n"
    for i in range(0,len(array_names)):
        html += '<h2>' + array_names[i] + '</h2>' + "\n"
        html += '<div id="' + array_names[i] + '_div' + '"></div>' + "\n"

    html += '</body>' + "\n"
    html += '</div>' + "\n"
    html += '</html>' + "\n"
    fh.write(html)
    fh.close()


def _get_google_scripts (options) :
    array_names, java_scripts = _get_google_js(options)
    html=''
    for i in range(0, len(array_names)):
        html += 'google.charts.setOnLoadCallback(' + array_names[i].replace('-', '_') + ');' + "\n"
        html += 'function ' + array_names[i].replace('-', '_') + '() {'  + "\n"
        html += java_scripts[i] + "\n"
        html += '}' + "\n"
    html += '</script>' + "\n"
    return array_names, html

def _get_google_js (options):
    java_scripts = []
    array_names = []
    for i in range(0, len(options.csv_list)):
        _read_csv_file(options.csv_list[i], options.id_list[i],java_scripts,array_names)
    return array_names, java_scripts

def _read_csv_file (file,id,java_scripts,array_names):
    reader = csv.reader(file,delimiter=",")
    csv_data = list(reader)
    cmp=0
    for i in range(0,len(csv_data)):
        if len(csv_data[i]) == 0:
            js = js[:-1]
            js += ']);'  + "\n"
            js += 'var table = new google.visualization.Table(document.getElementById(\'' + (id + '_' + str(cmp)) + '_div' + '\'));' + "\n"
            js += 'table.draw(data, {showRowNumber: false, width: \'70%\', height: \'70%\'});' + "\n"
            java_scripts.append(js)
            array_names.append(id + '_' + str(cmp))
            cmp+=1
            continue
        if csv_data[i][0].startswith('#'):
            js = 'var data = new google.visualization.DataTable();' + "\n"
            for el in csv_data[i]:
                if el == '#step':
                    js += 'data.addColumn(\'string\', \'' + el + '\');' + "\n"
                elif el == 'file':
                    js += 'data.addColumn(\'string\', \'' + el + '\');' + "\n"
                elif el == '#index_id':
                    js += 'data.addColumn(\'string\', \'' + el + '\');' + "\n"
                else:
                    js += 'data.addColumn(\'number\', \'' + el + '\');' + "\n"
            js += 'data.addRows([' + "\n"
        else:
            js += '['
            for j in range(0,len(csv_data[i])):
                if(j == len(csv_data[i])-1):
                    if csv_data[i][j].isdigit():
                        js += ',' + csv_data[i][j]
                    else:
                        js += ',\'' + csv_data[i][j] + '\''
                    js += '],'
                elif(j == 0):
                    if csv_data[i][j].isdigit():
                        js += csv_data[i][j]
                    else:
                        js += '\'' + csv_data[i][j] + '\''
                else:
                    if csv_data[i][j].isdigit():
                        js += ',' + csv_data[i][j]
                    else:
                        js += ',\'' + csv_data[i][j] + '\''
    js = js[:-1]
    js += ']);'  + "\n"
    js += 'var table = new google.visualization.Table(document.getElementById(\'' + (id + '_' + str(cmp)) + '_div' + '\'));' + "\n"
    js += 'table.draw(data, {showRowNumber: false, width: \'70%\', height: \'70%\'});' + "\n"
    java_scripts.append(js)
    array_names.append(id + '_' + str(cmp))

def _set_options ():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--csv',help='A csv file from demultiplex.pl.',action='append',required=True,type=argparse.FileType('r'),dest='csv_list')
    parser.add_argument('-i','--id',help='The ID corresponding to the CSV file..',action='append',required=True,dest='id_list')
    parser.add_argument('-o','--out',help='The output the HTML page.',action='store',type=str,default='./')
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    main()
