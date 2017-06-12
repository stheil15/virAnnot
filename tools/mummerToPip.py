#!/usr/bin/python
import sys
from optparse import OptionParser
import csv
import matplotlib.pyplot as plt
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main ():
    options = _set_options()
    data = dict()
    data['sample'] = _read_tab_file(options.input_file,options.id)
    _create_plot(data['sample'],options.prefix,options.output,options.color,options.plot_title)


def _read_tab_file (file,identity):
    fh = open(file, "r")
    lines = fh.readlines()
    data = dict()
    print (lines[1].rstrip())
    if(lines[1].rstrip() == 'PROMER'):
        data = _parse_promer(lines,identity)
    elif(lines[1].rstrip() == 'NUCMER'):
        data = _parse_nucmer(lines,identity)
    else:
        sys.exit(1)
    return data


def _parse_promer(lines,identity):
    data = dict()
    for i in range(4,len(lines)):
        line = lines[i].split("\t")
        if(float(line[6]) < identity):
            continue
        if(line[13] in data):
            data[line[13]]['coords'].append(_create_line([line[0],line[1],line[6]]))
        else:
            data[line[13]] = dict()
            data[line[13]]['coords'] = []
            data[line[13]]['length'] = int(line[9])
            data[line[13]]['coords'].append(_create_line([line[0],line[1],line[6]]))
    return data


def _parse_nucmer(lines,identity):
    data = dict()
    for i in range(4,len(lines)):
        line = lines[i].split("\t")
        if(float(line[6]) < identity):
            continue
        if(line[9] in data):
            data[line[9]]['coords'].append(_create_line([line[0],line[1],line[6]]))
        else:
            data[line[9]] = dict()
            data[line[9]]['coords'] = []
            data[line[9]]['length'] = int(line[7])
            data[line[9]]['coords'].append(_create_line([line[0],line[1],line[6]]))
    return data


def _create_plot (data,prefix,output,color,title):
    total_length = 0
    plt.figure(figsize=[18,12], dpi=400)
    plt.xlabel('Reference in bp')
    plt.ylabel('% identity')
    for key in data:
        print(key)
        for i in range(len(data[key]['coords'])):
            plt.plot([data[key]['coords'][i][0][0]+total_length, data[key]['coords'][i][0][1]+total_length], data[key]['coords'][i][1],color=color)
        total_length += data[key]['length']
        plt.axvline(x=total_length,color='grey')
    plt.axis([0,total_length,0,101])
    plt.savefig(output + ".svg")
    plt.close


def _create_line (line):
    curve = []
    x = []
    x.append(int(line[0]))
    x.append(int(line[1]))
    curve.append(x)
    y = []
    y.append(float(line[2]))
    y.append(float(line[2]))
    curve.append(y)
    return curve


def _set_options ():
    usage = "usage: %prog -t mummer.tab.txt"
    parser = OptionParser(usage)
    parser.add_option('-p','--prefix',dest='prefix',help='A prefix for the output files.',action='store',metavar='FILE',default='pip')
    parser.add_option('-i','--input',dest='input_file',help='A tab file produce by mummer.',action='store',metavar='FILE')
    parser.add_option('--id',dest='id',help='Minimum identity value 0 < X < 100.',action='store',metavar='float',default=0)
    parser.add_option('-o','--output',dest='output',help='An output file name.',action='store',metavar='string')
    parser.add_option('-c','--color',dest='color',help='A color name.',action='store',metavar='string',default='blue')
    parser.add_option('-t','--title',dest='plot_title',help='A plot title. If spaces use "".',action='store',metavar='string',default='Percent identity plot')
    (options, args) = parser.parse_args()
    return options


def _help ():
    print 'mummerToPip.py -t mummer.tab.txt'


if __name__ == "__main__":
    main()
