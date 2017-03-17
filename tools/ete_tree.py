#!/home/stheil/softwares/miniconda2/bin/python
import sys
from optparse import OptionParser
from ete3 import PhyloTree, Tree, NodeStyle, TreeStyle, SeqGroup, SeqMotifFace
import csv

def main ():
    options = _set_options()
    _create_tree(options.tree_file,options.fasta_file,options.out_file,options.color_file)

def _create_tree (tree,fasta,out,color):
    seqs = SeqGroup(fasta, format="fasta")
    t = Tree(tree)
    colors = _parse_color_file(color)
    node_names = t.get_leaf_names()
    for name in node_names:
        seq = seqs.get_seq(name)
        seqFace = SeqMotifFace(seq, seq_format="()")
        node = t.get_leaves_by_name(name)
        for i in range(0,len(node)):
            if name in colors:
                ns = NodeStyle()
                ns['bgcolor'] = colors[name]
                node[i].set_style(ns)
            node[i].add_face(seqFace,0,'aligned')
    t.render(out)

def _parse_color_file (file):
    fh = open(file)
    reader = csv.reader(fh,delimiter="\t")
    data = list(reader)
    colors = {}
    for i in range(0,len(data)):
        colors[data[i][0]] = data[i][1]
    return colors

def _set_options ():
	usage = "usage: %prog -t newick.tree -f seq.fasta -c color.txt"
	parser = OptionParser(usage)
	parser.add_option('-t','--tree',dest='tree_file',help='A newick file containing the tree.',action='store',metavar='FILE')
	parser.add_option('-f','--fasta',dest='fasta_file',help='A fasta file containing sequences to align.',action='store',metavar='FILE')
	parser.add_option('-o','--out',dest='out_file',help='A fasta file containing sequences to align.',action='store',metavar='FILE')
	parser.add_option('-c','--color',dest='color_file',help='A csv file containing the coloring scheme for the tree.',action='store',metavar='FILE')
	(options, args) = parser.parse_args()
	if not options.fasta_file:
		parser.error('You must provide a fasta file.')
	return options


def _help ():
	print 'pfam2tree.py -f seq.fasta -c color.txt'


if __name__ == "__main__":
	main()
