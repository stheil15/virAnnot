# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path

class Select:
	"""
	Generate fasta file from blastx results
	With viral contig only
	"""

	def __init__ (self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd (self):
		cmd = 'select.py'
		cmd += ' -b ' + self.blastx
		cmd += ' -o ' + self.out
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		self.execution = 1
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_select_cmd.txt'
		if 'blastx' in args:
			self.blastx = self.check_file(self.wd + '/' + args['blastx'])
		else:
			self.execution = 0
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
		else:
			self.out = os.getcwd() + '/' + self.sample + '_viral-contig.fa'
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'out' in args:
			self.out = args['out']

	def check_file (f):
		try:
			open(f)
			return f
		except IOError:
			print('File not found ' + f)

