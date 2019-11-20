import os.path
from subprocess import call
import logging as log
from collections import defaultdict
import csv

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


	def check_args(self, args:dict):
		self.execution = 1
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_select_cmd.txt'
		if 'blastx' in args:
			self.blastx = self._check_file(self.wd + '/' + args['blastx'])
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

	def _check_file (self,f):
		try:
			open(f)
			return f
		except IOError:
			print('File not found ' + f)

