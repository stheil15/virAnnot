# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys


class Filter:
	"""
    Filter step of virAnnot module
	From assembly scaffold fasta output, filter sequences over a given length threshold
	Author: Eden Darnige
	"""

	def __init__(self, args):
		self.execution = 1
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd(self):
		"""
		Create command
		"""
		cmd = 'filter.py'
		cmd += ' -i ' + self.i
		cmd += ' -l ' + self.len
		cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	# Verify that all mandatory parameters are present
	def check_args(self, args:dict):
		self.execution = 1
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_filter_cmd.txt'
		if 'i' in args:
			self.i = self._check_file(self.wd + '/' + args['i'])
		else:
			self.i = ''
		if 'len' in args:
			self.len = str(args['len'])
		else:
			self.len = '1'
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
        #keep
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
			if 'sge' in args:
			    self.sge = bool(args['sge'])


	# Exisiting file
	def _check_file(self, input_file):
		try:
			open(input_file)
			return input_file
		except IOError:
			print('File not found ' + input_file)
			self.execution = 0
