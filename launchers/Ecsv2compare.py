# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys

class Ecsv2compare:

	def __init__ (self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()

	def create_cmd (self):
		cmd = 'ecsv2compare.py'
		for c in self.blast_files:
			cmd += ' -b ' + str(c)
		if self.rps_file != '':
			cmd += ' -r ' + self.rps_file
		cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args (self, args=dict):
		self.execution=1
		self.sample = args['sample']
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + 'ecsv2compare_cmd.txt'
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
		self.blast_files = []
		for i in range(1, 10, 1):
			opt_name = 'b' + str(object=i)
			if opt_name in args:
				if os.path.exists(self.wd + '/' + args[opt_name]):
					self.blast_files.append(self.wd + '/' + args[opt_name])
		if 'r' in args:
			if os.path.exists(self.wd + '/' + args['r']):
				self.rps_file = self._check_file(self.wd + '/' + args['r'])
			else:
				self.rps_file = ''
		else:
			self.rps_file = ''
		if len(self.blast_files) == 0:
			self.execution=0
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'


	def _check_file(f):
		try:
			open(f)
			return f
		except IOError:
			print('File not found ' + f)
			sys.exit(1)
