"""
Authors: Marie Lefebvre
"""
# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys


class Ecsv2merge:
	"""
	This class is a part of the virAnnot module.
	Merge blastx results for singletons ans contigs. CSV format.
	"""

	def __init__(self, args):
		self.execution = 1
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd(self):
		"""
		Create command using rps2merge submodule (tool)
		"""
		cmd = 'ecsv2merge.py'
		cmd += ' -c ' + self.contigs
		cmd += ' -s ' + self.singletons
		cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)
		


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
		self.wd = os.getcwd()
		self.params = args['params']
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.cmd_file = self.wd + '/' + self.sample + '/' + self.sample + '_ecsv2merge_cmd.txt'
		if 'out' in args:
			self.out = self.wd + '/' + self.sample + '/' + args['out']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		# contigs files
		if 'contigs' in args:
			self.contigs = self.wd + '/' + self.sample + '/' + args['contigs']
		# singletons files
		if 'singletons' in args:
			self.singletons = self.wd + '/' + self.sample + '/' + args['singletons']

