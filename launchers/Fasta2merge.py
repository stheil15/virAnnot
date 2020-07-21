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


class Fasta2merge:
	"""
	This class is a part of the virAnnot module.
	It merges fasta signletons + contings
	"""

	def __init__(self, args):
		self.execution = 1
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd(self):
		"""
		Create command to merge two fasta files
		"""
		if os.path.exists(self.contigs):
			cmd = 'cat ' + self.contigs + ' > ' + self.out
			log.debug(cmd)
			self.cmd.append(cmd)
		if os.path.exists(self.singletons):
			cmd = 'cat ' + self.singletons + ' >> ' + self.out
			log.debug(cmd)
			self.cmd.append(cmd)
		if os.path.exists(self.contigs) != True & os.path.exists(self.singletons) != True:
			cmd = 'echo "" > ' + self.out
			log.debug(cmd)
			self.cmd.append(cmd)
			log.debug('No files available')


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
		self.wd = os.getcwd()
		self.params = args['params']
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.cmd_file = self.wd + '/' + self.sample + '/' + self.sample + '_fasta2merge_cmd.txt'
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

