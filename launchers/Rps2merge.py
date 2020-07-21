# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log


class Rps2merge:
	"""
	This class is a part of the virAnnot module.
	It merges OTU assignation results with blastx and rps results
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
		cmd = 'rps2merge.py'
		for s_id in self.blast_files:
			cmd += ' -b ' + self.blast_files[s_id]['blastx']
			cmd += ' -p ' + self.blast_files[s_id]['pfam']
		cmd += ' -r ' + self.rps_path
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
		self.cmd_file = self.wd + '/' + 'rps2merge_cmd.txt'
		if 'out' in args:
			self.out = args['out']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		# blastx results files
		if 'iter' in args:
			if args['iter'] == 'global':
				self.iter = 'global'
				self.blast_files = {}
				for s_id in args['args']:
					if s_id not in self.blast_files:
						if os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['pfam']) and os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['blastx']):
							self.blast_files[s_id] = {}
							self.blast_files[s_id]['pfam'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['pfam']
							self.blast_files[s_id]['blastx'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['blastx']
							self.blast_files[s_id]['id'] = args['args'][s_id]['id']
						else:
							self.execution = 0
		else:
			log.critical('No iter parameters.')
		if 'rps_folder' in args:
			if os.path.exists(self.wd + '/' + args['rps_folder']):
				self.rps_path = self.wd + '/' + args['rps_folder']
			else:
				log.critical('Rps2tree folder not found.')
				self.execution = 0
		else:
			self.rps_path = self.wd + '/Rps2tree_global'

