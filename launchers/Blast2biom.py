"""
	This class is a part of the virAnnot module
	============
	Authors: Sebastien Theil, Marie Lefebvre
"""

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log

class Blast2biom:
	"""
	This module is part of virAnnot module
	It creates the command that will generate biom (with OTUs)
	from BLast results
	"""


	def __init__(self, args):
		self.check_args(args)
		self.create_cmd()

	def create_cmd(self):
		"""
		Create command
		"""
		cmd = 'blast2biom.py --out ' + self.out
		if self.iter == 'global':
			for s_id in self.blast_files:
				for i in range(0, len(self.blast_files[s_id]['id'])):
					cmd += ' -i ' + self.blast_files[s_id]['id'][i]
					cmd += ' -b ' + self.blast_files[s_id]['csv_file'][i]
		cmd += ' --levels 5 --regex Viruses --regex Viroids'
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
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
		self.wd = os.getcwd()
		self.cmd_file = self.wd + '/' + 'blast2biom_cmd.txt'
		self.cmd = []
		if 'iter' in args:
			if args['iter'] == 'global':
				self.iter = 'global'
				self.blast_files = {}
				for s_id in args['args']:
					for i in range(1, 100, 1):
						id_name = 'id' + str(object=i)
						opt_name = 'b' + str(object=i)
						if id_name not in args['args'][s_id] and opt_name not in args['args'][s_id]:
							continue
						if opt_name in args['args'][s_id]:
							if os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id][opt_name]):
								if s_id not in self.blast_files:
									self.blast_files[s_id] = {}
									self.blast_files[s_id]['csv_file'] = []
									self.blast_files[s_id]['id'] = []
								self.blast_files[s_id]['csv_file'].append(self.wd + '/' + s_id + '/' + args['args'][s_id][opt_name])
								self.blast_files[s_id]['id'].append(args['args'][s_id][id_name])
		if len(self.blast_files.keys()) == 0:
			self.execution = 0
