import os.path
from subprocess import call
import logging as log

class BlastCompare:

	def __init__ (self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd (self):
		cmd = 'compare_annot.py'
		for i in range(0,len(self.blast_files['id'])):
			cmd += ' --id ' + self.blast_files['id'][i]
			cmd += ' --blast ' + self.blast_files['csv_file'][i]
		cmd += ' --out ' + self.out
		self.cmd.append(cmd)
		log.debug(cmd)


	def check_args (self, args=dict):
		self.execution=1
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_compare_cmd.txt'
		self.blast_files = {}
		self.blast_files['csv_file'] = []
		self.blast_files['id'] = []
		for i in range(1, 100, 1):
			id_name = 'id' + str(object=i)
			opt_name = 'b' + str(object=i)
			if id_name not in args and opt_name not in args:
				continue
			else:
				if os.path.exists(self.wd + '/' + args[opt_name]):
					self.blast_files['csv_file'].append(self.wd + '/' + args[opt_name])
					self.blast_files['id'].append(args[id_name])
		if len(self.blast_files['csv_file']) == 0:
			self.execution=0
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
