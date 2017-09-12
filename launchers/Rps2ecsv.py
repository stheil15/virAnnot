import os.path
from subprocess import call
import logging as log

class Rps2ecsv:

	def __init__ (self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd (self):
		cmd = 'rps2ecsv.pl'
		cmd += ' -b ' + self.b
		cmd += ' -e ' + str(self.evalue)
		cmd += ' -o ' + self.out
		self.cmd.append(cmd)


	def check_args (self, args: dict):
		self.execution=1
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_blast2ecsv_cmd.txt'
		if 'contigs' in args:
			self.contigs = self.wd + '/' + args['contigs']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'out' in args:
			self.out = args['out']
		if 'params' in args:
			self.params = args['params']
		if 'evalue' in args:
			self.evalue = args['evalue']
		else:
			self.evalue = 0.0001
		if 'b' in args:
			if os.path.exists(self.wd + '/' + args['b']):
				self.b = self.wd + '/' + args['b']
			else:
				self.b = ''
				self.execution=0
		else:
			log.critical('You must provide a blast file.')
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'

	def _check_file (self,f):
		try:
			open(f)
			return f
		except IOError:
			print('File not found ' + f)
