import os.path
from subprocess import call
import logging as log

class Normalization:

	def __init__ (self, args):
		self.check_args(args)
		self._create_cmd()


	def check_args (self, args: dict):
		self.wd = os.getcwd()
		self.params = args['params']
		self.cmd = []
		self.execution=1
		if 'iter' in args:
			if args['iter'] == 'sample':
				self.sample = args['sample']
				self.cmd_file = self.wd + '/' + 'Normalization_' + self.sample + '_cmd.txt'
				self.iter = 'sample'
			elif args['iter'] == 'library':
				self.library = args['library']
				self.cmd_file = self.wd + '/'  + 'Normalization_' + self.library + '_cmd.txt'
				self.iter = 'library'
		if 'i1' in args:
			self.i1 = self._check_file(self.wd + '/' + args['i1'])
		else:
			log.critical('Need r1 file.')
			self.execution=0
		if 'i2' in args:
			self.i2 = self._check_file(self.wd + '/' + args['i2'])
		else:
			log.critical('Need r2 file.')
			self.execution=0
		if 'o1' in args:
			self.o1 = args['o1']
		else:
			log.critical('Need o1 file.')
		if 'o2' in args:
			self.o2 = args['o2']
		else:
			log.critical('Need o2 file.')
		if 'num' in args:
			self.num = args['num']
		else:
			log.critical('Need number of subsample reads.')
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			log.debug('n_cpu option not found. default 1')
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False


	def _create_cmd (self):
		cmd = ''
		cmd += self.params['bin']['seqtk'] + ' sample -s1000 ' + self.i1 + ' ' + str(self.num)
		cmd += ' > ' + self.o1 
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = ''
		cmd += self.params['bin']['seqtk'] + ' sample -s1000 ' + self.i2 + ' ' + str(self.num)
		cmd += ' > ' + self.o2 
		log.debug(cmd)
		self.cmd.append(cmd)


	def _check_file (self,f):
		try:
			open(f)
			return f
		except IOError:
			print('File not found ' + f)
