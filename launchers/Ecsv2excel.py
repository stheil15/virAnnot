#!/user/bin/pyton3.4
"""
Authors         :Sebastien Theil, Marie Lefebvre
usage           :python3.4 Ecsv2exel.py
python_version  :3.4
"""
import os.path
import logging as log

class Ecsv2excel:
	"""
	This class is part of virAnnot module
	It creates the command that will convert
	the dispatched csv files to a summary Excel file
	"""

	def __init__(self, args):
		self.execution = 0
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd(self):
		"""
		Create command
		"""
		cmd = 'ecsv2excel.pl'
		for bfile in self.blast_files:
			cmd += ' -b ' + str(bfile)
		if self.r != '':
			cmd += ' -r ' + self.r
		cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
		self.sample = args['sample']
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + 'ecsv2excel_cmd.txt'
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		self.blast_files = []
		for i in range(1, 10, 1):
			opt_name = 'b' + str(object=i)
			if opt_name in args:
				if os.path.exists(self.wd + '/' + args[opt_name]):
					self.blast_files.append(self.wd + '/' + args[opt_name])
		if 'r' in args:
			if os.path.exists(self.wd + '/' + args['r']):
				self.r = self._check_file(self.wd + '/' + args['r'])
			else:
				self.r = ''
		else:
			self.r = ''
		if len(self.blast_files) == 0:
			self.execution = 0
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'


	def _check_file(self, input_file):
		"""
		Verify that file exists
		"""
		try:
			open(input_file)
			return input_file
		except IOError:
			print 'File not found ' + input_file
			self.execution = 0
