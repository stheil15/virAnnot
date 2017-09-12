import os.path
from subprocess import call
import logging as log


class Ecsv2krona:

	def __init__ (self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd (self):
		cmd = 'ecsv2krona.pl'
		cmd += ' -data ' + self.data
		if self.r:
			cmd += ' -r '
		cmd += ' -c ' + self.c
		for s_id in self.blast_files:
			for i in range(0,len(self.blast_files[s_id]['id'])):
				cmd += ' -id ' + self.blast_files[s_id]['id'][i]
				cmd += ' -i ' + self.blast_files[s_id]['csv_file'][i]
				if len(self.blast_files[s_id]['xml_file']) > 0:
					cmd += ' -x ' + self.blast_files[s_id]['xml_file'][i]
		cmd += ' -mt 1 -m'
		if self.out != '':
			cmd += ' -o ' + self.out
		if self.outdir != '':
			cmd += ' -outdir ' + self.outdir
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args (self, args: dict):
		self.execution=1
		if 'out' in args:
			self.out = args['out']
		else:
			self.out = ''
		if 'outdir' in args:
			self.outdir = args['outdir']
		else:
			self.outdir = ''
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'data' in args:
			self.data = args['data']
		if 'r' in args:
			self.r = bool(args['r'])
		if 'c' in args:
			self.c = args['c']
		if 'iter' in args:
			if args['iter'] == 'global':
				self.wd = os.getcwd()
				self.cmd_file = self.wd + '/' + 'ecsv2krona_cmd.txt'
				if 'out' in args:
					self.out = args['out']
				self.iter = 'global'
				self.blast_files = {}
				for s_id in args['args']:
					if s_id not in self.blast_files:
						self.blast_files[s_id] = {}
						self.blast_files[s_id]['csv_file'] = []
						self.blast_files[s_id]['xml_file'] = []
						self.blast_files[s_id]['id'] = []
					for i in range(1, 100, 1):
						id_name = 'id' + str(object=i)
						opt_name = 'b' + str(object=i)
						xml_opt_name = 'x' + str(object=i)
						if id_name not in args['args'][s_id] and opt_name not in args['args'][s_id]:
							continue
						if opt_name in args['args'][s_id]:
							if os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id][opt_name]):
								self.blast_files[s_id]['csv_file'].append(self.wd + '/' + s_id + '/' + args['args'][s_id][opt_name])
								self.blast_files[s_id]['id'].append(args['args'][s_id][id_name])
								if xml_opt_name in args['args'][s_id]:
									self.blast_files[s_id]['xml_file'].append(self.wd + '/' + s_id + '/' + args['args'][s_id][xml_opt_name])
			elif args['iter'] == 'sample':
				self.sample = args['sample']
				self.wd = os.getcwd() + '/' + self.sample
				self.cmd_file = self.wd + '/' + 'ecsv2krona_cmd.txt'
				if 'out' in args:
					self.out = args['out']
				self.iter = 'sample'
				self.blast_files = {}
				if self.sample not in self.blast_files:
					self.blast_files[self.sample]={}
					self.blast_files[self.sample]['csv_file'] = []
					self.blast_files[self.sample]['xml_file'] = []
					self.blast_files[self.sample]['id'] = []
				for i in range(1, 100, 1):
					id_name = 'id' + str(object=i)
					opt_name = 'b' + str(object=i)
					xml_opt_name = 'x' + str(object=i)
					if id_name not in args and opt_name not in args:
						continue
					if os.path.exists(self.wd + '/' + args[opt_name]):
						if opt_name in args:
							self.blast_files[self.sample]['csv_file'].append(self._check_file(self.wd + '/' + args[opt_name]))
							self.blast_files[self.sample]['id'].append(args[id_name])
						if xml_opt_name in args:
							self.blast_files[self.sample]['xml_file'].append(self._check_file(self.wd + '/' + args[xml_opt_name]))
				if len(self.blast_files[self.sample]['csv_file']) == 0:
					self.execution=0
			# elif args['iter'] == 'library':
			# 	if 'out' in args:
			# 		self.out = args['out']
			# 	self.iter = 'library'
			# 	self.library = args['library']
			# 	print(self.library)
			# 	self.blast_files = {}
			# 	if self.library not in self.blast_files:
			# 		self.blast_files[self.library]={}
			# 		self.blast_files[self.library]['csv_file'] = []
			# 		self.blast_files[self.library]['xml_file'] = []
			# 		self.blast_files[self.library]['id'] = []
			# 	for i in range(1, 100, 1):
			# 		id_name = 'id' + str(object=i)
			# 		opt_name = 'b' + str(object=i)
			# 		xml_opt_name = 'x' + str(object=i)
			# 		if id_name not in args and opt_name not in args:
			# 			continue
			# 		print(args)
			# 		if os.path.exists(self.wd + '/' + args[self.library][opt_name]):
			# 			if opt_name in args:
			# 				self.blast_files[self.library]['csv_file'].append(self._check_file(self.wd + '/' + args[opt_name]))
			# 				self.blast_files[self.library]['id'].append(args[id_name])
			# 			if xml_opt_name in args:
			# 				self.blast_files[self.library]['xml_file'].append(self._check_file(self.wd + '/' + args[xml_opt_name]))
			else:
				log.error('iter argument not handled. ' + args['iter'])
				self.execution=0
		if len(self.blast_files.keys()) == 0:
			self.execution=0



	def _check_file (self,f):
		try:
			open(f)
			return f
		except IOError:
			print('File not found ' + f)
