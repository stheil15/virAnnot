import os.path
from subprocess import call
import logging as log

class Blast2hist:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()

    def create_cmd (self):
        cmd = 'blast2hist.py --out ' + self.out
        if self.iter == 'global':
            for s_id in self.blast_files:
                for i in range(0,len(self.blast_files[s_id]['id'])):
                    cmd += ' -i ' + self.blast_files[s_id]['id'][i]
                    cmd += ' -b ' + self.blast_files[s_id]['csv_file'][i]
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_args (self, args: dict):
        self.execution=1
        if 'out' in args:
            self.out = args['out']
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        self.wd = os.getcwd()
        self.cmd_file = self.wd + '/' + 'blast2hist_cmd.txt'
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
            self.execution=0

    def launch (self):
        if not self.sge:
            for el in self.cmd:
                log.debug(el)
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + __name__ + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
