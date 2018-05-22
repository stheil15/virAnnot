import os.path
from subprocess import call
import logging as log

class Getresults:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()

    def create_cmd (self):
        cmd = ''
        if not os.path.exists(self.out):
            cmd = 'mkdir ' + self.out + "\n"
        for d in self.directories:
            if os.path.exists(d):
                cmd += 'cp -r ' + d + ' ' + self.out + "\n"
        for f in self.files:
            if os.path.exists(f):
                cmd += 'cp ' + f + ' ' + self.out + "\n"
        for s_id in self.sample_files:
            if not os.path.exists(self.out + '/' + s_id):
                cmd += 'mkdir ' + self.out + '/' + s_id + "\n"
            for f in self.sample_files[s_id]:
                if os.path.exists(s_id + '/' + f):
                    cmd += 'cp ' + s_id + '/' + f + ' ' + self.out + '/' + s_id + "\n"
        for s_id in self.sample_dir:
            if not os.path.exists(self.out + '/' + s_id):
                cmd += 'mkdir ' + self.out + '/' + s_id + "\n"
            for d in self.sample_dir[s_id]:
                if os.path.exists(s_id + '/' + d):
                    cmd += 'cp -r ' + s_id + '/' + d + ' ' + self.out + '/' + s_id + "\n"
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_args (self, args: dict):
        self.execution=1
        self.wd = os.getcwd()
        self.cmd_file = self.wd + '/' + 'getResults_cmd.txt'
        if 'out' in args:
            self.out = self.wd + '/' + args['out']
        self.files = []
        self.directories = []
        self.sample_dir = {}
        self.sample_files = {}
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'global_dir' in args:
            self.directories = args['global_dir']
        if 'global_file' in args:
            self.files = args['global_file']
        if 'sample_dir' in args:
            for s in args['sample_dir']:
                if s not in self.sample_dir:
                    self.sample_dir[s]=[]
                self.sample_dir[s] = args['sample_dir'][s]
        if 'sample_files' in args:
            for s in args['sample_files']:
                if s not in self.sample_files:
                    self.sample_files[s]=[]
                self.sample_files[s] = args['sample_files'][s]

    def launch (self):
        if not self.sge:
            for el in self.cmd:
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + 'Getresults' + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
