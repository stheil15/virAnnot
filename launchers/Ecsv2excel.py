import os.path
from subprocess import call
import logging as log
from fabric.api import env
from fabric.operations import run as fabric_run
from fabric.context_managers import settings, hide
import sys

class Ecsv2excel:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()


    def create_cmd (self):
        cmd = 'ecsv2excel.pl'
        for c in self.blast_files:
            cmd += ' -b ' + str(c)
        cmd += ' -r ' + self.r
        cmd += ' -o ' + self.out
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_args (self, args: dict):
        self.sample = args['sample']
        self.wd = os.getcwd() + '/' + self.sample
        self.cmd_file = self.wd + '/' + self.sample + '/' + 'ecsv2excel_cmd.txt'
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
                self.blast_files.append(self._check_file(self.wd + '/' + args[opt_name]))
        if 'r' in args:
            self.r = self._check_file(self.wd + '/' + args['r'])



    def launch (self):
        if not self.sge:
            for el in self.cmd:
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + self.sample + '_rps2ecsv' + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
            sys.exit(1)
