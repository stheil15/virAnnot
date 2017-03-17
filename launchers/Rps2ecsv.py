import os.path
from subprocess import call
import logging as log
from fabric.api import env
from fabric.operations import run as fabric_run
from fabric.context_managers import settings, hide

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
            self.b = self.wd + '/' + args['b']
        else:
            log.critical('You must provide a blast file.')



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
