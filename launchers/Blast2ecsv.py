import os.path
from subprocess import call
import logging as log
from fabric.api import env
from fabric.operations import run as fabric_run
from fabric.context_managers import settings, hide

class Blast2ecsv:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()

    def create_cmd (self):
        cmd = 'blast2ecsv.pl'
        cmd += ' -t ' + self.type
        cmd += ' -b ' + self.b
        cmd += ' -if ' + self.in_format
        cmd += ' -e ' + str(self.evalue)
        cmd += ' -pm ' + self.pm
        if self.fhit:
            cmd += ' -fhit '
        cmd += ' -seq ' + self.contigs
        if self.vs:
            cmd += ' -vs '
        if self.r:
            cmd += ' -r '
        cmd += ' -o ' + self.out
        if self.rn != '':
            cmd += ' -rn ' + self.rn
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
        if 'fhit' in args:
            self.fhit = bool(args['fhit'])
        else:
            self.fhit = False
        if 'pm' in args:
            self.pm = args['pm']
        else:
            self.pm = 'local'
        if 'if' in args:
            self.in_format = args['if']
        else:
            self.in_format = 'xml'
        if 'r' in args:
            self.r = bool(args['r'])
        else:
            self.r = False
        if 'vs' in args:
            self.vs = bool(args['vs'])
        else:
            self.vs = False
        if 'b' in args:
            self.b = self.wd + '/' + args['b']
        else:
            log.critical('You must provide a blast file.')
        if 'rn' in args:
            self.rn = self.wd + '/' + args['rn']
        else:
            self.rn = ''
        if 'type' in args:
            self.type = args['type']


    def launch (self):
        if not self.sge:
            for el in self.cmd:
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + self.sample + '_blast2ecsv' + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
