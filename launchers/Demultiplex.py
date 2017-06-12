import os.path
import logging
from subprocess import call
import logging as log

class Demultiplex:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = self._create_cmd()


    def check_args (self, args: dict):
        self.wd = os.getcwd()
        self.params = args['params']
        self.sample = args['sample']

        if 'tmp_prefix' in args:
            self.tmp_prefix = args['tmp_prefix']
        else:
            if 'iter' in args:
                if args['iter'] == 'library' or args['iter'] == 'SampleID':
                    self.tmp_prefix = self.sample
                else:
                    self.tmp_prefix = 'prefix'
            else:
                self.tmp_prefix = 'prefix'

        if 'middle' in args:
            self.middle = args['middle']

        if 'min_qual' in args:
            self.min_qual = args['min_qual']

        if 'min_len' in args:
            self.min_len = args['min_len']

        if 'polyA' in args:
            self.polyA = args['polyA']

        if 'mid' in args:
            self.mid = args['mid']
        else:
            self.mid =''

        if 'common' in args:
            self.common = args['common']
        else:
            self.common = ''

        if 'i1' in args:
            self.i1 = self._check_file(self.wd + '/' + args['i1'])
        else:
            log.critical('Need r1 file.')
        if 'i2' in args:
            self.i2 = self._check_file(self.wd + '/' + args['i2'])
        else:
            log.critical('Need r2 file.')
        if 'adapters' in args:
            self.adapters = self._check_file(args['adapters'])
        else:
            self.adapters=''
        if 'o1' in args:
            self.o1 = args['o1']

        if 'o2' in args:
            self.o2 = args['o2']
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
        cmd = 'demultiplex.pl'
        if(self.polyA):
            cmd += ' -polyA'
        if(self.mid != ''):
            for m in self.mid:
                cmd += ' -index ' + m + '=' + self.mid[m]
        if(self.common != ''):
            cmd += ' -c ' + 'COMMON=' + self.common
        cmd += ' -1 ' + self.i1
        cmd += ' -2 ' + self.i2
        cmd += ' -tmp_prefix ' + self.tmp_prefix
        cmd += ' -middle ' + str(self.middle)
        cmd += ' -q ' + str(self.min_qual)
        cmd += ' -l ' + str(self.min_len)
        cmd += ' -a ' + self.adapters
        log.debug(cmd)
        return cmd


    def launch (self):
        if(self.sge):

            qsub_call = "echo '%s' | " % self.cmd
            qsub_call +=   "qsub -wd " + os.getcwd() + " -V -N " + self.sample
            log.debug(qsub_call)
            os.system(qsub_call)
        else:
            log.debug(self.cmd)
            os.system (self.cmd)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
