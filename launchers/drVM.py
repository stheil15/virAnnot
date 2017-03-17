import os.path
from subprocess import call
import logging as log

class drVM:

    def __init__ (self, args):
        log.info( 'drVM module')
        self.cmd = []
        self.check_args(args)
        self.cmd_file = self.sample + '_drVM_cmd.txt'
        self._create_cmd()

    def _create_cmd (self):
        cmd = 'drVM.py'
        cmd += ' -1 ' + str(self.i1)
        cmd += ' -2 ' + str(self.i2)
        cmd += ' -t ' + str(self.n_cpu)
        cmd += ' -bi ' + str(self.identity)
        cmd += ' -cl ' + str(self.min_len)
        cmd += ' -o ' + self.sample + '/' + self.sample + '_drVM'
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_seq_format (self, in_file):
        in_file = str(object=in_file)
        self._check_file(in_file)
        out = ''
        if in_file.lower().endswith('.fa') or in_file.lower().endswith('.fasta') or in_file.lower().endswith('.fas'):
            out = os.path.splitext(in_file)[0] + '.fq'
            if not self._check_file(out):
                cmd = 'fasta_to_fastq' + ' ' + in_file + ' > ' + out
                log.debug(str(cmd))
                self.cmd.append(cmd)
            return out
        else:
            log.debug('Format seems to be fasta.')
            return in_file



    def check_args (self, args: dict):
        if 'i1' in args:
            self.i1 = self.check_seq_format(args['i1'])
        if 'i2' in args:
            self.i2 = self.check_seq_format(args['i2'])
        if 'n_cpu' in args:
            self.n_cpu = str(args['n_cpu'])
        else:
            log.debug('n_cpu option not found. default 1')
            self.n_cpu = '1'
        if 'identity' in args:
            self.identity = str(args['identity'])
        else:
            log.debug('identity option not found. default 80')
            self.identity = '80'
        if 'min_len' in args:
            self.min_len = str(args['min_len'])
        else:
            log.debug('min_len option not found. default 1000')
            self.min_len = '1000'
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'sample' in args:
            self.sample = str(args['sample'])

    def launch (self):
        if(self.sge):
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + os.getcwd() + " -V -N " + self.sample + '_drVM' + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)
        else:
            for el in self.cmd:
                os.system (el)

    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
