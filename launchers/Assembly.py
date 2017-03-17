import os.path
from subprocess import call
import logging as log
import shutil

class Assembly:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()


    def create_cmd (self):
        if self.prog == 'idba':
            self._create_idba_cmd()


    def _create_idba_cmd (self):
        self.i1 = self.check_seq_format(self.i1)
        self.i2 = self.check_seq_format(self.i2)
        merged_read = self.i1.split('_')[0] + '_merged.fa'
        cmd = 'fq2fa --merge ' + self.i1 + ' ' + self.i2 + ' ' + merged_read
        log.debug(cmd)
        self.cmd.append(cmd)
        cmd = 'idba_ud -r ' + merged_read + ' --num_threads ' + self.n_cpu + ' -o ' + self.wd + '/' + self.sample + '_idba'
        log.debug(cmd)
        self.cmd.append(cmd)
        cmd = 'sed -i \'s,^>scaffold_\([0-9]*\),>' + self.sample + '_\\1,\' ' + self.wd + '/' + self.sample + '_idba' + '/scaffold.fa'
        log.debug(cmd)
        self.cmd.append(cmd)
        cmd = 'cp ' + self.wd + '/' + self.sample + '_idba' + '/scaffold.fa' + ' ' + self.wd + '/' + self.out
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
            log.debug('Format seems to be fastq.')
            return in_file


    def check_args (self, args: dict):
        if 'sample' in args:
            self.sample = str(args['sample'])
        self.wd = os.getcwd() + '/' + self.sample
        self.cmd_file = self.wd + '/' + self.sample + '_idba_cmd.txt'
        if 'i1' in args:
            self.i1 = self._check_file(self.wd + '/' + args['i1'])
        else:
            log.critical('Need r1 file.')
        if 'i2' in args:
            self.i2 = self._check_file(self.wd + '/' + args['i2'])
        else:
            log.critical('Need r2 file.')
        if 'prog' in args:
            if args['prog'] == 'idba':
                self.prog = args['prog']
            else:
                log.critical('Wrong assembly program name.')
        else:
            log.critical('Program name is mandatory.')
        if 'n_cpu' in args:
            self.n_cpu = str(args['n_cpu'])
        else:
            log.debug('n_cpu option not found. default 1')
            self.n_cpu = '1'
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'out' in args:
            self.out = args['out']


    def launch (self):
        if(self.sge):
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + self.sample + '_drVM' + ' ' + self.cmd_file
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
