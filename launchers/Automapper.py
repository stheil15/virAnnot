import os.path
from subprocess import call
import logging as log
from fabric.api import env
from fabric.operations import run as fabric_run
from fabric.context_managers import settings, hide

class Automapper:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()


    def create_cmd (self):
        cmd = 'autoMapper_v2.pl'
        cmd += ' -q ' + self.contigs
        cmd += ' -b ' + self.ecsv
        if self.reads != '':
            if not os.path.exists(self.reads):
                cmd = 'fq2fa --merge ' + self.i1 + ' ' + self.i2 + ' ' + self.reads
                log.debug(cmd)
                self.cmd.append(cmd)
            cmd += ' -reads ' + self.reads
        cmd += ' -o ' + self.out
        cmd += ' -p '
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_args (self, args: dict):
        if 'sample' in args:
            self.sample = args['sample']
        self.wd = os.getcwd() + '/' + self.sample
        self.cmd_file = self.wd + '/' + 'autoMapper_cmd.txt'
        if 'out' in args:
            self.out = self.wd + '/' + args['out']
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'contigs' in args:
            self.contigs = self.wd + '/' + args['contigs']
        if 'ecsv' in args:
            self.ecsv = self.wd + '/' + args['ecsv']
        if 'i1' in args:
            self.i1 = self.check_seq_format(self.wd + '/' + args['i1'])
        else:
            self.i1 = ''
            log.debug('No r1 file.')
        if 'i2' in args:
            self.i2 = self.check_seq_format(self.wd + '/' + args['i2'])
        else:
            self.i2 = ''
            log.debug('No r2 file.')
        if self.i1 != '' and self.i2 != '':
            self.reads = self.i1.split('_')[0] + '_merged.fa'
        else:
            self.reads = ''


    def launch (self):
        if not self.sge:
            for el in self.cmd:
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + self.sample + '_rps2tree' + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)


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
