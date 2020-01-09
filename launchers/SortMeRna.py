import os.path
from subprocess import call
import logging as log

class SortMeRna:

    def __init__ (self, args):
        self.check_args(args)
        if self.execution==1:
            self.cmd = self._create_cmd()
        else:
            self.cmd = ''


    def check_args (self, args=dict):
        self.wd = os.getcwd()
        self.params = args['params']
        self.sample = args['sample']
        self.execution=1
        if 'i1' in args and os.path.exists(self.wd + '/' + args['i1']):
            self.i1 = self.wd + '/' + args['i1']
        else:
            self.i1 = ''
            log.critical('Need r1 file.')
            self.execution=0

        if 'i2' in args and os.path.exists(self.wd + '/' + args['i2']):
            self.i2 = self.wd + '/' + args['i2']
        else:
            self.i2 = ''
            log.critical('Need r2 file.')
            self.execution=0

        if 'db' in args:
            self.db = []
            self.db = args['db'].split(',')
        else:
            log.critical('Need Sortmerna db option.')
            self.execution=0
        if 'o1' in args:
            self.o1 = args['o1']
        else:
            self.o1 = self.wd + '/' + self.sample + '/' + 'sortmerna.r1'
        if 'o2' in args:
            self.o2 = args['o2']
        else:
            self.o2 = self.wd + '/' + self.sample + '/' + 'sortmerna.r2'

        if 'n_cpu' in args:
            self.n_cpu = str(args['n_cpu'])
        else:
            log.debug('n_cpu option not found. default 1')
            self.n_cpu = '1'

        if 'evalue' in args:
            self.evalue = args['evalue']
        else:
            self.evalue = 1

        if 'aligned' in args:
            self.aligned = args['aligned']
        else:
            self.aligned = self.wd + '/' + self.sample + '/' + self.sample + '_RNA'

        if 'm' in args:
            self.m = args['m']
        else:
            self.m = 1024

        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False


    def launch (self):
        if(self.sge):
            qsub_call = "echo '%s' | " % self.cmd
            qsub_call +=   "qsub -wd " + os.getcwd() + " -V -N " + self.sample + '_sortRna' + ' -pe multithread ' + self.n_cpu
            log.info(qsub_call)
            os.system(qsub_call)
        else:
            log.info('Launching local command:')
            log.info(self.cmd)
            os.system (self.cmd)


    def _create_cmd (self):
        cmd = ''
        cmd += self.params['bin']['merge-paired-reads'] + ' ' + self.i1 + ' ' + self.i2 + ' ' + self.wd + '/' + self.sample + '/merged-paired-reads.fq' + "\n"
        cmd += self.params['bin']['sortmerna'] + ' -a ' + str(self.n_cpu)
        cmd += ' --ref '
        for i in range(0,len(self.db)):
            cmd += self.params['SortMeRna']['db'][self.db[i]] + '.fasta,' + self.params['SortMeRna']['db'][self.db[i]] + '_accession_taxonomy.txt'
            if i != len(self.db)-1:
                cmd += ':'
        cmd += ' --reads ' + self.wd + '/' + self.sample + '/merged-paired-reads.fq'
        cmd += ' --num_alignments 1'
        cmd += ' --aligned ' + self.aligned
        cmd += ' --paired_out '
        cmd += ' --other ' + self.wd + '/' + self.sample + '/' + self.sample + '_nonRNA' + ' --fastx --log ' + "\n"
        cmd += self.params['bin']['unmerge-paired-reads']+ ' ' + self.wd + '/' + self.sample + '/' + self.sample + '_nonRNA.fq' + ' ' + self.o1 + ' ' + self.o2 + "\n"
        log.debug(cmd)
        return cmd
