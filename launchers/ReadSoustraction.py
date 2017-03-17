import os.path
import logging
from subprocess import call

class ReadSoustraction:

    def __init__ (self,i1,i2,db,o1,o2,n_cpu,sge,params,sample):
        self.i1 = self._check_file(i1)
        self.i2 = self._check_file(i2)
        self.db = self._check_bowtie_db(db)
        self.o1 = o1
        self.o2 = o2
        self.n_cpu = n_cpu
        self.sge = sge
        self.params = params
        self.sample = sample
        self.cmd = self._create_cmd()


    def launch (self):
        if(self.sge):
            qsub_call = "echo '%s' | " % self.cmd
            qsub_call +=   "qsub -wd " + os.getcwd() + " -V -N " + self.sample + '_rs_' + self.db
            os.system(qsub_call)
        else:
            os.system (self.cmd)


    def _create_cmd (self):
        cmd = ''
        cmd += self.params['bin']['bowtie'] + ' -p ' + str(self.n_cpu)
        cmd += ' -x ' + self.params['ReadSoustraction']['db']['vitis']
        cmd += ' -1 ' + self.i1 + ' -2 ' + self.i2 + ' | '
        cmd += self.params['bin']['samtools'] + ' view -bS - '
        cmd += ' | ' + self.params['bin']['samtools'] + ' view -u -f 12 -F 256 - | ' + self.params['bin']['bedtools'] + ' bamtofastq -i - -fq ' + self.o1 + ' -fq2 ' + self.o2
        return cmd


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)


    def _check_bowtie_db (self,db):
        return db
