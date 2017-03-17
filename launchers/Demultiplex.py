import os.path
import logging
from subprocess import call

class Demultiplex:

    def __init__ (self, i1: str, i2: str, adapters: str, tmp_prefix: str, middle:int, min_qual:int, polyA:bool, min_len:int, sge: bool, params: dict, sample:str, mid: None, common: None):
        self.i1 = self._check_file(i1)
        self.i2 = self._check_file(i2)
        self.adapters = self._check_file(adapters)
        self.tmp_prefix = tmp_prefix
        self.middle = middle
        self.min_qual = min_qual
        self.min_len = min_len
        self.polyA = polyA
        self.min_len = min_len
        self.sge = sge
        self.sample = sample
        if(mid is not None):
            self.mid = mid
        if(common is not None):
            self.common = common
        self.cmd = self._create_cmd()


    def _create_cmd (self):
        cmd = 'demultiplex_v2.pl'
        if(self.polyA):
            cmd += ' -polyA'
        if(self.mid):
            for m in self.mid:
                cmd += ' -index ' + m + '=' + self.mid[m]
        if(self.common):
            cmd += ' -c ' + 'COMMON=' + self.common
        cmd += ' -1 ' + self.i1
        cmd += ' -2 ' + self.i2
        cmd += ' -tmp_prefix ' + self.tmp_prefix
        cmd += ' -middle ' + str(self.middle)
        cmd += ' -q ' + str(self.min_qual)
        cmd += ' -l ' + str(self.min_len)
        cmd += ' -a ' + self.adapters
        return cmd


    def launch (self):
        if(self.sge):
            qsub_call = "echo '%s' | " % self.cmd
            qsub_call +=   "qsub -wd " + os.getcwd() + " -V -N " + self.sample
            os.system(qsub_call)
        else:
            os.system (self.cmd)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
