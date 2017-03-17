import os.path
from subprocess import call
import logging as log

class Getresults:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()

    def create_cmd (self):
        cmd = 'mkdir ' + self.out + "\n"
        cmd += 'cp -r ' + self.rps2tree + ' ' + self.out + "\n"
        cmd += 'cp -r ' + self.krona + ' ' + self.out + "\n"
        cmd += 'cp ' + self.krona + ' ' + self.out + "\n"
        for s_id in self.files:
            cmd += 'mkdir ' + self.out + '/' + s_id + "\n"
            cmd += 'cp ' + self.files[s_id]['contigs'] + ' ' + self.out + '/' + s_id + "\n"
            cmd += 'cp ' + self.files[s_id]['i1'] + ' ' + self.out + '/' + s_id + "\n"
            cmd += 'cp ' + self.files[s_id]['i2'] + ' ' + self.out + '/' + s_id + "\n"
            cmd += 'mkdir ' + self.out + '/' + s_id + '/autoMapper' + "\n"
            cmd += 'cp -r ' + self.files[s_id]['autoMapper'] + '/results/*' + ' ' + self.out + '/' + s_id + '/autoMapper' + "\n"
            cmd += 'cp ' + self.files[s_id]['excel'] + ' ' + self.out + '/' + s_id + "\n"
            cmd += 'cp -r ' + self.files[s_id]['drVM'] + ' ' + self.out + '/' + s_id + "\n"
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_args (self, args: dict):
        self.wd = os.getcwd()
        self.cmd_file = self.wd + '/' + 'getResults_cmd.txt'
        if 'out' in args:
            self.out = self.wd + '/' + args['out']
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'rps2tree' in args:
            self.rps2tree = args['rps2tree']
        if 'krona' in args:
            self.krona = args['krona']
        if 'iter' in args:
            if args['iter'] == 'global':
                self.iter = 'global'
                self.files = {}
                for s_id in args['args']:
                    if s_id not in self.files:
                        self.files[s_id] = {}
                    for p in args['args'][s_id]:
                        if p == 'contigs':
                            self.files[s_id]['contigs'] = self._check_file(self.wd + '/' + s_id + '/' + args['args'][s_id]['contigs'])
                        if p == 'i1':
                            self.files[s_id]['i1'] = self._check_file(self.wd + '/' + s_id + '/' + args['args'][s_id]['i1'])
                        if p == 'i2':
                            self.files[s_id]['i2'] = self._check_file(self.wd + '/' + s_id + '/' + args['args'][s_id]['i2'])
                        if p == 'autoMapper':
                            self.files[s_id]['autoMapper'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['autoMapper']
                        if p == 'drVM':
                            self.files[s_id]['drVM'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['drVM']
                        if p == 'excel':
                            self.files[s_id]['excel'] = self._check_file(self.wd + '/' + s_id + '/' + args['args'][s_id]['excel'])


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
