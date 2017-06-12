import os.path
from subprocess import call
import logging as log

class Rps2tree:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()


    def create_cmd (self):
        cmd = 'rps2tree.pl'
        if self.iter == 'global':
            for s_id in self.blast_files:
                cmd += ' -id ' + self.blast_files[s_id]['id']
                cmd += ' -s ' + self.blast_files[s_id]['contigs']
                cmd += ' -i ' + self.blast_files[s_id]['pfam']
                cmd += ' -e ' + self.blast_files[s_id]['ecsv']
        else:
            log.debug('msg')
        cmd += ' -mp ' + str(self.min_prot)
        cmd += ' -vp ' + str(self.viral_portion)
        cmd += ' -o ' + self.out
        log.debug(cmd)
        self.cmd.append(cmd)


    def check_args (self, args: dict):
        self.execution=1
        self.wd = os.getcwd()
        self.cmd_file = self.wd + '/' + 'rps2tree_cmd.txt'
        if 'out' in args:
            self.out = args['out']
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'viral_portion' in args:
            self.viral_portion = args['viral_portion']
        if 'min_prot' in args:
            self.min_prot = args['min_prot']
        if 'iter' in args:
            if args['iter'] == 'global':
                self.iter = 'global'
                self.blast_files = {}
                for s_id in args['args']:
                    if s_id not in self.blast_files:
                        if os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['pfam']) and os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['ecsv']) and os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['contigs']):
                            self.blast_files[s_id] = {}
                            self.blast_files[s_id]['pfam'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['pfam']
                            self.blast_files[s_id]['ecsv'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['ecsv']
                            self.blast_files[s_id]['contigs'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['contigs']
                            self.blast_files[s_id]['id'] = args['args'][s_id]['id']
        else:
            log.critical('No iter parameters.')
        if len(self.blast_files.keys()) == 0:
            self.execution=0


    def launch (self):
        if not self.sge:
            for el in self.cmd:
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " 'rps2tree' + ' ' + self.cmd_file
            log.debug(qsub_call)
            os.system(qsub_call)


    def _check_file (self,f):
        try:
            open(f)
            return f
        except IOError:
            print('File not found ' + f)
