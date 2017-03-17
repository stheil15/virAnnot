import os.path
from subprocess import call
import logging as log
from fabric.api import env
from fabric.operations import run as fabric_run
from fabric.context_managers import settings, hide

class Blast:

    def __init__ (self, args):
        self.check_args(args)
        self.cmd = []
        self.ssh_cmd = []
        self.create_cmd()

    def create_cmd (self):
        cmd=''
        cmd = 'scp ' + self.contigs + ' ' + self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch']
        log.debug(cmd)
        self.cmd.append(cmd)
        ssh_cmd = self._get_exec_script()
        fw =  open(self.remote_cmd_file, mode='w')
        fw.write(ssh_cmd)
        fw.close()
        cmd = 'ssh stheil@' + self.params['servers'][self.server]['adress'] + ' \'bash -s\' < ' + self.remote_cmd_file
        log.debug(cmd)
        self.cmd.append(cmd)
        cmd = 'scp stheil@' + self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch'] + '/' + os.path.basename(self.out) + ' ' + self.wd
        log.debug(cmd)
        self.cmd.append(cmd)


    def _get_exec_script (self):
        ssh_cmd = 'if [ -f ~/.bashrc ]; then' + "\n"
        ssh_cmd += 'source ~/.bashrc' + "\n"
        ssh_cmd += 'elif [ -f ~/.profile ]; then' + "\n"
        ssh_cmd += 'source ~/.profile' + "\n"
        ssh_cmd += 'elif [ -f ~/.bash_profile ]; then' + "\n"
        ssh_cmd += 'source /etc/profile' + "\n"
        ssh_cmd += 'source ~/.bash_profile' + "\n"
        ssh_cmd += 'else' + "\n"
        ssh_cmd += 'echo "source not found."' + "\n"
        ssh_cmd += 'fi' + "\n"
        ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + "\n"
        # ssh_cmd += 'mkdir ' + self.params['servers'][self.server]['scratch'] + '/' + self.sample + '_' + self.type + '_split' + "\n"
        ssh_cmd += 'echo "' + 'blast_launch.py -c ' + self.server + ' -n 100 --n_cpu ' + self.n_cpu + ' -d ' + self.params['servers'][self.server]['db'][self.db]
        ssh_cmd += ' -s ' + self.params['servers'][self.server]['scratch'] + '/' + os.path.basename(self.contigs) + ' --prefix ' + self.sample + '_' + self.type
        ssh_cmd += ' -p ' + self.type + ' -o ' + os.path.basename(self.out) + '"'
        qsub_cmd=''
        if self.server == 'enki':
            qsub_cmd = ' | ' + 'qsub -V -wd ' + self.params['servers'][self.server]['scratch'] + ' -N ' + self.sample
        elif self.server == 'avakas':
            # qsub_cmd = 'qsub -V -d ' + self.params['servers'][self.server]['scratch'] + ' -l walltime=48:00:00 -l nodes=1:ppn=1' + ' -N ' + self.sample
            qsub_cmd = ''
        elif self.server == 'genotoul':
            qsub_cmd = ' | ' + 'qsub -sync yes -V -wd ' + self.params['servers'][self.server]['scratch'] + ' -N ' + self.sample
        else:
            log.critical('unknown cluster.')
            log.debug(ssh_cmd)
        return ssh_cmd + qsub_cmd


    def check_args (self, args: dict):
        if 'sample' in args:
            self.sample = str(args['sample'])
        self.wd = os.getcwd() + '/' + self.sample
        self.cmd_file = self.wd + '/' + self.sample + '_blast_cmd.txt'
        self.remote_cmd_file = self.wd + '/' + self.sample + '_remote_blast_cmd.txt'
        accepted_type = ['blastx','blastn','blastp','blastn','rpstblastn']
        if 'contigs' in args:
            self.contigs = self.wd + '/' + args['contigs']
        if 'type' in args:
            if args['type'] in accepted_type:
                self.type = args['type']
            else:
                log.critical('Wrong blast type. ' + accepted_type)
        else:
            log.critical('Blast type is mandatory.')
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
        if 'params' in args:
            self.params = args['params']
        if 'server' in args:
            self.server = args['server']
        if 'db' in args:
            if args['db'] not in self.params['servers'][self.server]['db']:
                log.critical(arg['db'] + ' not defined in parameters file')
            else:
                self.db = args['db']
        else:
            log.critical('You must provide a database name.')



    def launch (self):
        if not self.sge:
            for el in self.cmd:
                os.system (el)
        else:
            fw =  open(self.cmd_file, mode='w')
            for el in self.cmd:
                fw.write(el + "\n")
            fw.close()
            qsub_call =   "qsub -wd " + self.wd + " -V -N " + self.sample + '_blast' + ' ' + self.cmd_file
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
