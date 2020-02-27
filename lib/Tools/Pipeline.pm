package Tools::Pipeline;

use Cwd 'abs_path', 'cwd' ;
use Data::Dumper;
use Logger::Logger;
use File::Basename;
use List::Util 'any';

my @assembler_prefix = ('idba.pairs','idba.both','abt');
my @blast_prefix = ('bltn','bltx','tbltx','rps');
my @blast_extension = ('xml','csv');
my @blast_db = ('nt','nr','pfam','vir');
my @servers = ('enki','genotoul','avakas');

sub init {
	# my ($self,$map_file)=@_;
	my ($class,$params)=@_;
  my $self={};
  bless $self;
  _set_options($self,$params);
	$self->{map} = parse_map_file($self->{_map_file});
	# my $samples->{map} = parse_map_file($map_file);
  $self->{assembler_prefix} = \@assembler_prefix;
  $self->{blast_prefix} = \@blast_prefix;
  $self->{blast_extension} = \@blast_extension;
  $self->{blast_db} = \@blast_db;
  $self->{info} = build_dir_and_files($self);
	return $self;
}

sub _set_options {
  my ($self,$params)=@_;
  if(defined($params->{map_file})){
    if(-e $params->{map_file}){
      $self->{_map_file} = abs_path($params->{map_file});
    }
    else{
      $logger->error('Map file not found.' . $params->{map_file});
    }
  }
  else{
    $logger->error('You have to provide a sample mapping file.');
    exit(1);
  }

  if(defined($params->{assembler})){
    if(any { /$params->{assembler}/ } @assembler_prefix){
      $self->{_assembler} = $params->{assembler};
    }
    else{
      $logger->error('Wrong assembler prefix. ' . $params->{assembler});
    }
  }

  if(defined($params->{blast_p})){
    if(any { /$params->{blast_p}/ } @blast_prefix){
      $self->{_blast_p} = $params->{blast_p};
    }
    else{
      $logger->error('Wrong Blast prefix.');
    }
  }

  if(defined($params->{blast_db})){
    if(any { /$params->{blast_db}/ } @blast_db){
      $self->{_blast_db} = $params->{blast_db};
    }
    else{
      $logger->error('Wrong Blast database.');
    }
  }
  if(defined($params->{server})){
    if(any { /$params->{server}/ } @servers){
      $self->get_server_conf($params->{server});
    }
  }
}

sub get_server_conf {
	my ($self,$server)=@_;
	if($server eq 'enki'){
		$self->{_server}->{_db}->{nr} = '/media/db/ncbi/nr/nr';
		$self->{_server}->{_db}->{nt} = '/media/db/ncbi/nt/nt';
		$self->{_server}->{_db}->{vir} = '/media/db/ncbi/nt/viruses_nt.fasta';
    $self->{_server}->{_db}->{pfam} = '/media/data/db/ncbi/cdd/Pfam';
	}
	elsif($server eq 'avakas'){
		$self->{_server}->{_db}->{nr} = '/home/stheil/db/nr/nr';
		$self->{_server}->{_db}->{nt} = '/home/stheil/db/nt/nt';
		$self->{_server}->{_db}->{vir} = '/home/stheil/db/vir/vir.fna';
		$self->{_server}->{_adress} = 'avakas.mcia.univ-bordeaux.fr';
		$self->{_server}->{_scratch_dir} = '/scratch/stheil';
		$self->{_server}->{_prog}->{bltx} = 'blastx';
		$self->{_server}->{_prog}->{tbltx} = 'tblastx';
		$self->{_server}->{_prog}->{bltn} = 'blastn';
	}
	elsif($server eq 'genotoul'){
		$self->{_server}->{_db}->{nr} = '/bank/blastdb/nr';
		$self->{_server}->{_db}->{nt} = '/bank/blastdb/nt';
		$self->{_server}->{_db}->{vir} = '/save/stheil/db/vir/viruses_nt.fasta';
		$self->{_server}->{_adress} = 'genotoul.toulouse.inra.fr';
		$self->{_server}->{_scratch_dir} = '/work/stheil';
		$self->{_server}->{_prog}->{bltx} = 'blastx+';
		$self->{_server}->{_prog}->{tbltx} = 'tblastx+';
		$self->{_server}->{_prog}->{bltn} = 'blastn+';
		$self->{_server}->{_qsub} = '';
	}
	else{
		$logger->error($server . 'not found in conf.');
	}
}


sub build_dir_and_files {
	my ($self) = @_;
	my $cwd = cwd();
	my $hash={};
	foreach my $sampleID (keys(%{$self->{map}})){
		# Sample dir path.
		# if(-e $cwd . '/' . $sampleID){
			$hash->{$sampleID}->{path} = $cwd . '/' . $sampleID;
		# }
		# else{
		# 	$logger->warn('Sample directory ' . $cwd . '/' . $sampleID . ' not found.');
		# }
		my @split = split('_',basename($self->{map}->{$sampleID}->{file1}));
		$self->{map}->{$sampleID}->{lib_name} = $split[0];

		$hash->{$sampleID}->{reads}->{pair_R1} = $cwd . '/' . $sampleID . '/' . $sampleID . '_R1.truePairs.fastq';
		$hash->{$sampleID}->{reads}->{pair_R2} = $cwd . '/' . $sampleID . '/' . $sampleID . '_R2.truePairs.fastq';
		$hash->{$sampleID}->{reads}->{merged_pairs_fastq} = $cwd . '/' . $sampleID . '/' . $sampleID . '_merged.truePairs.fastq';

		$hash->{$sampleID}->{reads}->{merged_pairs_fasta} = $cwd . '/' . $sampleID . '/' . $sampleID . '_merged.truePairs.fasta';
		$hash->{$sampleID}->{reads}->{merged_singletons_fasta} = $cwd . '/' . $sampleID . '/' . $sampleID . '_merged.singletons.fasta';
		$hash->{$sampleID}->{reads}->{merged_singletons_fastq} = $cwd . '/' . $sampleID . '/' . $sampleID . '_merged.singletons.fastq';
		$hash->{$sampleID}->{reads}->{merged_reads_fasta} = $cwd . '/' . $sampleID . '/' . $sampleID . '_merged.reads.fasta';
		# Blast files
		foreach my $ass_p (@assembler_prefix){
			foreach my $db (@blast_db){
				foreach my $blt_p (@blast_prefix){
					foreach my $blt_e (@blast_extension){
						if($db =~ /^nt|vir$/ && $blt_p =~ /^tbltx|bltn$/){
							$hash->{$sampleID}->{assembly}->{$ass_p}->{blast_files}->{$db}->{$blt_p}->{$blt_e} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' . $ass_p . '.' . $blt_p . '.' . $db . '.' . $blt_e;
						}
						elsif($db =~ /^nr$/ && $blt_p =~ /^bltx$/){
							$hash->{$sampleID}->{assembly}->{$ass_p}->{blast_files}->{$db}->{$blt_p}->{$blt_e} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' . $ass_p . '.' . $blt_p . '.' . $db . '.' . $blt_e;
						}
						elsif($db =~ /^pfam$/ && $blt_p =~ /^rps$/){
							$hash->{$sampleID}->{assembly}->{$ass_p}->{blast_files}->{$db}->{$blt_p}->{$blt_e} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' . $ass_p . '.' . $blt_p . '.' . $db . '.' . $blt_e;
						}
						elsif($db =~ /^merge$/ && $blt_p =~ /^blt$/){
							$hash->{$sampleID}->{assembly}->{$ass_p}->{blast_files}->{$db}->{$blt_p}->{$blt_e} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' . $ass_p . '.' . $blt_p . '.' . $db . '.' . $blt_e;
						}
						else{

						}
					}
				}
				# output files

				$hash->{$sampleID}->{output}->{$ass_p}->{$db}->{krona} = $cwd . '/' . $sampleID . '/' . $sampleID . '_krona.' .$ass_p . '.' . $db . '.html';
				$hash->{$sampleID}->{output}->{$ass_p}->{$db}->{autoMapper} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' .$ass_p . '.' . $db  . '.autoMapper';
			}
      $hash->{$sampleID}->{output}->{$ass_p}->{pfam}->{rps2tree} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' . $ass_p . '.pfam.nr.rps2tree';
			$hash->{$sampleID}->{output}->{$ass_p}->{'merge'}->{krona} = $cwd . '/' . $sampleID . '/' . $sampleID . '_krona.' .$ass_p . '.merged.html';
			$hash->{$sampleID}->{output}->{'all'}->{'merge'}->{krona} = $cwd . '/' . $sampleID . '/' . $sampleID . '_krona.all.merged.html';
			$hash->{$sampleID}->{output}->{$ass_p}->{'merge.blt'}->{blast2ecsv} = $cwd . '/' . $sampleID . '/' . $sampleID . '_' . $ass_p . '.' . 'merge.blt.csv';
      $hash->{$sampleID}->{output}->{$ass_p}->{scaffolds} = $cwd . '/' . $sampleID . '/' . $sampleID . '_scaffold.' .$ass_p  . '.fna';
			# ABT & IDBA scaffolds
			$hash->{$sampleID}->{assembly}->{$ass_p}->{seq_files}->{scaffolds} = $cwd . '/' . $sampleID . '/' . $sampleID . '_scaffold.' .$ass_p  . '.fna';
			$hash->{$sampleID}->{assembly}->{$ass_p}->{seq_files}->{scaffolds_base_name} = $sampleID . '_scaffold.' .$ass_p;

			# ABT & IDBA mapping reads
			$hash->{$sampleID}->{assembly}->{$ass_p}->{mapping_files}->{blat_psl} = $cwd . '/' . $sampleID . '/' . $sampleID . '_scaffold.' .$ass_p  . '.psl';
			$hash->{$sampleID}->{assembly}->{$ass_p}->{mapping_files}->{bowtie_bam} = $cwd . '/' . $sampleID . '/' . $sampleID . '_scaffold.' .$ass_p  . '.bam';
			$hash->{$sampleID}->{assembly}->{$ass_p}->{mapping_files}->{bowtie_bam_sort} = $cwd . '/' . $sampleID . '/' . $sampleID . '_scaffold.' .$ass_p  . '.bam.sort.bam';
			$hash->{$sampleID}->{assembly}->{$ass_p}->{mapping_files}->{blat_rn} = $cwd . '/' . $sampleID . '/' . $sampleID . '_scaffold.' .$ass_p  . '.rn';
		}
		$hash->{$sampleID}->{output}->{excel} = $cwd . '/' . $sampleID . '/' . $sampleID . '_excel.xlsx';
		$hash->{$sampleID}->{output}->{demultiplex} = $cwd . '/' . $self->{map}->{$sampleID}->{lib_name} . '_demultiplex.stats.csv';
	}
	return $hash;
}


sub parse_map_file {
	my ($file)=@_;
	my $map={};
	open(MAP,$file);
	my $f_line = <MAP>;
	chomp $f_line;
	if($f_line !~ /^#/){
		$logger->error('Missing headers in map file. First line must stat with an \'#SampleID\'');
		exit;
	}
	$f_line =~ s/#//;
	my @headers = split("\t",$f_line);
	while(my $line = <MAP>){
		chomp $line;
		my @line = split("\t",$line);
		for(my $i=1;$i<=$#line;$i++){
			if($headers[$i] =~ /file/){
				$map->{$line[0]}->{$headers[$i]} = abs_path($line[$i]);
			}
			else{
				$map->{$line[0]}->{$headers[$i]} = $line[$i];
			}
		}
	}
	close MAP;
	return $map;
}


1;
