#!/usr/bin/perl

use strict;
use Data::Dumper;
use Logger::Logger;
use Getopt::Long;
use Tools::Fastq;
use Cwd 'abs_path', 'cwd' ;
use File::Basename ;
use String::Random;

my $VERSION = '3.0' ;
my $lastmodif = '2013-4-30' ;

my $pair1='';
my $pair2='';
my $indexFile='';
my $adapter='';
my $common={};
my $minLength=0;
my $quality=0;
my $pairsOutput = 'split';
my $force;
my $out='';
my $middle=0;
my $minIndexMatchPercent=0.4;
my $verbosity=1;
my $testMode=0;
my $index={};
my $polyA;
my $tmp_file_prefix='';



GetOptions(
            "-1=s"             => \$pair1,
            "-2=s"             => \$pair2,
            "indexFile=s"      => \$indexFile,
						"i|index=s%"       => \$index,
						"c|common=s%"      => \$common,
            "a|adapters=s"     => \$adapter,
            "l|len=i"          => \$minLength,
            "q|qual=i"         => \$quality,
            "po|pairsOutput=s" => \$pairsOutput,
            "f|force"          => \$force,
            "middle=i"         => \$middle,
            "o|outdir=s"       => \$out,
            "polyA"            => \$polyA,
            "test"             => \$testMode,
            "tmp_prefix=s"     => \$tmp_file_prefix,
			      "v|verbosity=i"    => \$verbosity,
);
sub help() {
  my $prog = basename($0) ;
  print STDERR <<EOF ;
  ### $prog $VERSION ###
  #
  # AUTHOR:     Sebastien Theil
  # VERSION:    $VERSION ($lastmodif)
  #
  # This script is used to demultiplex illumina reads using cutadapt
  # and trim off the adaptator sequence. It support paired end reads
  # and single reads.

  USAGE:
  Paired End Reads:
  $prog -1 R1.fq -2 R2.fq -i index.fa -c common.fa

  ### OPTIONS ###
  -1               Paired End reads 1.
  -2               Paired End reads 2.
  -indexFile       Index file sequences in fasta format.
  -i|index         Index tag and sequences in TAG=SEQUENCE format.
  -a|adapters      Illumina adapters sequences to trim.
	-c|common        Common sequence adapter in TAG=SEQUENCE format.
  -l|len           Minimum length of read to keep, after trimming.
  -q|qual          Minimum quality value.
  -pf|pairsFormat  Pairs format (split|merge). Separates pairs in two file or merge them.
  -po|pairsOutput  Pairs output (split|merge). Separates pairs depending on where they come. truePairs: same index. falsePairs: another index and singleton.
	-f|force
	-middle [1|2]    Search for common tag in the middle of the read. 1: trim the read. 2: exclude the read.
	-test            Test mode, just print the commands and don't execute them.
	-tmp_prefix      The prefix to name the tmp files. If none random one will be generated.
  -v|verbosity
  -h

EOF
  exit(1) ;
}

if($verbosity > 1){
    Logger::Logger->changeMode($verbosity);
}


&main;


sub main {
	my $self = {};
	bless $self;
	_setOptions($self);
	_launchPairedReadStep($self);
	exit 0;
}


sub _step_01 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1}, $self->{index}, $tmp_file_prefix . "_step.01_R1",'d','k','-g','0','1','0.8');
  $self->{dispatch}->{$files->{1}} = _readInfoFile($self,$tmp_file_prefix . "_step.01_R1.info");
  _launchCutAdapt($self,$files->{2}, $self->{index}, $tmp_file_prefix . "_step.01_R2",'d','k','-g','0','1','0.8');
  $self->{dispatch}->{$files->{2}} = _readInfoFile($self,$tmp_file_prefix . "_step.01_R2.info");
	return ($tmp_file_prefix . '_step.01_R1.out',$tmp_file_prefix . '_step.01_R2.out');
}


sub _step_02 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1},$self->{_common},$tmp_file_prefix . "_step.02_R1",'k','k','-g','0.1',scalar(keys(%{$self->{_common}})),'0.7');
  _launchCutAdapt($self,$files->{2},$self->{_common},$tmp_file_prefix . "_step.02_R2",'k','k','-g','0.1',scalar(keys(%{$self->{_common}})),'0.7');
	return ($tmp_file_prefix . '_step.02_R1.out',$tmp_file_prefix . '_step.02_R2.out');
}


sub _step_021 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1},$self->{_common},$tmp_file_prefix . "_step.021_R1",'k','k','-g','0.2',scalar(keys(%{$self->{_common}})),'0.5');
  _launchCutAdapt($self,$files->{2},$self->{_common},$tmp_file_prefix . "_step.021_R2",'k','k','-g','0.2',scalar(keys(%{$self->{_common}})),'0.5');
	return ($tmp_file_prefix . '_step.021_R1.out',$tmp_file_prefix . '_step.021_R2.out');
}


sub _step_03 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1},$self->{illuminaAdapter},$tmp_file_prefix . "_step.03_R1",'k','k','-b','0.2',scalar(keys(%{$self->{illuminaAdapter}})),'0.6');
  _launchCutAdapt($self,$files->{2},$self->{illuminaAdapter},$tmp_file_prefix . "_step.03_R2",'k','k','-b','0.2',scalar(keys(%{$self->{illuminaAdapter}})),'0.6');
	return ($tmp_file_prefix . '_step.03_R1.out',$tmp_file_prefix . '_step.03_R2.out');
}


sub _step_04 {
  my ($self,$files,$tmp_file_prefix) = @_;
	my $h;
	foreach my $x (keys(%{$self->{_common}})){
		$h->{$x} = $self->{_common}->{$x} ;
		$h->{$x . '-REVCOMP'} = _reverseComplement($self,$self->{_common}->{$x});
	}
	if($self->{middle_mid} == 1){
		_launchCutAdapt($self,$files->{1},$h,$tmp_file_prefix . "_step.04_R1",'k','k','-b','0.1','1','0.5');
	  _launchCutAdapt($self,$files->{2},$h,$tmp_file_prefix . "_step.04_R2",'k','k','-b','0.1','1','0.5');
	}
  elsif($self->{middle_mid} == 2){
		_launchCutAdapt($self,$files->{1},$h,$tmp_file_prefix . "_step.04_R1",'k','d','-b','0.1','1','0.5');
	  _launchCutAdapt($self,$files->{2},$h,$tmp_file_prefix . "_step.04_R2",'k','d','-b','0.1','1','0.5');
	}
	else{
		$logger->error('Wrong middle_mid value. ' . $self->{middle_mid});
		exit 1 ;
	}
	return ($tmp_file_prefix . '_step.04_R1.out',$tmp_file_prefix . '_step.04_R2.out');
}

sub _step_05 {
  my ($self,$files,$tmp_file_prefix) = @_;
  my $h->{polyA} = 'AAAAAAAAAA';
  _launchCutAdapt($self,$files->{1}, $h, $tmp_file_prefix . "_step.05_R1",'k','k','-a','0','1','0.8');
  _launchCutAdapt($self,$files->{2}, $h, $tmp_file_prefix . "_step.05_R2",'k','k','-a','0','1','0.8');
  return ($tmp_file_prefix . '_step.05_R1.out',$tmp_file_prefix . '_step.05_R2.out');
}


sub _reverseComplement {
	my ($self,$seq) = @_;
	my $revComp = reverse($seq);
	$revComp =~ tr/atcgnATCGN/tagcnTAGCN/;
	return $revComp;
}


sub _launchPairedReadStep {
	my ($self) = @_;
	my $i=0;
  my $files;
  $files->{1} = $self->{readFiles}->{1};
  $files->{2} = $self->{readFiles}->{2};

  # STEP 01: demultiplex from index file.
  if(defined($self->{index})){
    ($files->{1},$files->{2}) = _step_01($self,$files,$self->{tmp_file_prefix});
  }
  # STEP 02: 5' trimming common part.
	if(defined($self->{_common})) {
    ($files->{1},$files->{2}) = _step_02($self,$files,$self->{tmp_file_prefix});
		($files->{1},$files->{2}) = _step_021($self,$files,$self->{tmp_file_prefix});
	}

  # STEP 03: 5' and 3' (reverse complement) trimming illumina adapters
	if(defined($self->{illuminaAdapter})) {
    ($files->{1},$files->{2}) = _step_03($self,$files,$self->{tmp_file_prefix});
	}

  # STEP 04: 3' trimming or excluding reads with common part (reverse complemented). (trimming chimeric reads)
	if(defined($self->{_common}) && $self->{middle_mid}>0) {
    ($files->{1},$files->{2}) = _step_04($self,$files,$self->{tmp_file_prefix});
	}

  if(defined($polyA)){
    ($files->{1},$files->{2}) = _step_05($self,$files,$self->{tmp_file_prefix});
  }

  #STEP 05: printing last fastq
	my $last_R1_file = $files->{1};
	my $last_R2_file = $files->{2};
	$last_R1_file = abs_path($last_R1_file);
	$last_R2_file = abs_path($last_R2_file);
	if($self->{_testMode}==0){
		$self->{sequences}->{$last_R1_file} = Tools::Fastq->new(file => $last_R1_file);
		$self->{sequences}->{$last_R2_file} = Tools::Fastq->new(file => $last_R2_file);

		if(! defined($self->{index})){
			$self->{dispatch}->{$self->{readFiles}->{1}} = _create_dispatch_hash($self,$self->{sequences}->{$last_R1_file});
			$self->{dispatch}->{$self->{readFiles}->{2}} = _create_dispatch_hash($self,$self->{sequences}->{$last_R2_file});
		}

		_linkPairs($self,$self->{readFiles}->{1},$self->{readFiles}->{2},$last_R1_file,$last_R2_file);
		_linkPairs($self,$self->{readFiles}->{2},$self->{readFiles}->{1},$last_R2_file,$last_R1_file);

		_printStats($self);

		_dispatchPairs($self,$last_R1_file,$last_R2_file);

	}
}


sub _create_dispatch_hash {
	my ($self,$file) = @_;
	my $h = {};
	foreach my $readName (keys($file->{index})){
		$h->{indexToRead}->{$self->{tmp_file_prefix}}->{$readName} = 1 ;
		$h->{readToIndex}->{$readName}->{$self->{tmp_file_prefix}} = 1;
	}
	return $h;
}


sub _printStats {
	my ($self) = @_;
  $logger->debug('Printing statistics...');
	my $printed=0;
	my $stat_file = $self->{tmp_file_prefix} . '_demultiplex.stats.csv';
	open(STATS,">$stat_file");
	foreach my $file (sort(keys(%{$self->{found}}))){
		print STATS $file . "\n";
		foreach my $key (sort(keys(%{$self->{found}->{$file}}))){
      if($key eq 'index'){
        foreach my $index (sort(keys(%{$self->{found}->{$file}->{index}}))){
          print STATS $index . "\t" . $self->{found}->{$file}->{index}->{$index}->{seq} . "\t" . $self->{found}->{$file}->{index}->{$index}->{nb} . "\n";
        }
      }
			else{
        print STATS $key . "\t" . $self->{found}->{$file}->{$key} . "\n";
      }
		}
		print STATS "\n";
		print STATS "\n";
	}

	# my @statusList = (sort (keys(%{$self->{stats}})));
	my @statusList;
	push(@statusList,'singleton');
	push(@statusList,'truePairs');

	foreach my $indexName (sort (keys(%{$self->{stats}}))){
		if($printed==0){
			print STATS "\t" . join ("\t", @statusList);
			print STATS "\n";
			$printed=1;
		}
		print STATS $indexName ;
		foreach my $status (@statusList){
			if(defined($self->{stats}->{$indexName}->{$status})){
				if($status eq 'singleton'){
					my $sum=0;
					foreach my $file (keys(%{$self->{stats}->{$indexName}->{$status}})) {
						$sum += scalar(keys(%{$self->{stats}->{$indexName}->{$status}->{$file}} )) ;
					}
					print STATS "\t" . $sum ;
				}
				else{
					print STATS "\t" . scalar(keys(%{$self->{stats}->{$indexName}->{$status}})) ;
				}
			}
			else{
				print STATS "\t" . '0';
			}
		}
		print STATS "\n";
	}
	close STATS;
}


sub _dispatchPairs {
  my ($self,$file1,$file2) = @_;
  my $file1_tool = Tools::Fastq->new(file => abs_path($file1));
  my $file2_tool = Tools::Fastq->new(file => abs_path($file2));
  $logger->info('Dispatching reads...');
  foreach my $indexName (keys(%{$self->{stats}})){
		$logger->info('Dispatching reads for index ' . $indexName);
		#TODO: change for outputDir option.
		my $path = $self->{_outDir} . '/' . $indexName;
		if(-e $path){
			$logger->warn('Directory ' . $path . ' already exists.');
			if(defined($self->{_force})){
				$logger->warn('Force mode...');
				`rm $path`;
			}
			else{
				next;
			}
		}
		else{
			`mkdir $path`;
		}
		chdir $path;
		my $pwd = abs_path('.');
		$logger->info($pwd);

		$self->{_outputFileHanle}->{truePairs} = _openFileHandles($self,$indexName, 'truePairs');
		$self->{_outputFileHanle}->{singleton} = _openFileHandles($self,$indexName, 'singletons');

		foreach my $status (keys(%{$self->{stats}->{$indexName}})){
			$logger->debug('Treating status ' . $status);
			if ($status eq 'truePairs') {
				my @readList = keys($self->{stats}->{$indexName}->{$status});
				if(scalar(@readList) > 0){
					_printPairs($self,\@readList,$self->{_outputFileHanle}->{truePairs},$file1_tool,$file2_tool);
				}
			}
			elsif($status eq 'singleton'){
				foreach my $file (keys(%{$self->{stats}->{$indexName}->{$status}})) {
          # print $file . "\n";
					my @readList = keys(%{$self->{stats}->{$indexName}->{$status}->{$file}});
					if(scalar(@readList) > 0){
            if($file eq $file1){
              _printSingleton($self,\@readList,$self->{_outputFileHanle}->{singleton}, $file1_tool);
            }
            else{
              _printSingleton($self,\@readList,$self->{_outputFileHanle}->{singleton}, $file2_tool);
            }
					}
				}
			}
			else{
				$logger->logdie('Unknown status ' . $status);
			}
		}
		_closeFileHandles($self,$self->{_outputFileHanle});
		chdir '..';
		$pwd = abs_path('.');
		$logger->info($pwd);
  }
}


sub _printSingleton {
	my ($self,$readList,$fhs,$file) = @_;
	$logger->debug('Retrieving singleton from ' . $file->{file} . '.');
	printf({$fhs->[0]} "%s",$file->retrieveFastqBlock($readList)) ;
}


sub _printPairs {
	my ($self,$readList,$fhs,$file1,$file2) = @_;
	if(scalar(@{$fhs}) == 2){
		if(defined($self->{_separator})){
			$logger->debug('Retrieving pairs from ' . $file1->{file} . '.');
			my @a = map {$_ . $self->{_separator} . '1'} @{$readList};
			printf({$fhs->[0]} "%s",$file1->retrieveFastqBlock(\@a)) ;
			$logger->debug('Retrieving pairs from ' . $file2->{file} . '.');
			@a = map {$_ . $self->{_separator} . '2'} @{$readList};
			printf({$fhs->[1]} "%s",$file2->retrieveFastqBlock(\@a)) ;
		}
		else{
			$logger->debug('Retrieving pairs from ' . $file1->{file} . '.');
			printf({$fhs->[0]} "%s",$file1->retrieveFastqBlock($readList)) ;
			$logger->debug('Retrieving pairs from ' . $file2->{file} . '.');
			printf({$fhs->[1]} "%s",$file2->retrieveFastqBlock($readList)) ;
		}
	}
	else{
		foreach my $id (@{$readList}){
			$logger->debug('Retrieving pairs from ' . $file1->{file} . '.');
			printf({$fhs->[0]} "%s",$file1->retrieveFastqBlock($id)) ;
			$logger->debug('Retrieving pairs from ' . $file2->{file} . '.');
			printf({$fhs->[0]} "%s",$file2->retrieveFastqBlock($id)) ;
		}
	}
}


sub _closeFileHandles {
	my ($self,$h) = @_;
	foreach my $k (keys(%{$h})){
		foreach my $fh (@{$h->{$k}}){
			close $fh;
		}
	}
}


sub _openFileHandles {
	my ($self,$prefix,$suffix) = @_;
	if($self->{_pairsFormat} eq 'split' && $suffix ne 'singletons'){
		my ($fh1,$fh2) = ();
		open($fh1,">$prefix" . "_R1." . "$suffix" . ".fastq");
		open($fh2,">$prefix" . "_R2." . "$suffix" . ".fastq");
		return [$fh1,$fh2];
	}
	else{
		my $fh1;
		open($fh1,">$prefix" . "_merged." . "$suffix" . ".fastq");
		return [$fh1];
	}
}


sub _linkPairs {
  my ($self,$file1,$file2,$last_fq1,$last_fq2) = @_;
  my $fqIndex1 = $self->{sequences}->{$last_fq1}->{index};
  my $fqIndex2 = $self->{sequences}->{$last_fq2}->{index};
  $logger->info('Linking pairs between ' . $file1 . ' and ' . $file2);
  foreach my $indexName (keys(%{$self->{dispatch}->{$file1}->{indexToRead}})){
		foreach my $readName (keys(%{$self->{dispatch}->{$file1}->{indexToRead}->{$indexName}})){
			my $read_base_name = '';
			my $pair1 = '';
			my $pair2 = '';
			my $read1 = '';
			my $read2 = '';
			if($readName =~ /(.*)([\.\/])([12])/){
				$read_base_name = $1;
				$self->{_separator} = $2;
				if($3 eq '1'){
					$pair1 = 1;
					$pair2 = 2;
				}
				else{
					$pair1 = 2;
					$pair2 = 1;
				}
				$read1 = $read_base_name . $self->{_separator} . $pair1;
				$read2 = $read_base_name . $self->{_separator} . $pair2;
			}
			else{
				$read_base_name = $readName;
				$read1 = $readName;
				$read2 = $readName;
			}
			if (defined($self->{stats}->{$indexName}->{truePairs}->{$read_base_name})) {
				next;
			}
			# print $read_base_name . "\t" . $separator . "\t" . $pair1 . "\t" . $pair2 . "\n";

			# if the read exist in the second fastq file.
			if (defined($self->{dispatch}->{$file2}->{readToIndex}->{$read2})) {
				# if it has the same indexName.
				# if($indexName eq (keys(%{$self->{dispatch}->{$file2}->{readToIndex}->{$readName}}))[0] ){
				if( defined($self->{dispatch}->{$file2}->{readToIndex}->{$read2}->{$indexName}) ){
					# if the read still exist in both last fastq file.
					if( defined($fqIndex1->{$read1}) &&  defined($fqIndex2->{$read2})){
						$self->{stats}->{$indexName}->{truePairs}->{$read_base_name} = 1;
					}
					else{
						if(defined($fqIndex1->{$read1})){
							$self->{stats}->{$indexName}->{singleton}->{$last_fq1}->{$read1} = 1;
						}
						if(defined($fqIndex2->{$read2})){
							$self->{stats}->{$indexName}->{singleton}->{$last_fq2}->{$read2} = 1;
						}
					}
				}
				# else pairs come from two different indexes.
				else{
					# if the read still exist in both last fastq file.
					if(defined($fqIndex1->{$read1})){
						$self->{stats}->{$indexName}->{singleton}->{$last_fq1}->{$read1} = 1;
					}
					if(defined($fqIndex2->{$read2})){
						$self->{stats}->{ (keys(%{$self->{dispatch}->{$file2}->{readToIndex}->{$read2}}))[0] }->{singleton}->{$last_fq2}->{$read2} = 1;
					}
				}
			}
			else{
				$self->{stats}->{$indexName}->{singleton}->{$last_fq1}->{$read1} = 1;
			}
		}
  }
}


sub _readInfoFile {
	my ($self,$file) = @_;
	my $h={};
	if($self->{_testMode}== 1){
		return $h;
	}
	open(FILE,$file) || $logger->logdie('Info file not found ' . $file );
	$logger->debug('Reading info file ' . $file);
	while(<FILE>){
		chomp;
		my @line = split(/\t/,$_);
		if(scalar(@line) <= 4){
			next;
		}
		my $match = { readName   => $line[0],
		              errorNb    => $line[1],
		              start      => $line[2],
		              end        => $line[3],
		              bLength    => length($line[4]),
		              indexMatch => $line[5],
		              aLength    => length($line[6]),
		              indexName  => $line[7]
		             };
		if($match->{readName} =~ /(\S+)\s\S+/){
			$match->{readName} = $1;
		}

		$h->{indexToRead}->{ $match->{indexName} }->{ $match->{readName} } = 1;
		$h->{readToIndex}->{ $match->{readName}  }->{ $match->{indexName} } = 1;
	}
	return $h
}


sub _launchCutAdapt {
	my ($self,$file,$index,$prefix,$discardUntrimmed,$discardTrimmed,$indexLocation,$error,$count,$minIndexMatchPercent) = @_;
	if(-e $prefix . '.out' && -e $prefix . '.info'){
		$logger->info('Skipping cutadapt execution...');
		_readLogFile($self,$prefix . '.log');
		return 0;
	}
	else{
		my $max_length=999;
		$self->{_cutAdaptCmd} = 'cutadapt -f fastq' ;
		foreach my $p (keys(%{$index})){
			$self->{_cutAdaptCmd} .= ' ' . $indexLocation . ' ' . $p . '=' . $index->{$p} ;
			if(length($index->{$p}) < $max_length){
				$max_length = length($index->{$p});
			}
		}
		$self->{_cutAdaptCmd} .= ' -o ' . $prefix . '.out' .
		                         ' --info-file ' . $prefix . '.info';
		if($discardUntrimmed eq 'd'){
			$self->{_cutAdaptCmd} .= ' --discard-untrimmed ' ;
		}
		if($discardTrimmed eq 'd'){
			$self->{_cutAdaptCmd} .= ' --discard-trimmed ' ;
		}
		if(defined($self->{_quality})){
			$self->{_cutAdaptCmd} .= ' -q ' . $self->{_quality} . ' --quality-base=' . $self->{_qualityBase} ;
		}
		if(defined($self->{_minLength})){
			$self->{_cutAdaptCmd} .= ' -m ' . $self->{_minLength} ;
		}
		if(defined($count) && $count > 0){
			$self->{_cutAdaptCmd} .= ' -n ' . $count;
		}
		my $overlapLength = int($max_length*$minIndexMatchPercent) ;
		$self->{_cutAdaptCmd} .= ' -O ' . $overlapLength;
		$self->{_cutAdaptCmd} .= ' -e ' . $error .
		                         ' ' . $file;
		$self->{_cutAdaptCmd} .= ' > ' . $prefix . '.log';
		$logger->debug($self->{_cutAdaptCmd});
		if(defined($self->{sge})){
			return $self->{_cutAdaptCmd};
		}
		else{
			if($self->{_testMode}==0){
				system($self->{_cutAdaptCmd});
				_readLogFile($self,$prefix . '.log');
			}
		}
	}
}


sub _readLogFile {
	my ($self,$file) = @_;
	$logger->debug('Reading log file ' . $file);
	open(LOG,$file) || $logger->logdie('Cannot open file ' . $file);
  my $seqName;
	while(<LOG>){
		chomp;
    if(/^Total reads processed:\s+(\S+)$/){
      $self->{found}->{$file}->{total_reads} = $1;
      $self->{found}->{$file}->{total_reads} =~ s/,//g;
    }
    if(/^Total basepairs processed:\s+(\S+)\sbp$/){
      $self->{found}->{$file}->{total_bases} = $1;
      $self->{found}->{$file}->{total_bases} =~ s/,//g;
    }
    if(/^Reads with adapters:\s+(\S+)\s\(\S+%\)$/){
      $self->{found}->{$file}->{reads_with_adapters} = $1;
      $self->{found}->{$file}->{reads_with_adapters} =~ s/,//g;
    }
    if(/^Reads written \(passing filters\):\s+(\S+)\s\(\S+%\)$/){
      $self->{found}->{$file}->{reads_ok} = $1;
      $self->{found}->{$file}->{reads_ok} =~ s/,//g;
    }
    if(/^Quality-trimmed:\s+(\S+)\sbp\s\(\S+%\)$/){
      $self->{found}->{$file}->{bases_quality_trimmed} = $1;
      $self->{found}->{$file}->{bases_quality_trimmed} =~ s/,//g;
    }
    if(/^Total written \(filtered\):\s+(\S+)\sbp\s\(\S+%\)$/){
      $self->{found}->{$file}->{bases_ok} = $1;
      $self->{found}->{$file}->{bases_ok} =~ s/,//g;
    }
    if(/^Reads that were too short:\s+(\S+)\s\(\S+\)$/){
      $self->{found}->{$file}->{reads_too_short} = $1;
      $self->{found}->{$file}->{reads_too_short} =~ s/,//g;
    }
    if(/^=== Adapter (\S+) ===$/){
      $seqName = $1;
    }
		if(/^Sequence: (\S+); Type: .*; Length: \d+; Trimmed: (\d+) times\./){
			$self->{found}->{$file}->{index}->{$seqName}->{nb} = $2;
			$self->{found}->{$file}->{index}->{$seqName}->{seq} = $1;
		}
	}
	close LOG;
}


sub _readIndexFile {
	my ($self,$file) = @_;
	open(FILE,$file) || $logger->logdie('File not found.');
	my $id ='';
	my $seq='';
	my $h;
	while(<FILE>){
		chomp;
		if(/^>(\S+)/){
			chomp;
			if($seq ne ''){
				$h->{$id} = $seq;
				$seq='';
			}
			$id=$1;
		}
		else{
			$seq.=$_;
		}
	}
	$h->{$id} = $seq;
	return $h;
}


#~ ------- TOOLS ------- ~#


sub _setOptions {
	my ($self) = @_;
	$self->{wd} = cwd;
	$self->{_qualityBase} = 33;

	if($pair1 ne ''){
		if(-e $pair1){
			$self->{readFiles}->{1} = abs_path($pair1);
		}
		else{
			$logger->error('Pair1 file ' . $pair1 . ' do not exists.');
			&help;
		}
		$self->{_pairs} = 1;
	}

	if($pair2 ne ''){
		if(-e $pair2){
			$self->{readFiles}->{2} = abs_path($pair2);
		}
		else{
			$logger->error('Pair2 file ' . $pair2 . ' do not exists.');
			&help;
		}
		$self->{_pairs} = 1;
	}

	if($pair1 eq '' && $pair2 eq ''){
		$logger->error('You must provide at least one read file.');
		&help;
	}

	$self->{_quality} = $quality;
	$self->{_minLength} = $minLength;

	if($adapter ne ''){
		$self->{illuminaAdapter} = $self->_readIndexFile($adapter);
		$self->{_illuminaAdapterFile} = abs_path($indexFile);
	}

	if(scalar(keys(%{$common})) > 0){
		$self->{_common} = $common;
	}

	if(scalar(keys(%{$index})) > 0){
		$self->{index} = $index;
	}
	elsif($indexFile ne '' && scalar(keys(%{$index})) != 0){
		$self->{_indexFile} = abs_path($indexFile);
		$self->{index} = $self->_readIndexFile($self->{_indexFile});
	}
	else{
		$logger->warn('No index file provided.');
	}

	if($pairsOutput !~ /split|merge/ ){
		$logger->error('Wrong pairsOutput option. split|merge');
		&help;
	}
	else{
		$self->{_pairsFormat} = $pairsOutput;
	}

	if(defined($force)){
		$self->{_force} = 1;
	}

	if(defined($out)){
		$self->{_outDir} = abs_path($out);
	}
	else{
		$self->{_outDir} = abs_path('.');
	}

	if($middle != 0){
		$self->{middle_mid} = $middle;
	}
	else{
		$self->{middle_mid} = undef;
	}
	$self->{_minIndexMatchPercent} = $minIndexMatchPercent;
  if($testMode==1){
    $self->{_testMode} = 1 ;
  }
	else{
		$self->{_testMode} = 0 ;
	}
  if($tmp_file_prefix ne ''){
    $self->{tmp_file_prefix} = $tmp_file_prefix;
  }
  else{
    my $foo = new String::Random;
  	$self->{tmp_file_prefix} = $foo->randregex('\w\w\w\w');
  }
}
