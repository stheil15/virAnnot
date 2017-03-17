#!/usr/bin/perl
use strict;
use Data::Dumper;
use Logger::Logger;
use File::Basename;
use Getopt::Long;
use Cwd 'abs_path' ;
use Tools::Fasta ;
use Tools::Fastq ;


my $ref='';
my @p1='';
my @p2='';
my @read='';
my $identity = 95;
my $coverage = 60;
my $verbosity=1;
my $format;
my $stat='nb_reads';
my $blatOutput='';
my $statOutput='';
my $map_prog='';

GetOptions(
	"r=s"              => \$ref,
	"id=i"             => \$identity,
	"cov=i"            => \$coverage,
	"b|blat=s"	       => \$blatOutput,
	"o|output=s"       =>	\$statOutput,
	"f=s"		           => \$format,
	"s=s"		           => \$stat,
	"p|prog=s"         => \$map_prog,
	"v|verbosity=i"    => \$verbosity
);

Logger::Logger->changeMode($verbosity);

&main;

sub main {
	my $self = {};
	bless $self;
	_setOptions($self);

	$self->{ref_tool} = Tools::Fasta->new(file => $self->{reference});

	if($map_prog eq 'blat'){
		_parseBlat($self);
	}
	elsif($map_prog eq 'bowtie'){
		_parseBam($self);
	}
	else{
		$logger->logdie('Wrong mapping program name.');
	}

	_print($self,$self->{_statistics}, $self->{statOutput});
}


sub _parseBam {
	my ($self) = @_;
	$logger->debug('Parsing BAM file.');
	foreach my $ref (sort(keys(%{$self->{ref_tool}->{index}}))){
		my $cmd = 'samtools view ' . $self->{blatOutput} . ' ' . $ref . ' | wc -l';
		my $res = `$cmd`;
		chomp $res;
		$self->{readNumber}->{$ref} = $res;
	}
}


sub _parseBlat {
	my ($self) = @_;
	open(PSL,$self->{blatOutput}) || $logger->logdie('Cant read psl file ' . $self->{blatOutput});
	my $already_seen={};
	while(my $line = <PSL>){
		chomp $line;
		my @line = split(/\t/,$line);
		if(defined($already_seen->{$line[13]})){
			next;
		}
		$already_seen->{$line[13]}++;
		my $tSize = $line[14];
		my $tStart = $line[15];
		my $tEnd = $line[16];
		my $identity = (100 *($line[0]+$line[2])/($line[0]+$line[1]+$line[2]));
		my $coverage = (($tEnd-$tStart)/$tSize) * 100;
		$logger->debug($line[13] . ' id: ' . $identity . ' cov: ' . $coverage);
		if($coverage >= $self->{_covThreshold} && $identity >= $self->{_idThreshold}){
			$self->{readNumber}->{$line[9]}++;
			$self->{totalReadNumber}++;
			$self->{size}->{$line[9]} = $line[10];
			$self->{reads_size_sum}->{$line[9]} += $tSize;
		}
		else{
			next;
		}
	}
	close PSL;
}


sub _print {
	my ($self, $mode, $statOutput) = @_;
	$logger->info('Printing results ...');

	if(!defined $mode){
		$mode = 'nb_reads';
	}

	open(OUT,">$statOutput");
	foreach my $ref (sort(keys(%{$self->{ref_tool}->{index}}))){
		my $value = 0;
		if(defined($self->{readNumber}->{$ref})){
			if($mode eq 'nb_reads'){
				$value = $self->{readNumber}->{$ref};
			}
			elsif($mode eq 'rpkm'){
				$value = $self->{readNumber}->{$ref} / ($self->{size}->{$ref}/1000 * $self->{totalReadNumber}/1000000);
			}
			elsif($mode eq 'mean_coverage'){
				$value =  $self->{reads_size_sum}->{$ref} / $self->{size}->{$ref};
			}
			elsif($mode eq 'bp_size'){
				$value =  $self->{reads_size_sum}->{$ref} ;
			}
		}
		print OUT $ref . "\t" . $value . "\n";
	}
	close OUT;
}


sub _setOptions {
	my ($self) = @_;
	if($map_prog !~ /blat|bowtie/){
		$logger->error('You must provide a valid mapping program name. blat|bowtie');
		&help;
	}
	if(-e $ref){
		$self->{reference} = abs_path($ref);
	}
	else{
		$logger->error('You must provide a reference.');
		&help;
	}
	if(-e $blatOutput){
		$self->{blatOutput} = abs_path($blatOutput);
	}
	else{
		$logger->error('You must provide a blat file.' . $blatOutput);
		&help;
	}

	if($identity >= 0 && $identity <= 100){
		$self->{_idThreshold} = $identity;
	}
	else{
		$logger->logdie('ERROR: identity must be between 0 and 100.');
	}
	if($coverage >=0 && $coverage <= 100){
		$self->{_covThreshold} = $coverage;
	}
	else{
		$logger->logdie('ERROR: coverage must be between 0 and 100.');
	}
	if($stat =~ /nb_reads|rpkm|mean_coverage|bp_size/){
		$self->{_statistics} = $stat;
	}
	else{
		$logger->logdie('Wrong option. Please choose between nb_reads|rpkm|mean_coverage');
	}
	if($statOutput eq ''){
		$self->{statOutput} = $self->{reference} . '.rn';
	}
	else{
		$self->{statOutput} = $statOutput;
	}
}


sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
#### $prog ####
#
# AUTHOR:     Sebastien THEIL
# LAST MODIF: 16/07/2015
# PURPOSE:    This script is basically used to map reads (using blat) on contigs (or scaffolds) and count the
              number of read per contig.
              Input can be multiple fasta or fastq.
              Output is a 2 column tabular file.

USAGE:
      $prog  -i singl.fastq -i singl.fasta -1 R1.fastq -2 R2.fastq ....

      ### OPTIONS ###
      -r          <string>   Rererence sequences file.
      -id         <int>      Identity threshold. (Default: $identity)
      -cov        <int>      Coverage threshold. (Default: $coverage)
      -p|prog     <blat|bowtie>
      -b          <string>   Blat output file path
      -o          <string>   Output file path
      -s          <string>   Statistics to print : nb_reads, rpkm or mean_coverage.
      -v          <int>      Verbosity level. (0 -> 4).
EOF
exit(1);
}
