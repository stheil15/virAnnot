#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

# CREATE TABLE `cdd_taxo` (
# `id` int unsigned NOT NULL DEFAULT '0',
# `family` varchar(32) DEFAULT NULL,
# `genus` varchar(32) DEFAULT NULL,
# `no rank` varchar(32) DEFAULT NULL,
# `specie` varchar(32) DEFAULT NULL,
# `superkingdom` varchar(32) DEFAULT NULL
# );

my $fof='';
my $verbosity=1;
my $help;


GetOptions(
  "i|fof=s"               => \$fof,
	"v|verbosity=i"         => \$verbosity,
  "h|help"                => \$help,
);

if($help){
  &help;
}
if($verbosity > 1){
  Logger::Logger->changeMode($verbosity);
}

&main;

sub main {
	my $self={};
  bless $self;
	open(FOF,$fof);
  print 'begin transaction;' . "\n";
	while(<FOF>){
		chomp;
		my $file = $_;
		$self->_compute_frequency($file);
	}
	close FOF;
  print 'end transaction;' . "\n";
}


sub _compute_frequency {
	my ($self,$file)=@_;
	my @ranks = ('superkingdom','no rank','family','genus','specie');
	$file =~ /.*\/(.*)\.tax\.txt/;
	my $profile_name = $1;
	open(FILE,$file);
	my $h;
	while(<FILE>){
		chomp;
		my @line = split(/\t/,$_);
		if($line[1] ne 'unknown'){
			my @taxo = split(/;/,$line[1]);
			if(scalar(@taxo) == 4 || scalar(@taxo) == 5){
				for(my $i=0;$i<=$#taxo;$i++){
					$h->{$ranks[$i]}->{$taxo[$i]}++;
					$h->{$ranks[$i]}->{sum}++;
				}
			}
		}
		else{

		}
	}
	if(!defined($h)){
		return 0;
	}
	my $line1;
	my $line2;
	$line1 = 'insert into cdd_taxo ("id"';
	$line2 = ') values ("' . $profile_name . '"';

	foreach my $rank (sort {$a cmp $b} (keys(%{$h}))){
		my $x;
		foreach my $name (keys(%{$h->{$rank}})){
			if($name eq 'sum'){
				next;
			}
			my $ratio = $h->{$rank}->{$name}/$h->{$rank}->{sum};
			if($ratio >= 0.01){
				$x->{$name} = $h->{$rank}->{$name}/$h->{$rank}->{sum};
			}

		}
		my $nb=0;
		$line1 .= ',"' . $rank . '"';
		$line2 .= ',"';
		foreach my $n (sort {$x->{$b} <=> $x->{$a}}(keys(%{$x}))){
			$nb++;
			if($nb<=5){
				$line2 .= sprintf("%s(%.2f);", $n , $x->{$n});

			}
		}
		$line2 .= '"';
	}
	print $line1 . $line2 . ");\n";
	# print Dumper $h;
	close FILE;
}
