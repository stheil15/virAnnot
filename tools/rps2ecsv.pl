#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use DBI;
use Logger::Logger;
use File::Basename;
use Tools::Taxonomy;
use Bio::SearchIO;

my $VERSION = '0.1';
my $lastmodif = '2016-01-11';

my $blastFile='';
my $outputFile = 'rpsblast_ecsv.csv';
my $max_evalue=0.001;
my $parse_description=0;
my $help;
my $verbosity = 1;
my $db = '/media/data/db/taxonomy/taxonomy.sqlite';

GetOptions(
	"b|blast=s"             => \$blastFile,
	"o|output=s"            => \$outputFile,
	"e|evalue=f"            => \$max_evalue,
	"pd"                    => \$parse_description,
	"v|verbosity=i" 	=> \$verbosity
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
	_set_options($self);
	_readInputFile($self,$blastFile);
	my $filtered_matches = _check_overlapping_domains($self,$self->{_matches});
	$self->{_matches} = $filtered_matches;
	_print($self);
}


sub _print {
	my ($self) = @_;
	open(OUT,">$outputFile");
	my @fields = ('query_id','query_length','cdd_id','hit_id','evalue','startQ','endQ','frame','description','superkingdom','no rank','family','genus');
	print OUT '#' . join("\t",@fields) . "\n";
	foreach my $query_id (sort(keys(%{$self->{_matches}}))){
		my $found=0;
		foreach my $m (@{$self->{_matches}->{$query_id}}){
			if(! defined($m->{no_hit})){
				$found=1;
				for(my $i=0;$i<=$#fields;$i++){
					if(!defined($m->{$fields[$i]})){
						$logger->error($fields[$i] . ' not defined.');
						$logger->error(Dumper $self->{_matches}->{$query_id});
						exit;
					}
					print OUT '"' . $m->{$fields[$i]} . '"';
					if($i<$#fields){
						print OUT "\t";
					}
					else{
						print OUT "\n";
					}
				}
			}
		}
		if($found==0){
			print OUT '"' . $query_id . '"' . "\t" . 'no_hit' . "\n";
		}
	}
}


sub _check_overlapping_domains {
	my ($self,$matches) = @_;
	my $locus;
	my $selected;
	my $treated;
	my $grabed;
	foreach my $q_id (keys(%{$matches})){
		# $logger->debug('treating ' . $q_id);
		if(scalar(@{$matches->{$q_id}}) == 1){
			# $logger->debug('only one match...');
			push(@{$grabed->{$q_id}},$matches->{$q_id}->[0]);
			if(defined($matches->{$q_id}->[0]->{no_hit})){
				$selected->{$q_id . '-no_hit'}++;
			}
			else{
				$selected->{$q_id . '-' . $matches->{$q_id}->[0]->{startQ} . '-' . $matches->{$q_id}->[0]->{endQ} . '-' . $matches->{$q_id}->[0]->{cdd_id}}++;
			}
			next;
		}
		# print Dumper $matches->{$q_id} ;
		for(my $i=0;$i<scalar(@{$matches->{$q_id}});$i++){
			my $found=0;
			my $overlapped;
			if(defined($matches->{$q_id}->[$i]->{no_hit})){
				next;
			}

			if(defined($selected->{$q_id . '-' . $matches->{$q_id}->[$i]->{startQ} . '-' . $matches->{$q_id}->[$i]->{endQ} . '-' . $matches->{$q_id}->[$i]->{cdd_id}})){
				next;
			}
			for(my $j=0;$j<scalar(@{$matches->{$q_id}});$j++){
				my $pair1 = $q_id . '-' . $matches->{$q_id}->[$i]->{startQ} . '-' . $matches->{$q_id}->[$i]->{endQ} . '-' . $matches->{$q_id}->[$i]->{cdd_id} . '-' . $matches->{$q_id}->[$j]->{startQ} . '-' . $matches->{$q_id}->[$j]->{endQ} . '-' . $matches->{$q_id}->[$j]->{cdd_id};
				my $pair2 = $q_id . '-' . $matches->{$q_id}->[$j]->{startQ} . '-' . $matches->{$q_id}->[$j]->{endQ} . '-' . $matches->{$q_id}->[$j]->{cdd_id} . '-' . $matches->{$q_id}->[$i]->{startQ} . '-' . $matches->{$q_id}->[$i]->{endQ} . '-' . $matches->{$q_id}->[$i]->{cdd_id};
				if($pair1 eq $pair2){
					next;
				}
				if(defined($treated->{$pair1}) || $treated->{$pair2}){
					# $logger->debug($pair1 . ' already compared');
					next;
				}
				$logger->debug('comparing ' . $matches->{$q_id}->[$i]->{cdd_id} . ' ' . $matches->{$q_id}->[$i]->{startQ} . ' ' . $matches->{$q_id}->[$i]->{endQ} . ' ' . $matches->{$q_id}->[$j]->{cdd_id} . ' ' . $matches->{$q_id}->[$j]->{startQ} . ' ' . $matches->{$q_id}->[$j]->{endQ});
				if(defined($matches->{$q_id}->[$i]->{'no_hit'})){
					push(@{$grabed->{$q_id}},$matches->{$q_id}->[$i]);
					next;
				}

				if(defined($selected->{$q_id . '-' . $matches->{$q_id}->[$j]->{startQ} . '-' . $matches->{$q_id}->[$j]->{endQ} . '-' . $matches->{$q_id}->[$j]->{cdd_id}})){
					next;
				}
				$treated->{$pair1}++;
				if($matches->{$q_id}->[$i]->{startQ} <= $matches->{$q_id}->[$j]->{endQ} && $matches->{$q_id}->[$i]->{endQ} >= $matches->{$q_id}->[$j]->{startQ}){
					$found=1;
					$selected->{$q_id . '-' . $matches->{$q_id}->[$i]->{startQ} . '-' . $matches->{$q_id}->[$i]->{endQ} . '-' . $matches->{$q_id}->[$i]->{cdd_id}}++;
					push(@{$overlapped},$matches->{$q_id}->[$i]);
					$selected->{$q_id . '-' . $matches->{$q_id}->[$j]->{startQ} . '-' . $matches->{$q_id}->[$j]->{endQ} . '-' . $matches->{$q_id}->[$j]->{cdd_id}}++;
					push(@{$overlapped},$matches->{$q_id}->[$j]);
					$logger->debug('overlapping');
				}
			}
			if($found==0){
				push(@{$grabed->{$q_id}},$matches->{$q_id}->[$i]);
				$logger->debug('no overlap found.');
			}
			else{
				my @sorted = sort {$a->{evalue} <=> $b->{evalue}} @{$overlapped};
				# print Dumper @sorted;
				push(@{$grabed->{$q_id}},$sorted[0]);
				$logger->debug('overlap found.');
			}
		}
	}
	return $grabed;
	# print Dumper $selected;
}


sub _readInputFile {
	my ($self,$file) = @_;
  $logger->debug('Reading input xml file.');
	$self->{_searchioObj} = Bio::SearchIO->new(   '-file'   => $file,
	'-format' => 'BLASTXML');
	my @matches;
	while ($self->{_result} = $self->{_searchioObj}->next_result()){
		my $oldQueryId='';
		if($self->{_result}->num_hits() == 0){
			next;
		}
		while ($self->{_currentHit} = $self->{_result}->next_hit() || ! $self->{_result}->num_hits()){
			my $match = {};
			if($oldQueryId ne $self->{_result}->query_name()){
				$oldQueryId = $self->{_result}->query_name();
			}
			if($parse_description == 1){
				$match->{query_id} = $self->{_result}->query_description();
			}
			else{
				$match->{query_id} = $self->{_result}->query_name();
			}
			$match->{query_length} = $self->{_result}->query_length();
			$match->{description} = $self->{_currentHit}->description();
			my @s = split(/,/,$match->{description});
			$match->{cdd_id} = $s[0];
			$match->{hit_id} = $self->{_currentHit}->name;
			my $startQ=9999999999;
			my $endQ=0;
			my @hsps = $self->{_currentHit}->hsps();

			foreach ( @hsps ) {
				$self->{_hsp} = $_;
        if($self->{_hsp}->evalue > $max_evalue){
          next;
        }
        else{
          my %match_hsp = %$match;
          $match_hsp{frame} = $self->{_hsp}->{QUERY_FRAME};
          $match_hsp{evalue} = $self->{_hsp}->evalue;
          $match_hsp{startQ} = $self->{_hsp}->start('query');
          $match_hsp{endQ} = $self->{_hsp}->end('query');
          $match_hsp{score} = $self->{_hsp}->score;
          _get_domain_taxonomy($self,\%match_hsp);
          push(@{$self->{_matches}->{$match_hsp{query_id}}},\%match_hsp);
        }
			}
      if(!defined($self->{_matches}->{$match->{query_id}})){
        $match->{no_hit}=1;
        push(@{$self->{_matches}->{$match->{query_id}}},$match);
      }
      # if(scalar(@{$self->{_matches}->{$match->{query_id}}}) == 0){
      #   $match->{no_hit}=1;
      # }
		}
	}
  $logger->info('End reading file.');
	return \@matches;
}

sub _get_domain_taxonomy {
  my ($self,$match)=@_;
  my $cmd = 'select * from cdd_taxo where id="' . $match->{cdd_id}. '"';
  my $sth=$self->{taxoTools}->{dbh}->prepare($cmd);
  $sth->execute();
  # print $cmd . "\n";
  my $res = $sth->fetchall_hashref('id');
  if(scalar(keys(%{$res})) > 0){
    $match->{'superkingdom'} = $res->{$match->{cdd_id}}->{'superkingdom'};
    $match->{'no rank'} = $res->{$match->{cdd_id}}->{'no rank'};
    $match->{'specie'} = $res->{$match->{cdd_id}}->{'specie'};
    $match->{'family'} = $res->{$match->{cdd_id}}->{'family'};
    $match->{'genus'} = $res->{$match->{cdd_id}}->{'genus'};
  }
  else{
    $match->{'superkingdom'} = 'unknown';
    $match->{'no rank'} = 'unknown';
    $match->{'specie'} = 'unknown';
    $match->{'family'} = 'unknown';
    $match->{'genus'} = 'unknown';
  }
}

sub _set_options {
  my ($self) = @_;

	if(! -e $blastFile){
    $logger->error('You must provide at least one Blast file.');
    &help;
  }
  my %taxonomyParams;
  if(! -e $db){
    $logger->error('Taxonomy sqlite database not found. ' . $db);
		&help();
  }
  else{
		$taxonomyParams{'dbh'} = DBI->connect("dbi:SQLite:dbname=$db","","");
		$self->{taxoTools} = Tools::Taxonomy->new(%taxonomyParams);
  }
}


sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    This script is used to parse blast+ output and produce a m8 output
#             extended with taxonomy information.
#

USAGE: perl $prog -i blast_1 -i blast_2 ... -i blast_n [OPTIONS]

          ### OPTIONS ###
          -b|blast                 <BLAST>  CSV or XML blast files
          -o|output                <ECSV>  Output ECSV file. Default: blast_ecsv.csv

          ### Blast parsing options ###
          -e|evalue                <FLOAT>  E-value threshold
          -pd                     Parse description.
          -v|verbosity            Verbosity level
          -h|help                 Print this help and exit

EOF
exit(1) ;
}
