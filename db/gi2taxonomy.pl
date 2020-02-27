#!usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use DBI;
use Tools::Blast;
use Tools::Taxonomy;
use Tools::Fasta;
use Pod::Usage;
use Logger::Logger;
use File::Basename;

my $VERSION = '0.1' ;
my $lastmodif = '2016-01-09' ;

my %blastOptions;
my $db = '/media/scratch/taxo/taxonomy.tmp.sqlite';
my $dbType = 'ncbi';
my $verbosity = 1;
my @inputFiles;
my $outputFile = 'blast_ecsv.csv';
my $help;
my $reduceTaxo;

GetOptions(
  "i|input=s"             => \@inputFiles,
  "o|output=s"            => \$outputFile,
  "db=s"                  => \$db,
  "r|reduce"              => \$reduceTaxo,
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
  open(OUT,">$outputFile");
  $self->_set_options();
    foreach my $file (@{$self->{inputFiles}}){
        open(FILE,$file);
        while (<FILE>) {
            chomp;
      my $acc;
      # if line start with >gi
      if ($_ =~ />gi(.*)/){
        my @id = split('\|', $_);
        $acc = $id[1];
        print OUT $acc . "\t" . $self->{taxoTools}->gi2taxonomy($acc,'BLASTX') . "\n";
      }
        }
    exit;
    }
  close OUT;
}


sub _set_options {
  my ($self) = @_;
  @{$self->{inputFiles}} = map{glob $_} @inputFiles;
  if(scalar(@{$self->{inputFiles}}) > 0){
    foreach my $file (@{$self->{inputFiles}}){
      if(! -e $file){
        $logger->error('gi file not found. ' . $file);
        &help;
      }
    }
  }
  else{
    $logger->error('You must provide at least one Blast file.');
    &help;
  }


  my %taxonomyParams;
  if(! -e $db){
    $logger->error('Taxonomy sqlite database not found. ' . $db);
  }
    else{
        $taxonomyParams{'dbh'} = DBI->connect("dbi:SQLite:dbname=$db","","");
    }
  if($reduceTaxo){
    $taxonomyParams{_dbFormat} = 'reduced';
  }
  else{
    $taxonomyParams{_dbFormat} = 'normal';
  }
  $self->{taxoTools} = Tools::Taxonomy->new(%taxonomyParams);

}
