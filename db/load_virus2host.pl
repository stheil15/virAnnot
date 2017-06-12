#!/usr/bin/perl
use strict;
print 'begin transaction;' . "\n";
open(FILE,$ARGV[0]);
my $old_cmd = '';
while(<FILE>){
  chomp;
  if(/^#/){next;}
  my @l = split("\t",$_);
  my $cmd = 'insert into virus2hosts (representative, neighbor, taxonomy, lineage, host) values ("' . $l[0] . '", "' . $l[1] . '", "' . $l[4] . '", "' . $l[3] . '", "' . $l[2]  . '");';
  if($old_cmd eq $cmd){
    # print $old_cmd . "+++++\n";
    next;
  }
  else{
    print $cmd . "\n";
    $old_cmd = $cmd;
  }
}
print 'end transaction;' . "\n";
