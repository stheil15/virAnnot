#!/usr/bin/perl

use strict ;
use Data::Dumper ;
use Cwd ;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;

my $PATH = '' ;

GetOptions(
	"dir=s" 	=> \$PATH,
);

&_readDir($PATH);
	
sub _readDir {
	my $path = $_[0] ;
	opendir(PATH, $path) || die "Can't open directory $path" ;
	my @files = readdir(PATH) ;
	foreach my $file (sort(@files)){
		if($file =~ /^\.\.?$/){next;}
		$file = $path . "/" . $file ;
		if(-l $file){
			print $file . "\n" ;
			next;
		}
        $file  = abs_path($file);
		if( -f $file ){ 
			print $file . "\n" ;
		}
		if( -d $file ){
			print $file . "\n" ;
			&_readDir($file);
		}
	}
	close PATH ;
}
