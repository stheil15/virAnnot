#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Logger::Logger;
use File::Basename;
use Cwd 'cwd';
use Tools::Mummer;
use Tools::Setter;

my $VERSION = "V3.0" ;
my $lastmodif = "2013-01-08" ;


my $ref='';
my $qry='';
my $mode='mum';
my $png =0;
my $length=0;
my $identity=0;
my $promer=0;
my $format='graph';
my $deltaMode = 'm';
my $matrix=2;
my $minMatch=0;
my $maxGap=0;
my $minCluster=0;
my $breaklen=0;
my $verbosity=1;
my $diagFactor = 0;
my $optimize = 1;

GetOptions( "r|ref=s"         => \$ref,
		        "q|qry=s"         => \$qry,
		        "l|length=i"      => \$length,
		        "i|identity=i"    => \$identity,
		        "p|promer"        => \$promer,
		        "m|mode=s"	      => \$mode,
		        "f|format=s"      => \$format,
		        "d|delta=s"       => \$deltaMode,
		        "x|matrix=i"      => \$matrix,
		        "s|minmatch=i"    => \$minMatch,
		        "g|maxgap=i"      => \$maxGap,
		        "c|mincluster=i"  => \$minCluster,
		        "b|breaklen=i"    => \$breaklen,
            "diagfactor=f"    => \$diagFactor,
            "optimize!"       => \$optimize,
		        "v|verbosity=i"   => \$verbosity,
		);


&main();


sub main {
	if($verbosity > 1){
		Logger::Logger->changeMode($verbosity);
	}
	if($ref eq '' || $qry eq ''){
		&help;
	}
	my $self = {} ;
	bless $self;
	$self->_setPrefix();

	$self->_checkOptions();

	$self->_executeAlignment();

	$self->_executeFiltering();

	$self->_executeRendering();
}

sub _executeRendering{
	my ($self) = @_ ;
	if($self->{format} =~ /^(graph|png)$/){
		$self->_executeMummerPlot();
	}
	elsif($self->{format} =~ /^(btab|tab|coords)$/){
		$self->_executeShowCoords();
	}
	elsif($self->{format} eq 'align'){
		$self->_executeShowAlign();
	}
	elsif($self->{format} eq 'fig'){
		$self->_executeMapView();
	}
	elsif($self->{format} eq 'fasta'){
		$self->_generateMultiFastaAlignmentFiles();
	}
}

sub _executeMummerPlot {
	my ($self) = @_ ;

	if($self->{format} eq 'graph'){
		$self->_setMummerPlot(0);

		`$self->{mummerPlotCmd}`;
	}
	if($self->{format} eq 'png'){
		$self->{pngFile} = $self->{prefix} . '.png';
		if(-e $self->{pngFile}){
			$logger->warn('PNG file with the same prefix already exist, skipping this step.');
		}
		else{
			$self->_setMummerPlot(1);
			`$self->{mummerPlotCmd}`;
		}
	}
}


sub _setMummerPlot {
	my ($self,$png) = @_ ;
	my $cmd = 'mummerplot ' ;
	$cmd .= ' -l ' ;
	$cmd .= ' -p ' . $self->{prefix};
	if($png == 1){
		$cmd .= ' -s large --png ';
	}
	if($self->{delta} eq '0'){
		$cmd .= ' ' . $self->{deltaFile} ;
	}
	else{
		$cmd .= ' ' . $self->{filterFile} ;
	}
	$logger->debug($cmd);
	$self->{mummerPlotCmd} = $cmd;
	$logger->debug($self->{mummerPlotCmd});
}


sub _executeShowCoords {
	my ($self) = @_ ;
	$self->{showCoordsFile} = $self->{prefix} . '.coords' ;
	if(-e $self->{showCoordsFile}){
		$logger->warn('Coords file already exist, skipping this step.');
	}
	else{
		$self->_setShowCoords();
		`$self->{showCoordsCmd}`;
	}
}

sub _setShowCoords {
	my ($self) = @_ ;
	my $cmd = 'show-coords -rl ';
	if($self->{format} eq 'coords'){
		$cmd .= ' -c ';
	}
	elsif($self->{format} eq 'btab'){
		$cmd .= ' -B ';
	}
	elsif($self->{format} eq 'tab'){
		$cmd .= ' -T ';
	}
	if($self->{delta} eq '0'){
		$cmd .= $self->{deltaFile}
	}
	else{
		$cmd .= $self->{filterFile}
	}
	$cmd .= ' > ' . $self->{showCoordsFile};
	$logger->debug($cmd);
	$self->{showCoordsCmd} = $cmd;
}

sub _executeShowAlign {
	my ($self) = @_ ;
	$self->_setShowAligns();
}

sub _setShowAligns {
	my ($self) = @_ ;
	$self->_setAllIds();
	my $res = '';
	foreach my $refId (@{$self->{refIds}}){
		foreach my $qryId (@{$self->{qryIds}}){
			my $cmd;
			if($self->{delta} eq '0'){
				$cmd = 'show-aligns -x ' . $self->{matrix} . ' -w 100 -r ' . $self->{deltaFile} . ' "' . $refId . '" "' . $qryId . '" 2> /dev/null';
			}
			else{
				$cmd = 'show-aligns -x ' . $self->{matrix} . ' -w 100 -r ' . $self->{filterFile} . ' "' . $refId . '" "' . $qryId . '" 2> /dev/null';
			}
			$res .= `$cmd`;
		}
	}
	open(OUT,'>' . $self->{prefix} . '.align');
	print OUT $res ;
	close OUT;
}

sub _generateMultiFastaAlignmentFiles {
	my ($self) = @_ ;
  if ($self->{format} eq 'fasta'){
		my $deltaFile;
		if($self->{delta} eq '0'){
			$deltaFile = $self->{deltaFile} ;
		}
		else{
			$deltaFile = $self->{filterFile} ;
		}
		if(! -z $deltaFile){
			my $mummer = Tools::Mummer -> new ({query => $self->{query}, reference => $self->{reference}, alignments =>  [$deltaFile]});
			my @files = $mummer->generateMultiFastaAlignment('fasta');
			open(FILE, '> fastaList');
			print FILE join ("\n", @files);
			close FILE
		}
	}
}

sub _setAllIds {
	my ($self) = @_ ;
	$self->{refIds} = $self->_getIdsFromFasta($self->{reference});
	$self->{qryIds} = $self->_getIdsFromFasta($self->{query});
}

sub _getIdsFromFasta {
	my ($self,$file) = @_ ;
	my @t;
	open(FILE,$file) || $logger->logdie('Cannot open file '. $file);
	while(<FILE>){
		chomp;
		if(/^>(\S+)/){
			push(@t,$1);
		}
	}
	if(scalar(@t) > 0){
		return \@t ;
	}
	else{
		$logger->logdie('Void array of IDs.');
	}
	close FILE;
}


sub _executeMapView {
	my ($self) = @_ ;
	$self->{mapViewFile} = $self->{prefix} . '.fig' ;
}


sub _executeFiltering {
	my ($self) = @_ ;
	$self->{filterFile} = $self->{prefix} . '.delta.filter' ;
	if($self->{delta} eq '0'){
		$logger->warn('Filtering option set to "0", skipping this step.');
	}
	elsif(-e $self->{filterFile}){
		$logger->warn('Filtering file already exist, skipping this step.');
	}
	else{
		$self->_setFiltering();
		`$self->{fileteringCmd}`;
	}
}


sub _setFiltering {
	my ($self) = @_ ;
	my $cmd = 'delta-filter ' ;
	$cmd .= '-' . $self->{delta} ;
	$cmd .= ' -i ' . $self->{identity} ;
	$cmd .= ' -l ' . $self->{hspLength} ;
	$cmd .= ' ' . $self->{deltaFile} ;
	$cmd .= ' > ' . $self->{filterFile} ;
	$logger->debug($cmd) ;
	$self->{fileteringCmd} = $cmd ;
}


sub _executeAlignment {
	my ($self) = @_ ;
	$self->{deltaFile} = $self->{prefix} . '.delta' ;
	if(-e $self->{deltaFile}){
		$logger->warn('Delta file with same prefix already exist, skipping alignment step.');
	}
	else{
		$self->_setAlignmentCmd();
		`$self->{alignCmd}`;
	}
}


sub _setAlignmentCmd {
	my ($self) = @_ ;
	my $cmd ;
	if($self->{promer} == 1){
		$cmd = 'promer ' ;
	}
	else{
		$cmd = 'nucmer ' ;
	}
	$cmd .= '--' . $self->{mode} ;
	$cmd .= ' -prefix ' . $self->{prefix} ;
	$cmd .= ' ' . $self->{reference} ;
	$cmd .= ' ' . $self->{query} ;
	$cmd .= ' -minmatch ' . $self->{minMatch} ;
	$cmd .= ' -maxgap ' . $self->{maxGap} ;
	$cmd .= ' -mincluster ' . $self->{minCluster} ;
	$cmd .= ' -breaklen ' . $self->{breaklen};
    $cmd .= ' -diagfactor ' . $self->{diagFactor};

    if($self->{optimize}){
        $cmd .= ' --optimize ';
    }
    else { $cmd .= ' --nooptimize ' }
	$logger->debug($cmd);

	if($verbosity == 1){
		$cmd .= ' 2> /dev/null';
	}
	$self->{alignCmd} = $cmd ;
}


sub _setPrefix{
	my ($self) = @_ ;
	$self->{prefix} = fileparse($ref, "\.[^.]*") . '_' . fileparse($qry, "\.[^.]*") ;
	if($identity != 0){
		$self->{prefix} .= '_i' . $identity ;
	}
	if($length != 0){
		$self->{prefix} .= '_l' . $length ;
	}
	if($matrix != 2){
		$self->{prefix} .= '_m' . $matrix ;
	}
	if($minMatch != 0){
		$self->{prefix} .= '_s' . $minMatch ;
	}
	if($maxGap != 0){
		$self->{prefix} .= '_g' . $maxGap ;
	}
	if($minCluster != 0){
		$self->{prefix} .= '_c' . $minCluster ;
	}
	if($breaklen != 0){
		$self->{prefix} .= '_b' . $breaklen ;
	}
	if($deltaMode ne '1' ){
		$self->{prefix} .= '_d' . $deltaMode ;
	}
	$logger->debug('Prefix for this alignment will be: ' . $self->{prefix});
}


sub _checkOptions {
	my ($self) = @_ ;
	$self->{reference} = setFile($ref);
	$self->{query} = setFile($qry);
	$self->{hspLength} = setPositiveInteger($length);
	$self->{identity} = setPourcentage($identity);
	$self->{promer} = $promer;
	$self->{mode} = $self->_setMode($mode);
	$self->{format} = $self->_setFormat($format);
	$self->{delta} = $self->_setDelta($deltaMode);
	$self->{matrix} = $self->_setPromerMatrix($matrix);
	$self->{minMatch} = $self->_setMinMatch($minMatch);
	$self->{maxGap} = $self->_setMaxGap($maxGap);
	$self->{minCluster} = $self->_setMinCluster($minCluster);
	$self->{breaklen} = $self->_setBreaklen($breaklen);
    $self->{diagFactor} = $self->_setDiagFactor($diagFactor);
    $self->{optimize} = $self->_setOptimize($optimize);
}


sub _setBreaklen {
	my ($self,$breaklen) = @_;
	$breaklen = setPositiveInteger($breaklen);
	if($self->{promer} == 1){
		if($breaklen == 0){
			return 60;
		}
		else{
			return $breaklen;
		}
	}
	else{
		if($breaklen == 0){
			return 200;
		}
		else{
			return $breaklen;
		}
	}
}


sub _setMinCluster {
	my ($self,$minCluster) = @_ ;
	$minCluster = setPositiveInteger($minCluster);
	if($self->{promer} == 1){
		if($minCluster == 0){
			return 20;
		}
		else{
			return $minCluster;
		}
	}
	else{
		if($minCluster == 0){
			return 65;
		}
		else{
			return $minCluster;
		}
	}
}


sub _setMaxGap {
	my ($self,$maxGap) = @_ ;
	$maxGap = setPositiveInteger($maxGap);
	if($self->{promer} == 1){
		if($maxGap == 0){
			return 30;
		}
		else{
			return $maxGap;
		}
	}
	else{
		if($maxGap == 0){
			return 90;
		}
		else{
			return $maxGap;
		}
	}
}


sub _setMinMatch {
	my ($self,$minMatch) = @_ ;
	$minMatch = setPositiveInteger($minMatch);
	if($self->{promer} == 1){
		if($minMatch == 0){
			return 6;
		}
		else{
			return $minMatch;
		}
	}
	else{
		if($minMatch == 0){
			return 20;
		}
		else{
			return $minMatch;
		}
	}
}


sub _setPromerMatrix {
	my ($self,$matrix) = @_ ;
	if($matrix < 1 || $matrix > 3){
		$logger->error('Matrix parameter out of range.');
		&help();
	}
	else{
		return $matrix ;
	}
}


sub _setDelta {
	my ($self,$deltaMode) = @_ ;
	if($deltaMode !~ /^(1|m|r|q|g|0)$/){
		$logger->error('Wrong delta mode selected.');
		&help;
	}
	else{
		return $deltaMode;
	}
}


sub _setFormat {
	my ($self,$format) = @_ ;
	if($format !~ /^(graph|png|coords|tab|btab|align|fasta)$/){
		$logger->error('Wrong format selected.');
		&help;
	}
	else{
		return $format;
	}
}


sub _setMode {
	my ($self,$mode) = @_ ;
	if($mode !~ /^(mum|mumreference|maxmatch)$/){
		$logger->error('Wrong mode selected.');
		&help;
	}
	else{
		return $mode;
	}
}

sub _setDiagFactor {
	my ($self,$diagFactor) = @_ ;
    if($diagFactor >= 0){
        if($self->{promer} == 1){
            if($diagFactor == 0){
                return 0.11;
            }
            else{
                return $diagFactor;
            }
        }
        else{
            if($diagFactor == 0){
                return 0.12;
            }
            else{
                return $diagFactor;
            }
        }
    }
    else{
        $logger->error('Wrong diagfactor value provided : '.$diagFactor.'. Must be a positive number');
        &help;
    }
}

sub _setOptimize {
	my ($self,$optimize) = @_ ;
    if($optimize == 0 || $optimize == 1){
        return $optimize;
    }
    else{
        $logger->error('Wrong optimize value provided : '.$optimize.'. Must be 0 or 1');
        &help;
    }
}

sub help {
	my $prog = basename($0) ;
print STDERR <<EOF ;
#### $prog ####
#
# CREATED:    2011-03-04
# LAST MODIF: $lastmodif
# AUTHOR:     Sebastien THEIL (INRA Bordeaux)
# MAIL:       sebastien theil at bordeaux inra fr
# VERSION:    $VERSION
#
# This script is just a pipeline that makes the use of the Mummer package easier.
# It can represente the alignment in different formats, like plots (graph, png), arrays (btab,
# tabular, coords), alignments (align) or figures (fig).
#
# For futher explainations see MUMMER website http://mummer.sourceforge.net/

USAGE:
       $prog  [OPTIONS] -r ref.fna -q query.fna

       -r|ref         Fasta file of the reference sequence(s).

       -q|qry         Fasta file of the query sequence(s).

       -p|promer      Amino acid alignments between two mutli-FASTA DNA input.

       ### ALIGNMENT OPTIONS ###

       -m|mode        mum, mumreference, maxmatch. (Default $mode)

       -x|matrix      Set the alignment matrix number to 1 [BLOSUM 45], 2 [BLOSUM62]
                      or 3 [BLOSUM 80]. (default 2) Available only for promer!

       -s|minmatch    Set the minimum length of a single match.
                      (default 20 for nucmer, 6 for promer)

       -g|maxgap      Set the maximum gap between two adjacent matches in a cluster.
                      (default 90 for nucmer, 30 for promer)

       -c|mincluster  Sets the minimum length of a cluster of matches.
                      (default 65 for nucmer, 20 for promer)

       -b|breaklen    Set the distance an alignment extension will attempt to extend poor scoring regions before giving up.
                      (default 200 for nucmer, 60 for promer)

       -diagfactor    Set the clustering diagonal difference separation factor
                      (default 0.12 for nucmer, 0.11 for promer)

       -[no]optimize  Toggle alignment score optimization, i.e. if an alignment
                      extension reaches the end of a sequence, it will backtrack
                      to optimize the alignment score instead of terminating the
                      alignment at the end of the sequence (default --optimize)


       ### FILTERING OPTIONS ###

       -l|length      Set the minumum length of an HSP. (default $length)

       -i|identity    Set the minimum identity of an HSP. (default $identity)

       -d|delta       1 : 1-to-1 alignment allowing for rearrangements.
                          (intersection of -r and -q alignments)
                      g : 1-to-1 global alignment not allowing rearrangements.
                      m : Many-to-many alignment allowing for rearrangements. (default)
                          (union of -r and -q alignments)
                      q : Maps each position of each query to its best hit in the reference,
                          allowing for reference overlaps.
                      r : Maps each position of each reference to its best hit in the query,
                          allowing for query overlaps.
                      0 : No filtering. Skip delta-filter stage.

       ### OUTPUT OPTIONS ###

       -f|format      graph, png, btab, tab, coords, align, fasta, fig. (default graph)

EOF
exit(1) ;
}
