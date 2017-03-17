#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use DBI;
use Cwd qw (abs_path cwd);
use File::Path qw(make_path);
use File::Copy;
use Tools::Fasta;
use Tools::Taxonomy;
use Tools::GoogleChart::Chart;
use File::Basename;
use Logger::Logger;

my $VERSION = '1.2';
my $lastmodif = '2015-10-27';

my $help;
my $verbosity=1;

my %mapperOptions = (
	'queries' 			    => '',
	'taxonomyDatabase' 	=> '/media/data/db/taxonomy/taxonomy.sqlite',
	'blastDatabase' 	  => '/media/data/db/ncbi/nt/nt',
	'blast'				      => '',
	'outputDir' 		    => 'autoMapper_out',
	'filter' 			      => 'Viruses',
	'windowsSize' 		  => 50,
);


my %mummerOptions = (
	'm' 			=> 'mum',
	'format'	=> 'align',
);


GetOptions(
	\%mummerOptions,
	"help|h"            => \$help,
	"r|ref=s"           => \$mapperOptions{'blastDatabase'},
	"reads=s"           => \$mapperOptions{'reads_file'},
	"q|qry=s"           => \$mapperOptions{'queries'},
	"b|blast=s"         => \$mapperOptions{'blast'},
	"t|taxonomy=s"      => \$mapperOptions{'taxonomyDatabase'},
  "o|output=s"        => \$mapperOptions{'outputDir'},
  "w|windows=i"       => \$mapperOptions{'windowsSize'},
  "f|filter=s"        => \$mapperOptions{'filter'},
	"v|verbosity=i"	    => \$verbosity,
	"p|promer",
	"m|mode=s",
	"x|matrix=i",
	"s|minmatch=i",
	"g|maxgap=i",
	"c|mincluster=i",
	"breaklen=i",
	"diagfactor=f",
	"optimize!",
	"l|length=i",
	"i|identity=i",
	"d|delta=s",
	"format=s",
);


if($verbosity){
    Logger::Logger->changeMode($verbosity);
}

&main;

sub main {
	my $self={};
	bless $self;
	$self->{_mapping_file_headers} = ['accession','tax_id','fasta_file','read_map_file','mat_file','organism','taxonomy','description'];
	$self->processOptions(\%mapperOptions,\%mummerOptions);
	my $query = Tools::Fasta->new(file => $mapperOptions{queries});
	my $taxid_infos = $self->readExtentedBlast($mapperOptions{blast}, $mapperOptions{filter});
	foreach my $tax_id (keys(%{$taxid_infos})){
		$self->{_map}->{$tax_id}->{taxonomy} = $taxid_infos->{$tax_id}->{taxonomy};
		$self->{_map}->{$tax_id}->{organism} = $taxid_infos->{$tax_id}->{organism};

		$taxid_infos->{$tax_id}->{acc_file} = $mapperOptions{outputDirAcc} . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '.csv';

		if(! $self->accListFromTaxid($tax_id, $taxid_infos->{$tax_id}->{acc_file}, \%mapperOptions)) {
			$logger->warn('no accession found, skipping...');
			next;
		}

		$taxid_infos->{$tax_id}->{fasta_ref} = $mapperOptions{outputDirFasta} . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '.fasta';
		if(! $self->fastaFromAccList($mapperOptions{blastDatabase}, $taxid_infos->{$tax_id}->{fasta_ref}, $taxid_infos->{$tax_id}->{acc_file})){
			next;
		}
		$taxid_infos->{$tax_id}->{mummerDirectory} = $mapperOptions{outputDirMummer} . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism};
		$self->createFiles([$taxid_infos->{$tax_id}->{mummerDirectory}], 1);
		_get_sequence_info($self,$tax_id,$taxid_infos->{$tax_id}->{fasta_ref});
		$mummerOptions{format} = 'fasta';
		$self->mummerMapping($taxid_infos->{$tax_id}->{fasta_ref}, $mapperOptions{queries}, \%mummerOptions, $taxid_infos->{$tax_id}->{mummerDirectory});
		foreach my $file (glob ($taxid_infos->{$tax_id}->{mummerDirectory} . '/fasta/*.fasta')){
			if(-z $file){
				next;
			}
      my @t = split('\|',$file);
      my $acc=$t[3];

			if(!defined($acc)){
				$logger->error('Unable to parse accession from ' . $acc);
				exit;
			}

			move($file,$taxid_infos->{$tax_id}->{mummerDirectory} . '/' . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '_' . $acc .'.fasta');
			$self->{_map}->{$tax_id}->{$acc}->{fasta_file} = $taxid_infos->{$tax_id}->{mummerDirectory} . '/' . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '_' . $acc .'.fasta';
      print "#######################" . "\n";
      print $self->{_map}->{$tax_id}->{$acc}->{fasta_file} . "\n";
		}

		$mummerOptions{format} = 'align';
		$self->mummerMapping($taxid_infos->{$tax_id}->{fasta_ref}, $mapperOptions{queries}, \%mummerOptions, $taxid_infos->{$tax_id}->{mummerDirectory});

		foreach my $aln_file (glob ($taxid_infos->{$tax_id}->{mummerDirectory} . '/*.align')){
			if(-z $aln_file){
				next;
			}
			my $hash = _googleChartIdentityFromMummerDirectory($self, $aln_file, $mapperOptions{windowsSize},$taxid_infos->{$tax_id},$tax_id );
			# $self->{_map}->{$tax_id}->{$gi}->{mat_file} = $self->googleChartIdentityFromMummerDirectory( $aln_file, $mapperOptions{windowsSize},$taxid_infos->{$tax_id},$tax_id );
			# $self->{_map}->{$gi}->{aln_file} = $aln_file;
		}
		$mummerOptions{format} = 'btab';
		$self->mummerMapping($taxid_infos->{$tax_id}->{fasta_ref}, $mapperOptions{'reads_file'}, \%mummerOptions, $taxid_infos->{$tax_id}->{mummerDirectory});
		foreach my $btab_file (glob ($taxid_infos->{$tax_id}->{mummerDirectory} . '/*.coords')){
			if(-z $btab_file){
				next;
			}
			# $self->{_map}->{$tax_id}->{$gi}->{btab_file} = $btab_file;
			my $mapping_matrix = $taxid_infos->{$tax_id}->{mummerDirectory} . '/' . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '_reads' . '.mat';
			my $matches = _read_mapping_coords($self,$btab_file,$mapping_matrix);
			_print_mapping_matrix($self,$matches,$taxid_infos->{$tax_id}->{mummerDirectory},$tax_id,$taxid_infos->{$tax_id}->{organism});
		}
	}
	if(scalar(keys(%{$self->{_map}})) != 0){
		my $map_file_name = $mapperOptions{'outputDirResults'} . '/' . 'map.txt';
		_print_mapping_file($self,$map_file_name);
		my $cmd = 'autoMapper_matrix_to_html.py';
		$cmd .= ' -m ' . $map_file_name;
		$cmd .= ' > ' . $self->{_sample_id} . '_autoMapper.html';
		$logger->debug($cmd);
		chdir($mapperOptions{'outputDirResults'});
		system($cmd);
	}
}


sub _get_sequence_info {
	my ($self,$tax_id,$file)=@_;
	my $ref_file = Tools::Fasta->new(file => $file);
	if(scalar(keys(%{$ref_file->{index}})) > 0){
		foreach my $id (keys(%{$ref_file->{index}})) {
      my @t = split('\|',$id);
      my $acc=$t[3];
			# my ($gi) = ($id =~/gi\|([^\|]*)\|/);
			$self->{_map}->{$tax_id}->{$acc}->{description} = $ref_file->{index}->{$id}->{description};
			my $data = $ref_file->retrieveFastaSequence($id);
			$self->{_reference_length}->{$id} = length($data->{$id});
		}
	}
}


sub _print_mapping_file {
	my ($self,$out_file)=@_;
	open(OUT,">$out_file");
	print OUT '#' . join("\t",@{$self->{_mapping_file_headers}}) . "\n";
	foreach my $tax_id (keys(%{$self->{_map}})){
		foreach my $key (keys($self->{_map}->{$tax_id})){
			if($key =~ /taxonomy|organism/){
				next;
			}
			else{
				if(!defined($self->{_map}->{$tax_id}->{$key}->{read_map_file}) || !defined($self->{_map}->{$tax_id}->{$key}->{mat_file})){
					next;
				}
				print OUT $key . "\t" . $tax_id . "\t" ;
        print $key . "\t" . $tax_id . "\n" ;
        print Dumper $self->{_map}->{$tax_id}->{$key};
        print $self->{_map}->{$tax_id}->{$key}->{fasta_file} . "\n";
				if(defined($self->{_map}->{$tax_id}->{$key}->{fasta_file})){
					my $out = $mapperOptions{'outputDirResultsData'} . '/' . basename($self->{_map}->{$tax_id}->{$key}->{fasta_file});
					$logger->debug('Moving ' . $self->{_map}->{$tax_id}->{$key}->{fasta_file} . ' in ' . $out);
					move($self->{_map}->{$tax_id}->{$key}->{fasta_file},$out);
					print OUT '"' . 'data/' . basename($self->{_map}->{$tax_id}->{$key}->{fasta_file}) . '"' . "\t" ;

				}
				else{
					print OUT '.' . "\t";
				}
				if(defined($self->{_map}->{$tax_id}->{$key}->{read_map_file})){
					my $out = $mapperOptions{'outputDirResultsData'} . '/' . basename($self->{_map}->{$tax_id}->{$key}->{read_map_file});
					$logger->debug('Moving ' . $self->{_map}->{$tax_id}->{$key}->{read_map_file} . ' in ' . $out);
					move($self->{_map}->{$tax_id}->{$key}->{read_map_file},$out);
					print OUT '"' . 'data/' . basename($self->{_map}->{$tax_id}->{$key}->{read_map_file}) . '"' . "\t" ;
				}
				else{
					print OUT '.' . "\t";
				}
				if(defined($self->{_map}->{$tax_id}->{$key}->{mat_file})){
					my $out = $mapperOptions{'outputDirResultsData'} . '/' . basename($self->{_map}->{$tax_id}->{$key}->{mat_file});
					$logger->debug('Moving ' . $self->{_map}->{$tax_id}->{$key}->{mat_file} . ' in ' . $out);
					move($self->{_map}->{$tax_id}->{$key}->{mat_file},$out);
					print OUT '"' . 'data/' . basename($self->{_map}->{$tax_id}->{$key}->{mat_file}) . '"' . "\t" ;
				}
				else{
					print OUT '.' . "\t";
				}
				if(defined($self->{_map}->{$tax_id}->{organism})){
					print OUT '"' . $self->{_map}->{$tax_id}->{organism} . '"' . "\t" ;
				}
				else{
					print OUT '.' . "\t";
				}
				if(defined($self->{_map}->{$tax_id}->{taxonomy})){
					print OUT '"' . $self->{_map}->{$tax_id}->{taxonomy} . '"' . "\t" ;
				}
				else{
					print OUT '.';
				}
				if(defined($self->{_map}->{$tax_id}->{$key}->{description})){
					print OUT '"' . $self->{_map}->{$tax_id}->{$key}->{description} . '"';
				}
				else{
					print OUT '.';
				}

			}
			print OUT "\n";
		}
	}
}


sub _print_mapping_matrix {
	my ($self,$matches,$path,$tax_id,$orga) = @_;
	foreach my $m (keys(%{$matches})){

    my @t = split('\|',$m);
    my $acc=$t[3];

		my $out_file = $path . '/' . $tax_id . '_' . $orga . '_' . $acc . '_reads' . '.mat';
		open(OUT,">$out_file");
		print OUT '#x;coverage' . "\n";
		my @last_line;
		if(! defined($self->{_reference_length}->{$m})){
			$logger->error('Reference length for ' . $m . ' not found.');
			exit;
		}
		for(my $i=1;$i<$self->{_reference_length}->{$m};$i++){
			my @line;
			if(defined($matches->{$m}->{$i})){
				$line[0] = $i ;
				$line[1] = $matches->{$m}->{$i};
			}
			else{
				$line[0] = $i ;
				$line[1] = 0;
			}
			my $l = join(';',@line[1..$#line]);
			my $ll = join(';',@last_line[1..$#last_line]);
			if($l ne $ll){
				print OUT join(';',@last_line) . "\n";
				print OUT join(';',@line) . "\n";
			}
			@last_line = @line;
		}
		$self->{_map}->{$tax_id}->{$acc}->{read_map_file} = $out_file;
		close OUT;
	}
}


sub _read_mapping_coords {
	my ($self,$file) = @_;
	open(COORDS,$file) || $logger->logdie('File not found: ' . $file);
	my $coverage;
	while(my $line =<COORDS>){
		chomp $line;
		my @line = split("\t",$line);
		for(my $i=$line[8];$i<$line[9];$i++){
			$coverage->{$line[5]}->{$i}++;
		}
	}
	close COORDS;
	return $coverage;
}


sub _googleChartIdentityFromMummerDirectory {
	my ($self, $aln_file, $windowsSize, $hash_info, $tax_id) = @_;

	my $googleChartObject = Tools::GoogleChart::Chart ->new('line_chart');
	$googleChartObject->series(0)->lineWidth(20);
	$googleChartObject->hAxis->viewWindowMode('maximized');
	$googleChartObject->divHeight('800');

	if(! defined $windowsSize || ! $windowsSize >= 1){
		$logger->warn("Invalid windows size : $windowsSize. Set to 1");
		$windowsSize = 1;
	}
	my @chartFiles;
	my @matrix_files;
	$self->{_alignments} = _read_align_file($self,$aln_file);
	if(defined($self->{_alignments})){
		_get_identity_matrix($self,$self->{_alignments},$mapperOptions{'windowsSize'});
		_print_matrix($self,$self->{_alignments}, $hash_info, $tax_id);
	}
	else{
		return 0;
	}
}


sub _print_matrix {
	my ($self,$align_hash,$hash_info,$tax_id) = @_;
	my $file_list;
	foreach my $ref (keys(%{$align_hash})){

    my @t = split('\|',$ref);
    my $acc=$t[3];

		my $file_name = $hash_info->{mummerDirectory} . '/' . $tax_id . '_' . $hash_info->{organism} . '_' . $acc . '_scaffolds' . '.mat';
		$self->{_map}->{$tax_id}->{$acc}->{mat_file} = $file_name;
		open(MATRIX,">$file_name");
		my @sorted_qry = keys(%{$align_hash->{$ref}});
		my $col_names = _get_col_names($self,\@sorted_qry,$align_hash->{$ref});
		print MATRIX $ref . ";" . join(";",@{$col_names}) . "\n";
		my @last_line;
		if(!defined($self->{_reference_length}->{$ref})){
			$logger->error('Ref length for ' . $ref . ' not found.' );
			exit;
		}
		for(my $i=0;$i<$self->{_reference_length}->{$ref};$i++){
			my $void=1;
			my @line;
			push(@line,$i);
			foreach my $qry (@sorted_qry){
				my @sorted_hsp = sort(keys(%{$align_hash->{$ref}->{$qry}}));
				foreach my $hsp (@sorted_hsp){
					if(!defined($align_hash->{$ref}->{$qry}->{$hsp}->{reference_based_matrix}->[$i])){
						print $ref . "\n" . $qry . "\n" . $hsp . "\n";
						# print Dumper $align_hash->{$ref}->{$qry}->{$hsp}->{reference_based_matrix};
						print scalar(@{$align_hash->{$ref}->{$qry}->{$hsp}->{reference_based_matrix}}) . "\n" . $i . "\n";
						exit;
					}
					push(@line,$align_hash->{$ref}->{$qry}->{$hsp}->{reference_based_matrix}->[$i]);
					if($align_hash->{$ref}->{$qry}->{$hsp}->{reference_based_matrix}->[$i] ne ''){
						$void=0;
					}
				}
			}
			if($i==0){
				print MATRIX join(';',@line) . "\n";
			}
			if($void==0 || $i == $self->{_reference_length}->{$ref}-1){
				my $line = join(';',@line[1..$#line]);
				my $last_line = join(';',@last_line[1..$#last_line]);
				if($line ne $last_line){
					print MATRIX join(';',@last_line) . "\n";
					print MATRIX join(';',@line) . "\n";
				}
			}
			@last_line = @line;
		}
		close MATRIX;
	}
	return $file_list;
}


sub _get_col_names {
	my ($self,$qry,$align_hash) = @_;
	my @col_names;
	for(my $i=0;$i<scalar(@{$qry});$i++){
		my @hsp = sort(keys(%{$align_hash->{$qry->[$i]}}));
		for(my $j=0;$j<scalar(@hsp);$j++){
			push(@col_names,$qry->[$i] . '_' . $hsp[$j]);
		}
	}
	return \@col_names;
}


sub _get_identity_matrix {
	my ($self,$align_hash,$window_size) = @_;
	foreach my $ref (keys(%{$align_hash})){
		foreach my $qry (keys(%{$align_hash->{$ref}})){
			$logger->trace($ref . ' ' . $qry);
			foreach my $hsp (keys(%{$align_hash->{$ref}->{$qry}})){
				$align_hash->{$ref}->{$qry}->{$hsp}->{identity_matrix} = _compute_identiy_v2($align_hash->{$ref}->{$qry}->{$hsp}->{ref},$align_hash->{$ref}->{$qry}->{$hsp}->{qry},$window_size);
				$align_hash->{$ref}->{$qry}->{$hsp}->{reference_based_matrix} = _get_reference_based_matrix_v2($self,$align_hash->{$ref}->{$qry}->{$hsp},$ref);
				delete $align_hash->{$ref}->{$qry}->{$hsp}->{identity_matrix};
			}
		}
	}
}


sub _get_reference_based_matrix_v2 {
	my ($self,$hash,$ref)=@_;
	my $matrix;
	if(!defined($self->{_reference_length}->{$ref})){
		$logger->error('Ref length for ' . $ref . ' not found.');
		$logger->debug(Dumper $self->{_reference_length});
		exit;
	}
	for(my $i=0;$i<$self->{_reference_length}->{$ref};$i++){
		if($i >= $hash->{r_s}-1 && $i<=$hash->{r_e}-1){
			foreach my $k (@{$hash->{identity_matrix}}){
				for(my $j=$k->{'s'};$j<=$k->{'e'};$j++){
					$matrix->[$i] = $k->{'i'};
					$i++;
				}
			}
			$i--;
		}
		else{
			$matrix->[$i] = '';
		}
	}
	return $matrix;
}


sub _compute_identiy_v2 {
	my ($ref,$qry,$window_size) = @_;
	my @ref_seq = split('',$ref);
	my @qry_seq = split('',$qry);
	my $matrix;
	for(my $i=0;$i<scalar(@ref_seq);$i+=$window_size+1){
		# print $i . "\n";
		my $range = $i + $window_size-1;
		if ($range > scalar(@ref_seq)){
			$range = scalar(@ref_seq)-1;
		}
		my @ref_chunk = @ref_seq[$i..$range];
		my @qry_chunk = @qry_seq[$i..$range];
		# print join('',@ref_chunk) . "\n";
		# print join('',@qry_chunk) . "\n";
		my $identic=0;
		# print scalar(@ref_chunk) . "\n";
		for(my $j=0;$j<scalar(@ref_chunk);$j++){
			if(!defined($ref_chunk[$j]) || !defined($qry_chunk[$j])){
				next;
			}
			if($ref_chunk[$j] eq $qry_chunk[$j]){
				$identic++;
			}
		}
		my $identity = $identic / scalar(@ref_chunk);
		# print $identic . ' / ' . scalar(@ref_chunk) . ' : ' . $identity . "\n";
		if(defined($mummerOptions{'p'})){
			if($i==0){
				push(@{$matrix},{'s' => 1, 'e' => ($range+1)*3, 'i' => sprintf("%.2f", $identity)});
			}
			else{
				push(@{$matrix},{'s' => ($i-1)*3, 'e' => ($range+1)*3, 'i' => sprintf("%.2f", $identity)});
			}
		}
		else{
			push(@{$matrix},{'s' => $i, 'e' => $range, 'i' => sprintf("%.2f", $identity)});
		}

	}
	return $matrix;
}


sub _read_align_file {
	my ($self, $file)=@_;
	my $hash;
	open(ALIGN,$file) || $logger->logdie($file . ' not found.');
	my $void=0;
	my $r_s='';
	my $r_e='';
	my $q_s='';
	my $q_e='';
	my $ref='';
	my $query='';
	while(my $line =<ALIGN>){
		chomp $line;
		if($line =~ /^--\sAlignments\sbetween\s(.*)\sand\s(.*)$/){
			$ref = $1;
			$query = $2;
		}
		if($line =~ /^--\sBEGIN\salignment\s\[\s[+-]\d\s(\d+)\s-\s(\d+)\s\|\s[+-]\d\s(\d+)\s-\s(\d+)\s\]/){
			$r_s = $1;
			$r_e = $2;
			$q_s = $3;
			$q_e = $4;
			if($r_s > $r_e){
				my $tmp = $r_s;
				$r_s = $r_e;
				$r_e = $tmp;
			}
			if($q_s > $q_e){
				my $tmp = $q_s;
				$q_s = $q_e;
				$q_e = $tmp;
			}
			$hash->{$ref}->{$query}->{$r_s . '_' . $r_e . '_' . $q_s . '_' . $q_e}->{r_s} = $r_s;
			$hash->{$ref}->{$query}->{$r_s . '_' . $r_e . '_' . $q_s . '_' . $q_e}->{r_e} = $r_e;
			$hash->{$ref}->{$query}->{$r_s . '_' . $r_e . '_' . $q_s . '_' . $q_e}->{q_s} = $q_s;
			$hash->{$ref}->{$query}->{$r_s . '_' . $r_e . '_' . $q_s . '_' . $q_e}->{q_e} = $q_e;
		}
		if($line =~ /^$/){
			$void=1;
		}
		if($line =~ /^\d+\s+(.*)$/){
			my $seq = $1;
			if($void==1){
				$hash->{$ref}->{$query}->{$r_s . '_' . $r_e . '_' . $q_s . '_' . $q_e}->{ref} .= $seq;
			}
			else{
				$hash->{$ref}->{$query}->{$r_s . '_' . $r_e . '_' . $q_s . '_' . $q_e}->{qry} .= $seq;
			}
			$void=0;
		}
		# print $ref . "\n";
	}
	close ALIGN;
	return $hash;
}


sub mummerMapping {
  my ($self, $reference, $query, $mummerOptions, $outputDirectory) = @_;
	my $currentDirectory = cwd();
	my %options;

	if(defined $mummerOptions && ref $mummerOptions eq 'HASH'){
		%options = %$mummerOptions;
	}
	$options{'r'} = $reference;
	$options{'q'} = $query;
  chdir $outputDirectory;

  my $cmd = 'mummerGraph.pl';

  foreach my $option (keys %options){
		$cmd .= ' -' . $option . ' ' . $options{$option};
	}

	$logger->debug($cmd);
  `$cmd`;
  chdir $currentDirectory;
}


sub _checkReferenceSequence {
	my ($self, $file) = @_;
	$logger->info('Checking references...');
	my $ref_file = Tools::Fasta->new(file => $file);
	# my $idx = indexFastaFile($file);
	my $info_ref;
	if(scalar(keys(%{$ref_file->{index}})) > 1){
		foreach my $id (keys(%{$ref_file->{index}})) {
			if($ref_file->{index}->{$id}->{description} =~ /partial/){
				$logger->trace($id . ' partial ' . $ref_file->{index}->{$id}->{description});
				next;
			}
			elsif($ref_file->{index}->{$id}->{description} =~ /complete [genome|sequence]/){
				$logger->trace('complete genome ' . $id);
				$info_ref->{comp_gen}->{$id}=1;
			}
			elsif($ref_file->{index}->{$id}->{description} =~ /(.*)\ssegment\s(\S+)/){
				my $orga = $1;
				my $segment = $2;
				if($ref_file->{index}->{$id}->{description} =~ /complete sequence/){
					$logger->trace('multipartite ' . $orga . ' ' . $segment);
					$info_ref->{multipartite}->{$orga}->{$segment} = $id;
				}
			}
			elsif($ref_file->{index}->{$id}->{description} =~ /polyprotein/ && $ref_file->{index}->{$id}->{description} =~ /complete/){
				$logger->trace('polyprotein ' . $id);
				$info_ref->{polyprotein}->{$id}=1;
			}
			else{
				$logger->trace('Nothing revelant found for ' . $id . ' ' . $ref_file->{index}->{$id}->{description});
			}
		}
	}
	else{
		return 1;
	}
	my $ref_seq='';
	if(defined($info_ref->{comp_gen})){
		my @array = keys(%{$info_ref->{comp_gen}});
		if(scalar(@array) > 1) {
			my $id_maxLength = _seekLongestSeq($self,$info_ref->{comp_gen}, $ref_file);
			$logger->info('Reference sequence: ' . $id_maxLength);
			$ref_seq = $ref_file->retrieveFastaBlock($id_maxLength);
		}
		else{
			$logger->info('Reference sequence: ' . $array[0]);
			$ref_seq = $ref_file->retrieveFastaBlock($array[0]);
		}
	}
	elsif(defined($info_ref->{multipartite})){
		if(scalar(keys(%{$info_ref->{multipartite}})) > 1){
			my $max_orga='';
			my $max=0;
			foreach my $orga (keys(%{$info_ref->{multipartite}})) {
				if(scalar(keys(%{$info_ref->{multipartite}->{$orga}})) > $max){
					$max = scalar(keys(%{$info_ref->{multipartite}->{$orga}}));
					$max_orga = $orga;
				}
			}
			if($max==1){
				my $segments;
				foreach my $orga (keys(%{$info_ref->{multipartite}})) {
					foreach my $seg (keys(%{$info_ref->{multipartite}->{$orga}})) {
						$segments->{$seg}->{$info_ref->{multipartite}->{$orga}->{$seg}} = 1;
					}
				}
				foreach my $seg (keys($segments)){
					my $id_maxLength = _seekLongestSeq($self,$segments->{$seg}, $ref_file);
					$ref_seq .= $ref_file->retrieveFastaBlock($id_maxLength);
				}
			}
			else{
				foreach my $seg (keys(%{$info_ref->{multipartite}->{$max_orga}})) {
					$logger->info('multipartite sequence: ' . $max_orga . ' ' . $seg);
					$ref_seq .= $ref_file->retrieveFastaBlock($info_ref->{multipartite}->{$max_orga}->{$seg});
				}
			}
		}
		else{
			foreach my $orga (keys(%{$info_ref->{multipartite}})) {
				foreach my $seg (keys(%{$info_ref->{multipartite}->{$orga}})) {
					$logger->info('multipartite sequence: ' . $orga . ' ' . $seg);
					$ref_seq .= $ref_file->retrieveFastaBlock($info_ref->{multipartite}->{$orga}->{$seg});
				}
			}
		}
	}
	elsif(defined($info_ref->{polyprotein})){
		my @array = keys(%{$info_ref->{polyprotein}});
		if(scalar(keys(%{$info_ref->{polyprotein}})) > 1) {
			my $id_maxLength = _seekLongestSeq($self,$info_ref->{polyprotein}, $ref_file);
			$logger->info('polyprotein sequence: ' . $id_maxLength);
			$ref_seq = $ref_file->retrieveFastaBlock($id_maxLength);
		}
		else{
			$logger->info('polyprotein sequence: ' . $array[0]);
			$ref_seq = $ref_file->retrieveFastaBlock($array[0]);
		}
	}
	else{
		$logger->warn('no reference were found.')
	}
	if($ref_seq ne ''){
		open(FILE,">$file");
		print FILE $ref_seq;
		close FILE;
		return 1;
	}
	else{
		return 0;
	}
}


sub _seekLongestSeq {
	my ($self, $h, $fasta_tool) = @_;
	my $maxLength=0;
	my $id_maxLength='';
        my $seq_hash=$fasta_tool->retrieveFastaSequence([keys(%{$h})]);
	foreach my $id (keys(%{$h})) {
		if(length($seq_hash->{$id}) > $maxLength){
			$maxLength = length($seq_hash->{$id});
			$id_maxLength= $id;
		}
	}
	return $id_maxLength;
}


sub fastaFromAccList {
  my ($self, $database, $outputFile, $accListFile) = @_;
  my $cmd = 'fastacmd -d "' . $database . '" -D 1 -o ' . '"' . $outputFile . '"';
  if ($accListFile){
		$cmd = 'fastacmd -d "' . $database . '" -i "' . $accListFile . '" -o "' . $outputFile . '"' ;
	}
	$cmd .= ' 2> /dev/null';
  $logger->debug($cmd);
	system($cmd);
	if(_checkReferenceSequence($self,$outputFile)){
		return 1;
	}
	else{
		return 0;
	}
}


sub accListFromTaxid {
  my ($self, $taxid, $outputFile, $mapperOptions) = @_;
  my $dbh = DBI->connect("dbi:SQLite:dbname=" . $mapperOptions->{'taxonomyDatabase'},"","") or $logger->logdie("Error connecting to database" . $mapperOptions->{'taxonomyDatabase'});
  $logger->info("Retrieving accession from taxid...");

  my $cmd = 'SELECT accession FROM nucl_accession2taxid WHERE taxid=' . $taxid;
  $logger->debug($cmd);
  my $query = $dbh->prepare($cmd) or $logger->logdie("Error preparing query for retrieving accession from taxid $taxid");
  $query->execute();

  my $res ;
  $res = $query->fetchall_arrayref();

	if(scalar(@{$res}) == 0){
		return 0;
	}
	else{
  	$self->printgiFileList($outputFile,$res);
  }
  $dbh->disconnect;
}


sub printgiFileList {
	my ($self, $outputFile, $giList) = @_;
	open (GI, '>'. $outputFile) || $logger->logdie("Cannot create the gi list file " . $outputFile . " : $!");
    for(my $i=0; $i<= $#{$giList}; $i++){
	    foreach my $gi (@{$giList->[$i]}){
	    	print GI $gi . "\n";
	    }
	}
    close GI;
}


sub readExtentedBlast {
  my ($self,$blastFile, $filter) = @_;
  my %taxidInfos;
	my $file = Tools::Taxonomy->readCSVextended($blastFile, "\t");
	my $taxid_infos;
  foreach my $hit (@{$file}){
		if($hit->{taxonomy} =~ /$filter/){
			$hit->{organism} =~ s/[\.\-\s\/]+/_/g;
			$taxid_infos->{$hit->{tax_id}}->{organism} = $hit->{organism};
			$taxid_infos->{$hit->{tax_id}}->{taxonomy} = $hit->{taxonomy};
			$taxid_infos->{$hit->{tax_id}}->{description} = $hit->{description};
		}
	}
	return $taxid_infos;
}


sub processOptions {
	my ($self, $mapperOptions, $mummerOptions) = @_;
	if(! defined $mapperOptions->{'windowsSize'} || $mapperOptions->{'windowsSize'} !~ /[0-9]*/ || $mapperOptions->{'windowsSize'} < 1 ){
		$logger->error('Windows size must be an integer value >= 1');
		&help;
	}
	foreach my $option ( ('taxonomyDatabase', 'queries', 'blast', 'reads_file') ){
		if(! -e $mapperOptions->{$option} ){
			$logger->error($option . ' : file ' . $mapperOptions->{$option} . ' doesn\'t exist');
			&help;
		}
	}
	foreach my $option ( ('outputDir', 'taxonomyDatabase', 'queries', 'blast', 'taxonomyDatabase', 'blastDatabase', 'reads_file') ){
		$mapperOptions->{$option} = abs_path($mapperOptions->{$option});
		if($option eq 'queries'){
			my @a = split('_',$mapperOptions->{$option});
			$self->{_sample_id} = basename($a[0]);
		}
	}
	$mapperOptions->{'outputDirAcc'} = $mapperOptions->{'outputDir'} . '/acc/';
	$mapperOptions->{'outputDirFasta'} = $mapperOptions->{'outputDir'} . '/fasta/';
	$mapperOptions->{'outputDirMummer'} = $mapperOptions->{'outputDir'} . '/mummer/';
	$mapperOptions->{'outputDirResults'} = $mapperOptions->{'outputDir'} . '/results/';
	$mapperOptions->{'outputDirResultsData'} = $mapperOptions->{'outputDirResults'} . 'data';
	$self->createFiles([$mapperOptions->{'outputDirResultsData'},$mapperOptions->{'outputDirAcc'}, $mapperOptions->{'outputDirFasta'}, $mapperOptions->{'outputDirMummer'}, $mapperOptions->{'outputDirResults'}], 1);
}


sub createFiles {
  my ($self, $files, $isDirectory) = @_;
  if(! ref $files){
    $files = [$files];
  }
  foreach my $file (@$files){
    if(! -e $file){
      if ($isDirectory) {
        make_path ($file) or $logger->logdie ('Cannot create the directory ' . $file . " : $!");
      }
      else {
        open(FILE,">$file") or $logger->logdie('Cannot create the file ' . $file . " : $!");
        close(FILE);
      }
    }
    else{
    	$logger->warn($file . ' already exist');
    }
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
# PURPOSE:    Generate a FASTA reference for each taxid present in the blast file
#             and remap the queries on their references with mummer.
#

USAGE: perl $prog -q fasta_file -b extented_blast -r blast_database -t taxonomy_database [OPTIONS]

          ### OPTIONS ###
          -q|qry       [FASTA FILE]
          -b|blast     [EXTENDED BLAST FILE]
          -r|ref       [BLAST DATABASE]
          -t|taxonomy  [TAXONOMY DATABASE]
          -o|output    [OUTPUT DIRECTORY]
          -w|windows   [INTEGER] Windows size used to calculate de identity % of the alignment over the reference.
          -f|filter    [STRING] Keep only hit whose taxonomy match the filter string
          -v|verbosity [INTEGER]
          -help|h      Print this help and exit.


          ### ALIGNMENT OPTIONS ###
          -p|promer      Amino acid alignments between two mutli-FASTA DNA input.
          -m|mode        mum, mumreference, maxmatch. (Default mum)
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
          -l|length      Set the minumum length of an HSP. (default 0)
          -i|identity    Set the minimum identity of an HSP. (default 0)
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
          -f|format      graph, png, btab, tab, coords, align, fasta, fig. (default graph).

EOF
exit(1);
}
