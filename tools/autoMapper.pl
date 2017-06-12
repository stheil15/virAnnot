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
use File::Basename;
use Logger::Logger;

my $VERSION = '2.0';
my $lastmodif = '2017-05-30';

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
	$self->{_mapping_file_headers} = ['tax_id','fasta_file','svg_file','coords_file','organism','taxonomy','coverage'];
	$self->processOptions(\%mapperOptions,\%mummerOptions);
	my $taxid_infos = $self->readExtentedBlast($mapperOptions{blast}, $mapperOptions{filter});
	foreach my $tax_id (keys(%{$taxid_infos})){
		$self->{_map}->{$tax_id}->{taxonomy} = $taxid_infos->{$tax_id}->{taxonomy};
		$self->{_map}->{$tax_id}->{organism} = $taxid_infos->{$tax_id}->{organism};

		$taxid_infos->{$tax_id}->{acc_file} = $mapperOptions{outputDirAcc} . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '.csv';

		if(! $self->accListFromTaxid($tax_id, $taxid_infos->{$tax_id}->{acc_file}, \%mapperOptions)) {
			$logger->warn('no reference found, skipping...');
			next;
		}

		$taxid_infos->{$tax_id}->{fasta_ref} = $mapperOptions{outputDirFasta} . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '.fasta';
		if(! $self->fastaFromAccList($mapperOptions{blastDatabase}, $taxid_infos->{$tax_id}->{fasta_ref}, $taxid_infos->{$tax_id}->{acc_file})){
			next;
		}
		$taxid_infos->{$tax_id}->{mummerDirectory} = $mapperOptions{outputDirMummer} . $tax_id . '_' . $taxid_infos->{$tax_id}->{organism};

    # Extract sequences from query file to align versus the selected reference.
    _get_query_sequences($self, $mapperOptions{queries},$taxid_infos->{$tax_id}->{query_seq},$taxid_infos->{$tax_id}->{mummerDirectory});

		$self->createFiles([$taxid_infos->{$tax_id}->{mummerDirectory}], 1);
		_get_sequence_info($self,$tax_id,$taxid_infos->{$tax_id}->{fasta_ref});
		$mummerOptions{format} = 'fasta';
		$self->mummerMapping($taxid_infos->{$tax_id}->{fasta_ref}, $taxid_infos->{$tax_id}->{mummerDirectory} . "/queries.fasta", \%mummerOptions, $taxid_infos->{$tax_id}->{mummerDirectory});


    # Storing file paths of alignments created by mummerGraph.pl
    foreach my $file (glob ($taxid_infos->{$tax_id}->{mummerDirectory} . '/fasta/*.fasta')){
			if(-z $file){
				next;
			}
      $logger->debug('MummerGraph fasta alignment found ' . $file);
      my $acc = basename($file);
      $acc =~ s/\.fasta//;
      if(defined($self->{_map}->{$tax_id}->{fasta_file}->{file})){
        print $file . "\n";
      }
			push(@{$self->{_map}->{$tax_id}->{fasta_file}->{file}},$file);
			# $self->{_map}->{$tax_id}->{fasta_file}->{file} = $file;
			# $self->{_map}->{$tax_id}->{fasta_file}->{desc} = $tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '_' . $acc  ;
			push(@{$self->{_map}->{$tax_id}->{fasta_file}->{desc}},$tax_id . '_' . $taxid_infos->{$tax_id}->{organism} . '_' . $acc);
		}

		$mummerOptions{format} = 'tab';
		$self->mummerMapping($taxid_infos->{$tax_id}->{fasta_ref}, $taxid_infos->{$tax_id}->{mummerDirectory} . "/queries.fasta", \%mummerOptions, $taxid_infos->{$tax_id}->{mummerDirectory});
    # exit;
		foreach my $aln_file (glob ($taxid_infos->{$tax_id}->{mummerDirectory} . '/*.coords')){
      my $cmd = 'wc -l < ' . $aln_file;
      my $res = `$cmd`;
			if($res == 4){
        $logger->debug('No alignments to plot... skip');
				next;
			}
      $self->{_map}->{$tax_id}->{cov_info} = _compute_genome_coverage($self,$aln_file);
      print Dumper $self->{_map}->{$tax_id}->{cov_info};
      my $under_threshold=0;
      foreach my $ref (keys(%{$self->{_map}->{$tax_id}->{cov_info}})){
        if($self->{_map}->{$tax_id}->{cov_info}->{$ref}->{cov} < 1 && $self->{_map}->{$tax_id}->{cov_info}->{$ref}->{depth} <= 2){
          $under_threshold++;
        }
      }
      if($under_threshold == scalar(keys(%{$self->{_map}->{$tax_id}->{cov_info}}))){
        $logger->debug('No reference coverage over threshold.');
        $self->{_map}->{$tax_id}->{under_threshold}=1;
      }
      _launchMummerToPip($self,$aln_file,$taxid_infos->{$tax_id},$tax_id);
      $self->{_map}->{$tax_id}->{svg_file} = $aln_file . '.svg';
      $self->{_map}->{$tax_id}->{aln_file} = $aln_file;
		}
	}
	if(scalar(keys(%{$self->{_map}})) != 0){
		my $map_file_name = $mapperOptions{'outputDirResults'} . '/' . 'map.txt';
		_print_mapping_file($self,$map_file_name,$mapperOptions{'outputDirResults'});
    _create_html($self,$map_file_name,$mapperOptions{'outputDirResults'});
	}
}


sub _create_html {
  my ($self,$map_file,$out_dir)=@_;
  my $cmd = 'autoMapper_html.py -m ' . $map_file . ' -t ' . $mapperOptions{'outputDir'} . ' -o ' . $out_dir ;
  $logger->debug($cmd);
  `$cmd`;
}


sub _compute_genome_coverage {
  my ($self,$aln_file) = @_;
  print $aln_file . "\n";
  open(FILE,$aln_file);
  my $ref_arrays;
  my $queries={};
  while(<FILE>){
    chomp;
    if(! /^\d+/){next;}
    else{
      my @line = split('\s+',$_);
      if(defined($queries->{$line[13]}->{$line[14]})){
        next;
      }
      $queries->{$line[13]}->{$line[14]}++;
      if(! $ref_arrays->{$line[13]}){
        @{$ref_arrays->{$line[13]}} = ('0') x $line[9];
      }
      for(my $i=$line[0];$i<=$line[1];$i++){
        $ref_arrays->{$line[13]}->[$i]++;
      }
    }
  }
  my $info;
  foreach my $ref (keys(%{$ref_arrays})) {
    my $nucl_covered=0;
    my $max_depth=0;
    for(my $i=0; $i<scalar(@{$ref_arrays->{$ref}});$i++){
      if($ref_arrays->{$ref}->[$i] != 0){
        $nucl_covered++;
        if($ref_arrays->{$ref}->[$i] > $max_depth){
          $max_depth = $ref_arrays->{$ref}->[$i];
        }
      }
    }
    my $cov = $nucl_covered / scalar(@{$ref_arrays->{$ref}});
    $info->{$ref}->{cov} = sprintf("%.2f",$cov*100);
    $info->{$ref}->{depth} = $max_depth;
  }
  return $info;
}


sub _get_query_sequences {
  my ($self,$fasta,$qry_lst,$dir)=@_;
  $logger->debug('Extracting query seq ' . Dumper $qry_lst);
  if(! -e $dir){
    mkdir $dir;
  }
  my $fasta_tool = Tools::Fasta->new(file=>$fasta);
  my $out_f = $dir  . '/queries.fasta';
  open(OUT,">$out_f") || die 'cannot create file ' . $out_f;
  foreach my $id (@{$qry_lst}){
    print OUT $fasta_tool->retrieveFastaBlock($id) . "\n";
  }
  close OUT
}


sub _launch_coverage_plot {
  my ($self,$file,$orga) = @_;
  my $cmd = 'coverage_plot.py -b ' . $file . ' -t ' . $orga . ' -o ' . $file ;
  `$cmd`;
}


sub _launch_bowtie {
  my ($self, $ref, $reads, $dir) = @_;
  my $cmd = 'cd ' . $dir . "\n";
  $cmd .= 'bowtie2-build ' . $ref . ' ' . $ref . "\n";
  $cmd .= 'bowtie2 --sensitive-local --no-unal -p 5 -f -U ' . $reads . ' -x ' . $ref;
  $cmd .= ' | samtools view -F 4 -bS - | samtools sort - ' . $ref . "\n";
  $cmd .= 'samtools index ' . $ref . '.bam' . "\n";
  `$cmd`;
}


sub _launchMummerToPip {
  my ($self, $aln_file, $hash_info, $tax_id) = @_;
  my $cmd = 'mummerToPip.py -p "" -i ' . $aln_file . ' -t ' . $hash_info->{organism} . ' -o ' . $aln_file ;
  $logger->debug($cmd);
  `$cmd`;
}


sub _get_sequence_info {
	my ($self,$tax_id,$file)=@_;
  $logger->debug('Getting sequence informations...');
	my $ref_file = Tools::Fasta->new(file => $file);
	if(scalar(keys(%{$ref_file->{index}})) > 0){
		foreach my $id (keys(%{$ref_file->{index}})) {
			$self->{_map}->{$tax_id}->{$id}->{description} = $ref_file->{index}->{$id}->{description};
			my $data = $ref_file->retrieveFastaSequence($id);
			$self->{_reference_length}->{$id} = length($data->{$id});
		}
	}
  $logger->debug('Done.');
}


sub _print_mapping_file {
	my ($self,$out_file)=@_;
	open(OUT,">$out_file");
	print OUT '#' . join("\t",@{$self->{_mapping_file_headers}}) . "\n";
	foreach my $tax_id (keys(%{$self->{_map}})){
    if((! defined($self->{_map}->{$tax_id}->{fasta_file}) && ! defined($self->{_map}->{$tax_id}->{svg_file})) || defined($self->{_map}->{$tax_id}->{under_threshold})){
      next;
    }
		print OUT $tax_id . "\t" ;
		if(defined($self->{_map}->{$tax_id}->{fasta_file})){
      print OUT '"';
      for(my $i=0;$i<scalar(@{$self->{_map}->{$tax_id}->{fasta_file}->{file}});$i++){
        my $acc = basename($self->{_map}->{$tax_id}->{fasta_file}->{file}->[$i]);
        $acc =~ s/\.fasta//;

        my $out = $mapperOptions{'outputDirResultsData'} . '/' . $self->{_map}->{$tax_id}->{fasta_file}->{desc}->[$i] . '.fasta';
        $logger->debug('Moving ' . $self->{_map}->{$tax_id}->{fasta_file}->{file}->[$i] . ' in ' . $out);
        move($self->{_map}->{$tax_id}->{fasta_file}->{file}->[$i],$out);
        if($i>0){
          print OUT ','
        }
        print OUT 'data/' . basename($out) ;
      }
      print OUT '"' . "\t";
		}
		else{
			print OUT '.' . "\t";
		}
    if(defined($self->{_map}->{$tax_id}->{svg_file})){
      my $out = $mapperOptions{'outputDirResultsData'} . '/' . basename($self->{_map}->{$tax_id}->{svg_file});
      move($self->{_map}->{$tax_id}->{svg_file},$out);
      print OUT '"' . 'data/' . basename($out) . '"' . "\t" ;
    }
    else{
			print OUT '.' . "\t";
		}
    if(defined($self->{_map}->{$tax_id}->{aln_file})){
      my $out = $mapperOptions{'outputDirResultsData'} . '/' . basename($self->{_map}->{$tax_id}->{aln_file});
      move($self->{_map}->{$tax_id}->{aln_file},$out);
      print OUT '"' . 'data/' . basename($out) . '"' . "\t" ;
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
			print OUT '.'. "\t" ;
		}
    if(defined($self->{_map}->{$tax_id}->{cov_info})){
      print OUT '"' ;
      my @keys = keys(%{$self->{_map}->{$tax_id}->{cov_info}});
      for(my $i=0;$i<scalar(@keys);$i++){
        print OUT $keys[$i] . ';coverage:' . $self->{_map}->{$tax_id}->{cov_info}->{$keys[$i]}->{cov} . ';max_depth:' . $self->{_map}->{$tax_id}->{cov_info}->{$keys[$i]}->{depth};
        if($i <= (scalar(@keys)-1)){
          print OUT ',';
        }
      }
      print OUT '"' ;
    }
	  print OUT "\n";
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
  if($verbosity >= 3){
    $cmd .= ' -v 3 ';
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
      $logger->trace($ref_file->{index}->{$id}->{description});
			if($ref_file->{index}->{$id}->{description} =~ /partial/){
				# $logger->trace('partial **');
				next;
			}
			elsif($ref_file->{index}->{$id}->{description} =~ /complete genome/){
				$logger->trace('** complete genome **');
				$info_ref->{comp_gen}->{$id}=1;
			}
			elsif($ref_file->{index}->{$id}->{description} =~ /(.*)\ssegment\s(\S+)/){
				my $orga = $1;
				my $segment = $2;
				if($ref_file->{index}->{$id}->{description} =~ /complete/){
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
  elsif(scalar(keys(%{$ref_file->{index}})) == 1){
    return 1;
  }
	else{
		return 0;
	}
	my $ref_seq='';
  if(defined($info_ref->{multipartite})){
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
  elsif(defined($info_ref->{comp_gen})){
    my @array = keys(%{$info_ref->{comp_gen}});
    if(scalar(@array) > 1) {
      my $id_maxLength = _seekLongestSeq($self,$info_ref->{comp_gen}, $ref_file);
      $logger->info('Reference complete genome: ' . $id_maxLength);
      $ref_seq = $ref_file->retrieveFastaBlock($id_maxLength);
    }
    else{
      $logger->info('Reference complete genome: ' . $array[0]);
      $ref_seq = $ref_file->retrieveFastaBlock($array[0]);
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

  my $cmd = 'SELECT "accession.version" FROM nucl_accession2taxid WHERE taxid=' . $taxid;
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
			push(@{$taxid_infos->{$hit->{tax_id}}->{query_seq}},$hit->{query_id});
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
    if(defined($mapperOptions->{$option})){
      if(! -e $mapperOptions->{$option} ){
        $logger->error($option . ' : file ' . $mapperOptions->{$option} . ' doesn\'t exist');
        &help;
      }
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


          ### MUMMER ALIGNMENT OPTIONS ###
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

          ### MUMMER FILTERING OPTIONS ###
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
