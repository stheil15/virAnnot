#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use File::Which;
use File::Copy;
use Cwd;
use Pod::Usage;
use Logger::Logger;
use Tools::Taxonomy;
use Tools::Blast;
use DBI;
use XML::Writer;
use Data::Dumper;
use Math::Round;

my $VERSION = '2.0' ;
my $lastmodif = '2017-05-30' ;

my @input_files;
my @xml_files;
my @id_samples;
my @group_file;
my $output = '';
my $outdir = '';
my $merge;
my $mergeType=1;
my $help;
my $man;
my $verbosity = 1;
my $color = 'none';
my $data = 'none';
my $db = '/media/data/db/taxonomy/taxonomy.sqlite';
my $dbType = 'ncbi';
my $reduceTaxo;

GetOptions(
          "i|input:s"      => \@input_files,
          "x|xml:s"        => \@xml_files,
          "g|group:s"      => \@group_file,
          "id:s"           => \@id_samples,
          "o|output=s"     => \$output,
          "outdir=s"       => \$outdir,
          "m|merge"        => \$merge,
          "mt|mergeType:i" => \$mergeType,
          "db=s"           => \$db,
          "dt|dbtype=s"    => \$dbType,
          "c|color=s"      => \$color,
          "data=s"         => \$data,
          "h|help"         => \$help,
          "r|reduce"       => \$reduceTaxo,
          "v|verbosity=i"  => \$verbosity,
) ;


if($verbosity > 1){
  Logger::Logger->changeMode( $verbosity );
}


&main;


sub main {
  my $self={};
  bless $self;

  _set_options($self);

  if(! -e $self->{_krona_out_dir}){
    $logger->info('Creating ' . $self->{_krona_out_dir} . ' folder.');
    mkdir $self->{_krona_out_dir};
    chdir $self->{_krona_out_dir};
    my $source = dirname(which('blast2html.py')) . "/sequences.html";
    my $destination = cwd() . "/sequences.html";
    copy($source, $destination) or die "Copy failed: $!, $source to $destination";
  }
  else{
    chdir $self->{_krona_out_dir};
    if(! -e (cwd() . "/sequences.html")){
      my $source = dirname(which('blast2html.py')) . "/sequences.html";
    my $destination = cwd() . "/sequences.html";
    copy($source, $destination) or die "Copy failed: $!, $source to $destination";
    }
  }
  if(defined($self->{xml_files})){
    if(! -e $self->{_out_dir}){
      mkdir $self->{_out_dir};
    }
  }

  for(my $i=0;$i<=$#{$self->{filesList}};$i++){
    $logger->debug('Loading file ' . $self->{filesList}->[$i]);
    $self->{_hits}->[$i] = $self->{taxoTools}->readCSVextended($self->{filesList}->[$i],"\t");
    if(defined($self->{xml_files})){
      $self->{_qry_lst}->[$i] = _get_qry_list($self,$self->{_hits}->[$i]);
      $self->{_extract_sequences}->[$i] = _extract_sequences($self,$i);
      $self->{_html_files}->[$i] = _extract_blast($self,$i);
    }
  }

  my $tree = [];
  if(defined $merge && $mergeType == 2){
    $logger->info('Merging hits in a single analysis...');
    for(my $i=0;$i<=$#{$self->{filesList}};$i++){
      push(@{$self->{_merged_hits}}, @{$self->{_hits}->[$i]});
      delete $self->{_hits}->[$i];
    }
    $logger->info('Building taxonomy tree...');
    $tree->[0] = $self->{taxoTools}->buildFullTaxo( $self->{_merged_hits} );
    printXML($self, $tree, [0], $color, $data);
  }
  else{
    for(my $i=0;$i<=$#{$self->{filesList}};$i++){
      if(defined($self->{group_file}->[$i]) && $self->{group_file}->[$i] ne ''){
        $tree->[$i] = $self->{taxoTools}->buildFullTaxo( $self->{_hits}->[$i], 1 );
      }
      else{
        $tree->[$i] = $self->{taxoTools}->buildFullTaxo( $self->{_hits}->[$i] );
      }
    }
  }
  if(! defined($merge)){
    for(my $i=0;$i<=$#{$self->{filesList}};$i++){
      printXML($self, $tree, [$i], $color, $data);
    }
  }
  else{
    my $index_list;
    for(my $i=0;$i<=$#{$self->{filesList}};$i++){
      push(@{$index_list},$i);
    }
    printXML($self, $tree, $index_list, $color, $data);
  }
}

sub _get_qry_list {
  my ($self,$hits)=@_;
  my $h={};
  foreach my $hit (@{$hits}){
    if($hit->{'taxonomy'} =~ /Viruses/ || $hit->{'taxonomy'} =~ /Viroids/){
      $h->{$hit->{'tax_id'}}->{$hit->{query_id}}++;
      # $logger->info('hits' . $hits);
    }
  }
  $logger->debug('Selected ' . scalar(%{$h}) . ' hits.');
  return $h;
}

# Extract sequence for each taxo
sub _extract_sequences {
  my ($self,$index)=@_;
  my $h={};
  foreach my $hit (@{$self->{_hits}->[$index]}){
    my $seq=$hit->{'sequence'};
    my $tax_id=$hit->{'tax_id'};
    $h->{$tax_id}->{$seq}++;
  }
  # # read $h
  # foreach my $key (keys %$h){
  #   $logger->info($key);
  #   foreach my $inner (keys %{$h->{$key}}){
  #     $logger->info($inner);
  #   }
  # }

  return $h;
}


sub getDeconvolutedHit {
  my($blastHit, $groupInfo) = @_;
  my @deconvolutedHit;
  my %newHit = %$blastHit;
  if(! $newHit{'nb_hit'}){
    $newHit{'nb_hit'} = 1;
  }
  if(! $newHit{'nb_reads'}){
    $newHit{'nb_reads'} = 1;
  }
  my $totalReads = $newHit{'nb_reads'};
  $logger->trace('Deconvoluting hit ' . $blastHit->{'query_id'} . '...');
  if(defined $groupInfo){
    $logger->trace('Searching group in ' . (scalar keys %$groupInfo) . ' keys ...');
    foreach my $group (keys %$groupInfo){
      if($blastHit->{'query_id'} =~ /^$group$/){
        $totalReads = 0;
        $logger->trace('Hit ' . $blastHit->{'query_id'} . ' matched with group ' .$group);
        foreach my $substitutionHit (keys %{$groupInfo->{$group}}){
          my %tempHit = %newHit;
          $tempHit{'query_id'} =~ s/$group/$substitutionHit/ee;
          $tempHit{'nb_reads'} = $groupInfo->{$group}{$substitutionHit};
          $totalReads += $tempHit{'nb_reads'};
          push (@deconvolutedHit, {%tempHit});
        }
      }
    }
  }
  if(! @deconvolutedHit){
    @deconvolutedHit = (\%newHit);
  }
  $logger->trace('Hit ' . $blastHit->{'query_id'} . ' deconvoluted in ' .(scalar @deconvolutedHit). ' hits');
  foreach my $hit (@deconvolutedHit){
    $hit->{'nb_hit'} = $hit->{'nb_hit'} * ( $hit->{'nb_reads'} / $totalReads );
  }
  return @deconvolutedHit;
}


# Index_list : array of file name e.g. MID-GENC03_nr_dmdx 
sub printXML {
  my ($self, $tree, $index_list, $color, $data) = @_;
  my $output = IO::File->new(">$self->{_krona_file_name}");
  my $writer = XML::Writer->new(OUTPUT => $output, DATA_INDENT => " ", DATA_MODE => 1);
  print $output '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">' . "\n";
  print $output '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">' . "\n";
  print $output ' <head>' . "\n";
  print $output '  <meta charset="utf-8"/>' . "\n";
  print $output '  <link rel="shortcut icon" href="http://krona.sourceforge.net/img/favicon.ico"/>' . "\n";
  print $output '  <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"http://krona.sourceforge.net\"."}</script>' . "\n";
  print $output '  <script src="http://krona.sourceforge.net/src/krona-2.0.js"></script>' . "\n";
  print $output ' </head>' . "\n";
  print $output ' <body>' . "\n";
  print $output '  <img id="hiddenImage" src="http://krona.sourceforge.net/img/hidden.png" style="display:none"/>' . "\n";
  print $output '  <img id="loadingImage" src="http://krona.sourceforge.net/img/loading.gif" style="display:none"/>' . "\n";
  print $output '  <img id="logo" src="http://krona.sourceforge.net/img/logo.png" style="display:none"/>' . "\n";
  print $output '  <noscript>Javascript must be enabled to view this page.</noscript>' . "\n";
  print $output '  <div style="display:none">' . "\n";
  $writer->startTag('krona', 'collapse' => "false", 'key' => "true");
  $writer->startTag('attributes', 'magnitude' => 'magnitude');
  $writer->startTag('attribute', 'display' => 'Total');
  $writer->characters("magnitude");
  $writer->endTag('attribute');
  if(defined($self->{_qry_lst})){
    $writer->startTag('attribute', 'display' => 'Blast results', 'hrefBase' => './');
    $writer->characters("blast");
    $writer->endTag('attribute');
  }
  if(defined($self->{_extract_sequences})){
    $writer->startTag('attribute', 'display' => 'Sequences', 'hrefBase' => './');
    $writer->characters("sequence");
    $writer->endTag('attribute');
  }
  if ($color eq 'taxid' ) {
    $writer->startTag('attribute', 'display' => 'taxId');
    $writer->characters("taxid");
    $writer->endTag('attribute');
    $writer->endTag('attributes');
    $writer->startTag('color', attribute => "taxid", valueStart => 0, valueEnd => 100000000, hueStart => 0, hueEnd => 100000000, default => 'true');
    $writer->endTag('color');
  }
  elsif( $color eq 'identity' ){
    $writer->startTag('attribute', 'display' => 'identity');
    $writer->characters("identity");
    $writer->endTag('attribute');
    $writer->endTag('attributes');
    $writer->startTag('color', attribute => "identity", valueStart => 0, valueEnd => 100, hueStart => 0, hueEnd => 100, default => 'true');
    $writer->endTag('color');
  }
  else{
    $writer->endTag('attributes');
  }


  $writer->startTag('datasets');
  if(defined $merge && $mergeType == 2){
    my $suffix = '';
    for(my $i=0;$i<=$#{$index_list};$i++){
      $suffix .= $self->{id_samples}->[$i] . ' ';
    }
    if($data eq 'contigs'){
      $writer->startTag('dataset');
      $writer->characters('contigs hits '.$suffix);
      $writer->endTag('dataset');
    }
    elsif($data eq 'reads'){
      $writer->startTag('dataset');
      $writer->characters('reads hits '.$suffix);
      $writer->endTag('dataset');
    }
    elsif($data eq 'both'){
      $writer->startTag('dataset');
      $writer->characters('contigs hits '.$suffix);
      $writer->endTag('dataset');
      $writer->startTag('dataset');
      $writer->characters('reads hits '.$suffix);
      $writer->endTag('dataset');
    }
    else{
      $writer->startTag('dataset');
      $writer->characters($suffix);
      $writer->endTag('dataset');
    }
  }
  else{
    for(my $i=0;$i<=$#{$index_list};$i++){
      my $suffix = $self->{id_samples}->[$i];
      if($data eq 'contigs'){
        $writer->startTag('dataset');
        $writer->characters('contigs hits '.$suffix);
        $writer->endTag('dataset');
      }
      elsif($data eq 'reads'){
        $writer->startTag('dataset');
        $writer->characters('reads hits '.$suffix);
        $writer->endTag('dataset');
      }
      elsif($data eq 'both'){
        $writer->startTag('dataset');
        $writer->characters('contigs hits '.$suffix);
        $writer->endTag('dataset');
        $writer->startTag('dataset');
        $writer->characters('reads hits '.$suffix);
        $writer->endTag('dataset');
      }
      else{
        $writer->startTag('dataset');
        $writer->characters($suffix);
        $writer->endTag('dataset');
      }
    }
  }
  $writer->endTag('datasets');


  $writer->startTag('node', 'name' => "all");
  $writer->startTag('magnitude');
  foreach my $i (@{$index_list}){
    my $total=0;
    my $total_r=0;
    foreach my $k (keys(%{$tree->[$i]})){
      if($data eq 'contigs'){
        $total += nearest(0.00001, $tree->[$i]->{$k}->{number});
      }
      elsif($data eq 'reads'){
        $total_r += $tree->[$i]->{$k}->{number_r};
      }
      elsif($data eq 'both'){
        $total += nearest(0.00001, $tree->[$i]->{$k}->{number});
        $total_r += $tree->[$i]->{$k}->{number_r};
      }
      else{
        $total += nearest(0.00001, $tree->[$i]->{$k}->{number});
      }
    }
    if($data eq 'contigs'){
      $writer->dataElement(val => $total);
    }
    elsif($data eq 'reads'){
      $writer->dataElement(val => $total_r);
    }
    elsif($data eq 'both'){
      $writer->dataElement(val => $total);
      $writer->dataElement(val => $total_r);
    }
    else{
      $writer->dataElement(val => $total);
    }
  }

  $writer->endTag('magnitude');
  XMLPrinter($self,$tree,$writer,$color,$data, $index_list);
  $writer->endTag('node');
  $writer->endTag('krona');
  print $output '</div></body></html>' . "\n";
  $writer->end();
}


# param index_list : list of taxonomic name
# param data : both | read | contig
sub XMLPrinter {
  my ($self,$tree,$writer,$color,$data, $index_list) = @_;
  my %keysList;
  foreach my $i (@{$index_list}){
    if(exists $tree->[$i]){
      foreach my $k (keys %{$tree->[$i]}){
        if($k !~/rank|number|taxid|identity|number_r|identity_r|sequence/){
          $keysList{$k} = $tree->[$i]->{$k}->{taxid};
        }
      }
    }
  }
  # $k correspond to the taxonomic name
  foreach my $k (keys %keysList){
    my $newHash;
    $writer->startTag('node', name => "$k");
    $writer->startTag('magnitude');
    foreach my $i (@{$index_list}){
      if(exists $tree->[$i]->{$k} && defined $tree->[$i]->{$k}){
        $newHash->[$i] = $tree->[$i]->{$k};
      }
      # $logger->info('data ' . $data);
      if($data eq 'contigs'){
        if(exists $tree->[$i]->{$k}->{number} && defined $tree->[$i]->{$k}->{number}){
          $writer->dataElement(val => nearest(0.00001, $tree->[$i]->{$k}->{number}));
        }
        else{
          $writer->dataElement(val => '');
        }
      }
      elsif($data eq 'reads'){
        if(exists $tree->[$i]->{$k}->{number_r} && defined $tree->[$i]->{$k}->{number_r}){
          $writer->dataElement(val => $tree->[$i]->{$k}->{number_r});
        }
        else{
          $writer->dataElement(val => '');
        }
      }
      elsif($data eq 'both'){
        if(exists $tree->[$i]->{$k}->{number} && defined $tree->[$i]->{$k}->{number}){
          $writer->dataElement(val => nearest(0.00001, $tree->[$i]->{$k}->{number}));
        }
        else{
          $writer->dataElement(val => '');
        }
        if(exists $tree->[$i]->{$k}->{number_r} && defined $tree->[$i]->{$k}->{number_r}){
          $writer->dataElement(val => $tree->[$i]->{$k}->{number_r});
        }
        else{
          $writer->dataElement(val => '');
        }
      }
      else{
        if(exists $tree->[$i]->{$k}->{number} && defined $tree->[$i]->{$k}->{number}){
          $writer->dataElement(val => nearest(0.00001, $tree->[$i]->{$k}->{number}));
        }
        else{
          $writer->dataElement(val => '');
        }
      }
    }
    $writer->endTag('magnitude');
    if(defined($self->{_qry_lst})){
      $writer->startTag('blast');
      my $name_s = $k;
      $name_s =~ s/[^A-Za-z0-9]/_/g;
      for(my $i=0;$i<=$#{$index_list};$i++){
        if(defined($self->{_html_files}->[$i]->{$name_s})){
          if($data eq 'contigs'){
            $writer->startTag('val','href' => $self->{_html_files}->[$i]->{$name_s});
            $writer->characters('link');
            $writer->endTag('val');
          }
          elsif($data eq 'reads'){
            $writer->dataElement(val => '');
          }
          elsif($data eq 'both'){
            $writer->startTag('val','href' => $self->{_html_files}->[$i]->{$name_s});
            $writer->characters('link');
            $writer->endTag('val');
            $writer->startTag('val');
            $writer->endTag('val');
          }
          else{

          }
        }
        else{
          if($data eq 'contigs' || $data eq 'reads'){
            $writer->dataElement(val => '');
          }
          elsif($data eq 'both'){
            $writer->dataElement(val => '');
            $writer->dataElement(val => '');
          }
        }
      }
      $writer->endTag('blast');
    }
    if(defined($self->{_extract_sequences})){
      $writer->startTag('sequence');
      #  $keysList{$k} is tax_id
      #  $k is taxonomic name
      for(my $i=0;$i<=$#{$index_list};$i++){
        my $sequences = "sequences.html";
        my $compt = 1;
        if (defined($keysList{$k})) {
          foreach my $inner (keys %{$self->{_extract_sequences}->[$i]->{$keysList{$k}}}){
            if(defined($inner) && !($inner eq "")){
              if($compt == 1){
                $sequences = $sequences . "?name=" . $k . "&?seq" . $compt . "=" . $inner;
              }else{
                $sequences = $sequences . "&?seq" . $compt . "=" . $inner;
              }
              $compt++;
            }
          }
          # $logger->info('inner ' . $sequences);
          foreach my $inner (keys %{$self->{_extract_sequences}->[$i]->{$keysList{$k}}}){
            # $inner is sequence
            if(defined($sequences) && !($sequences eq "sequences.html")){
              if($data eq 'contigs'){
                $writer->startTag('val','href' => $sequences);
                $writer->characters('link');
                $writer->endTag('val');
              }
              elsif($data eq 'reads'){
                $writer->dataElement(val => '');
              }
              elsif($data eq 'both'){
                $writer->startTag('val','href' => $sequences);
                $writer->characters('link');
                $writer->endTag('val');
              }
            }
          }
        }
      }
      $writer->endTag('sequence');
    }
    if ( $color eq 'taxid' ) {
      $writer->startTag('taxid');
      for(my $i=0;$i<=$#{$index_list};$i++){
        $writer->dataElement(val => $keysList{$k});
        if($data eq 'both'){
          $writer->dataElement(val => $keysList{$k});
        }
      }
      $writer->endTag('taxid');
    }
    elsif ( $color eq 'identity' ){
      # $logger->info('color ' . $color);
      $writer->startTag('identity');
      for(my $i=0;$i<=$#{$index_list};$i++){
        if($data eq 'contigs'){
          if(exists $tree->[$i]->{$k} && defined($tree->[$i]->{$k}->{identity})){
            $writer->dataElement(val => $tree->[$i]->{$k}->{identity} / $tree->[$i]->{$k}->{number});
          }
          else{
            $writer->dataElement(val => '');
          }
        }
        elsif($data eq 'reads'){
          if(exists $tree->[$i]->{$k} && defined($tree->[$i]->{$k}->{identity_r})){
            $writer->dataElement(val => $tree->[$i]->{$k}->{identity_r} / $tree->[$i]->{$k}->{number_r});
          }
          else{
            $writer->dataElement(val => '');
          }
        }
        elsif($data eq 'both'){
          if(exists $tree->[$i]->{$k} && defined($tree->[$i]->{$k}->{identity})){
            $writer->dataElement(val => $tree->[$i]->{$k}->{identity} / $tree->[$i]->{$k}->{number});
          }
          else{
            $writer->dataElement(val => '');
          }
          if(exists $tree->[$i]->{$k} && defined($tree->[$i]->{$k}->{identity_r})){
            $writer->dataElement(val => $tree->[$i]->{$k}->{identity_r} / $tree->[$i]->{$k}->{number_r});
          }
          else{
            $writer->dataElement(val => '');
          }
        }
        else{
          if(exists $tree->[$i]->{$k} && defined($tree->[$i]->{$k}->{identity})){
            $writer->dataElement(val => $tree->[$i]->{$k}->{identity} / $tree->[$i]->{$k}->{number});
          }
          else{
            $writer->dataElement(val => '');
          }
        }
      }
      $writer->endTag('identity');

    }
    XMLPrinter($self,$newHash,$writer,$color,$data, $index_list);
    $writer->endTag('node');
  }
}


sub _set_options {
  my ($self)=@_;
  if(scalar(@input_files) > 0){
    if(scalar(@id_samples) == 0){
      foreach my $file (@input_files){
        my($filename, $dirs, $suffix) = fileparse($file);
        my @p = split(/\./,$filename);
        push(@id_samples,$p[0]);
      }
    }
    $self->{id_samples} = \@id_samples;
    $self->{filesList} = \@input_files;
  }
  else{
    $logger->error('You must provide at least one csv file.');
    &help;
  }

  if(scalar(@xml_files) > 0){
    $self->{xml_files} = \@xml_files;
    $self->{_out_dir} = 'krona-data';
  }

  if(scalar(@group_file) > 0){
    $self->{group_file} = \@group_file;
  }

  if (!$color || ($color ne 'identity' && $color ne 'taxid' && $color ne 'none')){
    $color = 'identity';
    $logger->warn("Wrong color selected (must be identity, taxid or none). Set to identity");
  }

  if (!$data || ($data ne 'both' && $data ne 'reads' && $data ne 'contigs' && $data ne 'none')){
    $color = 'none';
    $logger->warn("Wrong data selected. Supported options are reads|contigs|both|none.");
  }
  my %taxonomyParams;
  if(! -e $db){
    $logger->error('Taxonomy sqlite database not found. ' . $db);
  }
  if($dbType eq 'ncbi'){
    $taxonomyParams{'dbh'} = DBI->connect("dbi:SQLite:dbname=$db","","");
  }
  elsif($dbType eq 'gg'){
    $taxonomyParams{'gg'} = 1;
  }
  else{
    $logger->error("Wrong database type. Must be ncbi or gg");
    &help();
  }
  if($reduceTaxo){
    $taxonomyParams{_dbFormat} = 'reduced';
  }
  else{
    $taxonomyParams{_dbFormat} = 'normal';
  }
  $self->{taxoTools} = Tools::Taxonomy->new(%taxonomyParams);
  if($outdir ne ''){
    $self->{_krona_out_dir} = $outdir;
  }
  else{
    $self->{_krona_out_dir} = 'krona_out_dir';
  }
  if($output ne ''){
    $self->{_krona_file_name} = $output;
  }
  else{
    if(defined($merge)){
      $self->{_krona_file_name} = 'merged.html';
    }
    else{
      $self->{_krona_file_name} = $self->{id_samples}->[0] . '.krona.html';
    }
  }

}


sub _extract_blast {
  my ($self,$index)=@_;
  if(! -e $self->{_out_dir} . '/' . $index){
    mkdir $self->{_out_dir} . '/' . $index;
  }
  my $krona_fof = $self->{_out_dir} . '/' . $index . '/' . 'krona_fof.txt';
  open(FOF,">$krona_fof");
  my @clean;
  for my $tax_id (keys(%{$self->{_qry_lst}->[$index]})){
    my $name_s = $self->{taxoTools}->retrieveNameFromTaxId($tax_id);
    $name_s =~ s/[^A-Za-z0-9]/_/g;
    my $html_out = $self->{_out_dir} . '/' . $index . '/' . $name_s . '.html';
    if(! -e $html_out){
      my $csv_out = $self->{_out_dir} . '/' . $index . '/' . $name_s . '.csv';
      if(! -e $csv_out){
        open(CSV,">$csv_out") || $logger->logdie('Cannot create file ' . $csv_out);
        foreach my $id (keys(%{$self->{_qry_lst}->[$index]->{$tax_id}})){
          print CSV $id . "\n";
        }
        close CSV;
        my $xml_out = $self->{_out_dir} . '/' . $index . '/' . $name_s . '.xml';
        print FOF $csv_out . "\t" . $xml_out . "\n";
        push(@clean,$csv_out);
        push(@clean,$xml_out);
      }
    }
  }
  close FOF;

  my $cmd = 'demultiplex-BLAST-results --fof ' . $krona_fof . ' -i ' . $self->{xml_files}->[$index];
  $logger->debug($cmd);
  if(! -z $krona_fof){
    `$cmd`;
  }
  my $h={};
  for my $tax_id (keys(%{$self->{_qry_lst}->[$index]})){
    my $name_s = $self->{taxoTools}->retrieveNameFromTaxId($tax_id);
    $name_s =~ s/[^A-Za-z0-9]/_/g;
    my $xml_out = $self->{_out_dir} . '/' . $index . '/' . $name_s . '.xml';
    my $html_out = $self->{_out_dir} . '/' . $index . '/' . $name_s . '.html';
    if(! -e $html_out){
      if( !(-z $xml_out) && !(-e $xml_out)){
        $cmd = 'blast2html.py -i ' . $xml_out . ' -o ' . $html_out;
        $logger->debug($cmd);
        `$cmd`;
      }
      else{
        $logger->debug('XML file does not exist or empty. ' . $xml_out);
      }

    }
    $h->{$name_s} = $html_out;
  }
  push(@clean,$krona_fof);
  foreach my $el (@clean){
    # `rm $el`;
  }
  return $h;
}


sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sandy Contreras & Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    This script is used to parse csv file containing tax_id field and creates Krona charts.
#

USAGE: perl $prog -i blast_csv_extended_1 -i blast_csv_extended_2 ... -i blast_csv_extended_n [OPTIONS]

       ### OPTIONS ###
       -i|input        <BLAST CSV>=<GROUP FILE>  Blast CSV extended file and CSV group file corresponding to blast (optional)
          Exemples of usage : $prog -i blastFile1=groupFile1 -i blastFile2 ... -i blastFilen=groupFilen
          Group file format is : sequence_regexp <TAB> substitution1 <SPACE> optional_weigth1 <TAB> ... substitutionx <SPACE> optional_weigthx
          If no weigth is provided, is set to 1
          If group file is provided but contig don't match with any regexp, the output number will be the read number for ecsv or 1 if no present (csv or xml blast)
          And the output name will be the contig name.
          Exemples : Having contig1_xxx has query name
          -group file : contig1_xxx <TAB> "read_indiv1" <TAB> "read_indiv2"
          will deconvolute each blast hit of contig1_xxx in 2 blast hits : read_indiv1 1 blast_res, read_indiv2 1 blast_res
          -group file : contig1_xxx <TAB> "read_indiv1" <SPACE> 2 <TAB> "read_indiv2"
          will deconvolute each blast hit of contig1_xxx in 2 blast hits : read_indiv1 2 blast_res, read_indiv2 1 blast_res
          -group file : ([^_]+)_.* <TAB> \$1 will deconvolute each blast hit of contig1_xxx in 1 blast hit : contig1 1 blast_res
          -group file : ([^_]+)_.* <TAB> \$1 <SPACE> 2 will deconvolute each blast hit of contig1_xxx in 1 blast hit : contig1 2 blast_res

       -f|format       <ecsv|csv|xml>   Blast format (ecsv, csv or xml)
       -m|merge        <NAME>  Merge the inputs csv files and generate a common Krona file
          If no name is provided, the default name will be merged.html

       -mt|mergeType   <1|2>  Merge type for blast.
          1 : Group blast results in the same Krona allowing user to switch between the different blast files.
          2 : Merge blast files into a single blast file and generate a Krona using the resulting file.

       -c|color        <STRING>  Coloration mode for Krona chart: identity, taxid or none (Default: $color).
       -data           <STRING>  Data to print in the Krona chart: reads, contigs or both.
       -o|output       <DIRECTORY>  Output directory
       -help|h         Print this help and exit
EOF
exit(1) ;
}