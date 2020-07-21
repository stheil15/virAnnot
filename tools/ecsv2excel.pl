#!/usr/bin/perl

use strict ;
use Getopt::Long;
use Tools::Taxonomy;
use Tools::Blast;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use File::Basename;
use File::Path qw(make_path);
use Pod::Usage;
use Data::Dumper;
use Logger::Logger;

my @blastFiles;
my @rpsFiles;
my @infoFiles;
my $verbosity=1;
my $clean;
my $outputFile="ecsv_excel.xlsx";
my $help,

my $VERSION = '1.2' ;
my $lastmodif = '2015-10-27' ;

GetOptions( "b|blast=s"				=> \@blastFiles,
						"r|rps=s"				  => \@rpsFiles,
						"c|clean"				  => \$clean,
						"o|output=s"			=> \$outputFile,
						"v|verbosity=i" 	=> \$verbosity,
						"h|help"				=> \$help,
) ;

if($verbosity > 1){
    Logger::Logger->changeMode($verbosity);
}

my @colors = (
	['taxonomy', 'Viruses;.*', '#FF0000'],
	['taxonomy', 'Viruses;dsRNA viruses.*', '#3366FF'], # blue
	['taxonomy', 'Viruses;ssRNA viruses.*', '#00CCFF'], # cyan
	['taxonomy', 'Viruses;Retro-transcribing viruses.*', 'green'],
	['taxonomy', 'Viruses;dsDNA viruses.*', 'orange'],
	['taxonomy', 'Viruses;ssDNA viruses.*', 'lime'],
	['taxonomy', 'Viroid', 'yellow'],
	['superkingdom', 'Viruses', '#FF0000'],
);

&main;

sub main {
	my $self={};
	bless $self;

	_set_options($self);
	my @blastExtentedHeaders;
	my $merge = [];
	my $workbook = Excel::Writer::XLSX->new( $outputFile );
	my @excelColors = $self->setExcelColor($workbook, \@colors);
	foreach my $file (@{$self->{blastFiles}}){
		my ($hits, $header) = Tools::Taxonomy->readCSVextended($file, "\t");

		if(scalar @$header > scalar @blastExtentedHeaders){
			@blastExtentedHeaders = @$header;
		}
		my($filename, $dirs, $suffix) = fileparse($file);
		my @p = split(/_/,$filename);
		$p[1] =~ s/\.csv//;
		my $worksheet = $workbook->add_worksheet($p[1]);
		$self->printCSVExcel($hits, $worksheet, \@blastExtentedHeaders, \@excelColors, ['algo', 'query_id']);

	}

	foreach my $rpsFile (@{$self->{rpsFiles}}){
		my ($rpsHits, $headers) = Tools::Taxonomy->readCSVextended($rpsFile, "\t");
		my($filename, $dirs, $suffix) = fileparse($rpsFile);
		my @p = split(/_/,$filename);
		$p[1] =~ s/\.csv//;
		my $worksheet = $workbook->add_worksheet($p[1]);
		$self->printCSVExcel($rpsHits, $worksheet, $headers, \@excelColors);
	}
}


sub printCSVExcel {
	my ($self, $hits, $worksheet, $headers, $colors, $groupByFields) = @_;

	my $nbLines = 1;
	my %lastFields;
	if( $headers && ! ref $headers ){
		$headers = [$headers];
	}
	if($groupByFields){
		@lastFields{@$groupByFields} = ();
	}
	foreach my $hit (@{$hits}){
	    # if(defined($hit->{'taxonomy'})){
	    #   if($hit->{'taxonomy'} eq 'unknown' || $hit->{'taxonomy'} eq ''){
	    #     next;
	    #   }
	    # }
		my ($color, $level) = (undef, defined $groupByFields);
		if(! $headers || ! @$headers){
			$headers = [keys %$hit];
		}
		if($nbLines == 1){
			$worksheet->write(0,0,$headers);
		}

		if(defined($hit->{taxonomy})){
			if($hit->{taxonomy} =~ /Viruses/){
				$self->{_virus_seq_id}->{$hit->{'query_id'}}++;
			}
		}

    my @fields;
		foreach my $field (@$headers){
			push(@fields,$hit->{$field});
			if(! defined $color){
        if($field eq 'taxonomy'){
          $color = $self->getRowColor($colors, $worksheet, $field, $hit->{$field},$hit->{'query_id'});
        }
			}
			if(defined $groupByFields && exists $lastFields{$field} && (defined $hit->{$field} || defined $lastFields{$field}) && (! defined $lastFields{$field} || $lastFields{$field} ne $hit->{$field})){
				$level = 0;
				$lastFields{$field} = $hit->{$field};
			}
		}
    if($color){
      $worksheet->write($nbLines,0,\@fields,$color);
    }
    else{
      $worksheet->write($nbLines,0,\@fields);
    }

		if($level){
			$worksheet->set_row($nbLines, undef, undef, $level, $level);
		}
		$nbLines++;
	}

	my $maxCell = xl_rowcol_to_cell($nbLines-1, scalar(@$headers) - 1);
	my $table = 'A1:' . $maxCell ;
	$worksheet->autofilter( $table );
	$worksheet->freeze_panes( 1, 1 );
}


sub setExcelColor {
	my ($self, $workbook, $colors) = @_;

	my %formatList;
	my @excelColors;

	foreach my $rule (@{$colors}){
    # print Dumper @$rule . "\n";
		my ($field, $regexp, $color) = @$rule;

		if(! exists $formatList{$color}){
		  $formatList{$colors} = $workbook->add_format(bg_color => $color);
		}
		push(@excelColors, [$field, $regexp, $formatList{$colors}]);
	}
	return @excelColors;
}


sub getRowColor {
	my ($self, $excelColors, $sheet, $field, $value, $query_id) = @_;
	my $color;
  # print $sheet->get_name() . "\t" . $value . "\n";
	foreach my $rule (@$excelColors){
		my ($fieldRule, $valueRule, $colorRule) = @$rule;
		if(defined $value && $field =~ /^$fieldRule$/ && $value =~ /$valueRule/){
			if($sheet->get_name() =~ /rps/){
				if(! defined($self->{_virus_seq_id}->{$query_id})){
					$color = $colorRule;
				}
			}
			else{
				$color = $colorRule;
			}
		}
	}
  # print $color . "\n";
	return $color;
}


sub _set_options {
	my ($self)=@_;
	@{$self->{blastFiles}} = map{glob $_} @blastFiles;
  if(scalar(@{$self->{blastFiles}}) > 0){
    foreach my $file (@{$self->{blastFiles}}){
      if(! -e $file){
        $logger->error('Blast file not found. ' . $file);
        &help;
      }
    }
  }
  else{
    $logger->error('You must provide at least one Blast file.');
    &help;
  }

	@{$self->{rpsFiles}} = map{glob $_} @rpsFiles;
  if(scalar(@{$self->{rpsFiles}}) > 0){
    foreach my $file (@{$self->{rpsFiles}}){
      if(! -e $file){
        $logger->error('Blast file not found. ' . $file);
        &help;
      }
    }
  }

}


# sub mergeBlastCSV {
# 	my ($self, $files, $removeDuplicates) = @_;
# 	my $mergedFile= [];
# 	my %seenSequence;
# 	my $nbFile = 1;
# 	foreach my $file ( @$files ){
# 		foreach my $line ( @$file ){
# 			if(! $removeDuplicates || (! exists $seenSequence{$line->{'query_id'}} && ! $line->{no_hit}) || $seenSequence{$line->{'query_id'}} == $nbFile){
# 				push ( @$mergedFile, $line );
# 				$seenSequence{$line->{'query_id'}} = 1;
# 				if($line->{taxonomy} =~ /Viruses/){
# 					$self->{_virus_seq_id}->{$line->{'query_id'}}++;
# 				}
# 			}
# 		}
# 		$nbFile++;
# 	}
# 	return $mergedFile;
# }


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

USAGE: perl $prog -b blast_csv_extended_1 -b blast_csv_extended_2 ... -b blast_csv_extended_n [OPTIONS]

          ### OPTIONS ###
          -b|blast  <BLAST ECSV>  Blast ECSV file.
					-o|output
          -r|rps    <CSV FILE>    Supplementary csv info file (rsp blast ...).
          -c|clean	      If multiple blast are provided and a query match, keep only one hit.
          -v|verbosity    Verbosity level.
          -help|h         Print this help and exit.


EOF
exit(1) ;
}
