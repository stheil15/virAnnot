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
						"i|info=s"				=> \@infoFiles,
						"c|clean"				  => \$clean,
						"o|output=s"			=> \$outputFile,
						"v|verbosity=i" 	=> \$verbosity,
						"h|help"				=> \$help,
) ;

if($verbosity > 1){
    Logger::Logger->changeMode($verbosity);
}

my @colors = (
	['merge', 'taxonomy', 'Viruses;.*', 'red'],
	['nr', 'taxonomy', 'Viruses;.*', 'red'],
	['merge', 'taxonomy', 'Viruses;dsRNA viruses.*', 'blue'],
	['nr', 'taxonomy', 'Viruses;dsRNA viruses.*', 'blue'],
	['merge', 'taxonomy', 'Viruses;ssRNA viruses.*', 'cyan'],
	['nr', 'taxonomy', 'Viruses;ssRNA viruses.*', 'cyan'],
	['merge', 'taxonomy', 'Viruses;Retro-transcribing viruses.*', 'green'],
	['nr', 'taxonomy', 'Viruses;Retro-transcribing viruses.*', 'green'],
	['merge', 'taxonomy', 'Viruses;dsDNA viruses.*', 'green'],
	['nr', 'taxonomy', 'Viruses;dsDNA viruses.*', 'green'],
	['merge', 'taxonomy', 'Viruses;ssDNA viruses.*', 'lime'],
	['nr', 'taxonomy', 'Viruses;ssDNA viruses.*', 'lime'],
	['merge', 'taxonomy', 'Viroid', 'yellow'],
	['nr', 'taxonomy', 'Viroid', 'yellow'],
	['pfam', 'superkingdom', 'Viruses', 'red'],
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

	# $fileInfos = $self->mergeBlastCSV($fileInfos, $clean);
	# my $workbook = Excel::Writer::XLSX->new( $outputFile );
	# my $worksheet = $workbook->add_worksheet('blast');


	# $self->printExcelTaxoStats($fileInfos, $workbook, 'blast_stats');

	foreach my $rpsFile (@{$self->{rpsFiles}}){
		my ($rpsHits, $headers) = Tools::Taxonomy->readCSVextended($rpsFile, "\t");
		my($filename, $dirs, $suffix) = fileparse($rpsFile);
		my @p = split(/_/,$filename);
		$p[1] =~ s/\.csv//;
		my $worksheet = $workbook->add_worksheet($p[1]);

		$self->printCSVExcel($rpsHits, $worksheet, $headers, \@excelColors);
	}
}


sub mergeBlastCSV {
	my ($self, $files, $removeDuplicates) = @_;
	my $mergedFile= [];
	my %seenSequence;
	my $nbFile = 1;
	foreach my $file ( @$files ){
		foreach my $line ( @$file ){
			if(! $removeDuplicates || (! exists $seenSequence{$line->{'query_id'}} && ! $line->{no_hit}) || $seenSequence{$line->{'query_id'}} == $nbFile){
				push ( @$mergedFile, $line );
				$seenSequence{$line->{'query_id'}} = 1;
				if($line->{taxonomy} =~ /Viruses/){
					$self->{_virus_seq_id}->{$line->{'query_id'}}++;
				}
			}
		}
		$nbFile++;
	}
	return $mergedFile;
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
		$worksheet->write($nbLines,0,\@fields);

		if($color || $level){
			$worksheet->set_row($nbLines, undef, $color, $level, $level);
		}
		$nbLines++;
	}
	# foreach my $file ( @$files ){

		# foreach my $line ( @$files ){
		#

		# }
	# }
	my $maxCell = xl_rowcol_to_cell($nbLines-1, scalar @$headers - 1);
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
		my ($sheet, $field, $regexp, $color) = @$rule;

		if(! exists $formatList{$color}){
		  $formatList{$colors} = $workbook->add_format(bg_color => $color);
		}
		push(@excelColors, [$sheet, $field, $regexp, $formatList{$colors}]);
	}
	return @excelColors;
}


sub getRowColor {
	my ($self, $excelColors, $sheet, $field, $value, $query_id) = @_;
	my $color;
  # print $sheet->get_name() . "\t" . $value . "\n";
	foreach my $rule (@$excelColors){
		my ($sheetRule, $fieldRule, $valueRule, $colorRule) = @$rule;
		if(defined $value && $sheet->get_name() =~ /$sheetRule/ && $field =~ /^$fieldRule$/ && $value =~ /$valueRule/){
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


sub printExcelTaxoStats {
	my ($self, $matches, $workbook, $filename) = @_;

	if(!$filename){
		$filename = '';
	}

	my $worksheet = $workbook->add_worksheet('stats_blast_taxo');
	my $format = $workbook->add_format(valign => 'vcenter', align  => 'center');

	my $lastLine = 3;

	$worksheet->merge_range( 'A1:A3', 'File', $format );
	$worksheet->merge_range( 'B1:B3', 'No_hit', $format );
	$worksheet->merge_range( 'C1:O1', 'Annotated', $format );
	$worksheet->merge_range( 'C2:C3', 'No taxonomy', $format );
	$worksheet->merge_range( 'D2:D3', 'Viroids', $format );
	$worksheet->merge_range( 'E2:G2', 'Cellular organisms', $format );
	$worksheet->merge_range( 'H2:L2', 'Viral taxonomy', $format );
	$worksheet->merge_range( 'M2:M3', 'Other', $format );
	$worksheet->merge_range( 'N2:N3', 'Total viruses', $format );
	$worksheet->merge_range( 'O2:O3', 'Total no viruses', $format );
	$worksheet->write(2,4,['Eukaryota', 'Archaea', 'Bacteria']);
	$worksheet->write(2,7,['dsRNA', 'ssRNA', 'dsDNA', 'ssDNA', 'Undeterminated']);

	my @keys = ('no_hit', 'no_taxonomy', 'viroids', 'eukaryota', 'archaea', 'bacteria', 'dsRNA_viruses', 'ssRNA_viruses', 'dsDNA_viruses', 'ssDNA_viruses', 'undeterminated_viruses', 'other', 'viruses', 'no_viruses');

	my $stats = Tools::Blast->getBlastStats($matches);
	my @line;
	push(@line,$filename);
	foreach my $field (@keys){
		push(@line,$stats->{$field});
	}
	$worksheet->write($lastLine,0,\@line);
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

	@{$self->{infoFiles}} = map{glob $_} @infoFiles;
	if(scalar(@{$self->{infoFiles}}) > 0){
		foreach my $file (@{$self->{infoFiles}}){
			if(! -e $file){
				$logger->error('Info file not found. ' . $file);
        &help;
			}
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
# PURPOSE:    This script is used to parse blast+ output and produce a m8 output
#             extended with taxonomy information.
#

USAGE: perl $prog -b blast_csv_extended_1 -b blast_csv_extended_2 ... -b blast_csv_extended_n [OPTIONS]

          ### OPTIONS ###
          -b|blast  <BLAST ECSV>  Blast ECSV file.
					-o|output
          -i|info   <CSV FILE>    Supplementary csv info file (rsp blast ...).
          Info files will be added in all the excel files and the sheet name will be equal to the filename
          -c|clean	      If multiple blast are provided and a query match, keep only one hit.
          -v|verbosity    Verbosity level.
          -help|h         Print this help and exit.


EOF
exit(1) ;
}
