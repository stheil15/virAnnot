package Tools::Tools;

use strict;
use Bio::SeqIO;
use Data::Dumper;
use Logger::Logger;
use Cwd 'abs_path';
use base 'Exporter';

our @EXPORT    = qw(retrieveFeaturesFromFile selectFeaturesOnPrimaryTag readDir retrieveSequence);


=head1 BIOPERL RELATED METHODS

=head2 retrieveFeaturesFromFile

=head3 Arguments

=over 4

=item 1.

A file path.

=item 2.

A file format 'EMBL|GENBANK'.

=back

=head3 Returns

=over 4

=item -

A hash in which the key is the sequence index and values are features and id.
ex: $r->{1}->{features} = \@features 
	$r->{1}->{id} = display_id
	$r->{1}->{seq} = Seq

=back

=cut

sub retrieveFeaturesFromFile {
	my ($file,$format) = @_ ;
	if($format !~ /EMBL|GENBANK/){
		$logger->logdie('Unauthorized format ' . $format);
	}
	my $res;
	my $c=0;
	my $seqIO = Bio::SeqIO->new( -format => $format, -file => $file, -verbose => -1) ;
	foreach my $s ($seqIO->next_seq){
		$res->{++$c}->{features} = [$s->all_SeqFeatures];
		$s->flush_SeqFeatures;
		$res->{$c}->{seq} = $s;
		$res->{$c}->{id} = $s->display_id;
	}
	$logger->debug("DEBUG: " . scalar(keys(%{$res})) . " sequences and features have been retrieve from " . $file . ".");
	return $res ;
}

=head2 selectFeaturesOnPrimaryTag

=head3 Arguments

=over 4

=item 1.

An array of Bio::SeqFeatures.

=item 2.

The primary tag to search for. ex: 'CDS' 

=back

=head3 Returns

=over 4

=item -

An array containing the selected features.

=back

=cut

sub selectFeaturesOnPrimaryTag {
	my ($features, $primaryTag) = @_ ;
	my @selectedFeatures ;
	foreach my $f (@{$features}){
		if($f->primary_tag =~ /$primaryTag/){
			push(@selectedFeatures,$f);
		}
	}
	$logger->debug( scalar(@selectedFeatures) . ' features with \'' . $primaryTag . '\' as primary tag have been found over '  . scalar(@{$features}) . '.');
	return \@selectedFeatures ;
}


=head2 checkNoteValue

=head3 Arguments

=over 4

=item 1.

A Bio:SeqFeatures object.

=item 2.

The tag in which to search for. ex: 'note', 'label'.

=item 3.

The value to search for. ex: 'partially sequenced'.

=back

=head3 Returns

=over 4

=item -

0 or 1. FALSE or TRUE

=back

=cut

sub checkNoteValue{
	my ($feature,$tag,$value) = @_;
	if($feature->has_tag($tag)){
		my @values = $feature->get_tag_values($tag);
		foreach my $value (@values){
			if($value =~ /$value/gi){
				$logger->debug($value . ' have been found in ' . $tag . '.');
				return 1;
			}
		}
	}
	return 0;
}


=head2 retrieveSequence

=head3 Arguments

=over 4

=item 1.

A file path.

=item 2.

A file format 'EMBL|GENBANK|FASTA'.

=back

=head3 Returns

=over 4

=item -

A Bio::Seq object.

=back

=cut

sub retrieveSequence {
	my ($file,$format) = @_ ;
	if($format !~ /^EMBL$|^GENBANK$|^FASTA$/){
		$logger->logdie($format . ' is not a supported format. (EMBL|GENBANK|FASTA)');
	}
	my $res;
	my $seqIO = Bio::SeqIO->new( -format => $format, -file => $file, -verbose => -1) ;
	while(my $s = $seqIO->next_seq){
		$res->{$s->display_id} = $s;
	}
	return $res ;
}

=head2 readDir

Reads a the content of a directory and returns an array of files.

=head3 Arguments

=over 4

=item 1.

A path.

=back

=head3 Returns

=over 4

=item -

An array of files path.

=back

=cut

#~ TODO: get the absolute path of $path.
sub readDir {
	my ($path) = @_ ;
	$path = abs_path($path);
	opendir(PATH, $path) || $logger->error("Can't open directory $path") ;
	my @files = readdir(PATH) ;
	my @tab;
	foreach my $file (@files){
			if($file =~ /^\.\.?$/){next;}
			$file = $path . "/" . $file ;
	$file  = abs_path($file);
			if(! -d $file){next;}
			push(@tab,$file);
	}
	close PATH ;
	return \@tab ;
}


1;
