#!/share/apps/bin/perl
package Tools::Setter;

use strict;
use Bio::SeqIO;
use Data::Dumper;
use Logger::Logger;
use Cwd 'abs_path';
use base 'Exporter';

our @EXPORT    = qw(setFile setPositiveInteger setPourcentage extractFilePathInfo);


=head1 BIOPERL RELATED METHODS

=head2 setFileExistence

Check if the file given in argument exists and not empty.

=head3 Arguments

=over 4

=item 1.

A file name.

=back

=head3 Returns

=over 4

=item -

The absolute path to the file or exit.

=back

=cut

sub setFile {
	my ($file) = @_;
	$file = abs_path($file);
	if(-e $file && ! -z $file ){
		return $file ;
	}
	else{
		$logger->logdie('The file ' . $file . ' do not exist or is empty.');
	}
}


=head2 setPositiveInteger

Check if the number given in argument is a positive integer.

=head3 Arguments

=over 4

=item 1.

A number.

=back

=head3 Returns

=over 4

=item -

The number or exit.

=back

=cut

sub setPositiveInteger {
	my ($int) = @_;
	if($int == abs($int) && $int >= 0){
		return $int
	}
	else{
		$logger->logdie('The number ' . $int . ' is not a positive integer.');
	}
}


=head2 setPourcentage

Check if the number given in argument is between 0 and 100.

=head3 Arguments

=over 4

=item 1.

A number.

=back

=head3 Returns

=over 4

=item -

The pourentage or exit.

=back

=cut

sub setPourcentage {
	my ($int) = @_;
	if($int >= 0 && $int <= 100){
		return $int
	}
	else{
		$logger->logdie('The number ' . $int . ' is not a correct pourcentage.');
	}
}

=head2 extractFilePathInfo

Extract informatins from an absolute file path.

=head3 Arguments

=over 4

=item 1.

An absolute path.

=back

=head3 Returns

=over 4

=item -

Returns a 3 key hash :
	file => the full path to the file.
	path => the full directory path.
	prefix => the prefix of the file name.
	suffix => the suffix of the file name.

=back

=cut


sub extractFilePathInfo{
    my ($path) = @_;
    $path = abs_path($path);
    if(! -e $path){
		$logger->logdie('File not found: ' . $path);
	}
	else{
		my @line = split('/',$path);
		my $file = pop @line ;
		$file =~/(.*)\.(.*)$/;
		my $h = { full => $path, file => $file , path => join('/',@line), prefix => $1, suffix => $2 };
		return $h;
	}
}

1;
