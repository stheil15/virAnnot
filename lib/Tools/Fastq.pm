package Tools::Fastq;

use strict;
use warnings;
use Logger::Logger;
use Storable;


=head1 INDEXED FASTQ RELATED METHODS

=head2

=head2 new

=head2

=head3 Description

Create a new Tools::Fastq object and index the FASTQ file

=head3 Arguments

=over 4

=item

A hash of parameters.

Currently accepted keys are :

'file' => FASTQ file path

=back

=head3 Returns

=over 4

=item

A Tools::Fastq object

=back

=cut

sub new {
	my ($class, %attrs) = @_;
	my $self = {};
	bless $self;
	if(defined($attrs{file})){
		$self->{file} = $attrs{file};
		open($self->{file_handle},$self->{file}) || $logger->logdie('Error opening file : '. $self->{file}.' : '.$!."\n");
		$self->indexFastqFile;
	}
	return $self;
}

=head2 indexFastqFile

=head2

=head3 Description

Index a FASTQ file creating a hash reference with the following structure :

$index -> {seq_id} = {'id_begin_position' => integer, 'id_length' => integer}

For each sequence id, the "@" symbol and all the text after space will be removed.

This cleaned id will be used as key for the index.

=head3 Arguments

=over 4

=item

None

=back

=head3 Returns

=over 4

=item

None

=back

=cut

sub indexFastqFile{

	  my ($self) = @_;
    $logger->info('Indexing file : '.$self->{file}."\n");
    my $index;
    my $id;
    my $id_begin_position = 0;
		my $fh = $self->{file_handle};
    while(my $line = <$fh>){

        if($line =~ /^@(\S+)/){

            $id = $1;
            chomp $id;
            $index -> {$id} = {'id_begin_position' => $id_begin_position, 'id_length' => length $line};
            $logger->trace('Indexing sequence' . $id . ' (position_begin_id : '. $index -> {$id}{'id_begin_position'} . ', id_length : '. $index -> {$id}{'id_length'} .') from ' . $self->{file} . "\n");
            <$fh>; <$fh>; <$fh>;
        }

        $id_begin_position = tell($fh);
    }

    $logger->info('File '.$self->{file}.' is now indexed (index contains '.(scalar keys %$index)." sequences)\n");
    $self->{index} = $index;
}

=head2 loadFastqIndexFile

=head2

=head3 Description

Retrieve index from file using Storable module

=head3 Arguments

=over 4

An index file

=back

=head3 Returns

=over 4

=item

A hash reference corresponding to the index of the input FASTQ file :

$index -> {seq_id} = {'id_begin_position' => integer, 'id_length' => integer}

=back

=cut

sub loadFastqIndexFile{

	my ($self, $file) = @_;
  $self->{index} = retrieve($file);
  $logger->info('File '.$file." is now loaded\n");
}

=he=head2 writeFastaIndexFile

=head2

=head3 Description

Write index to file using Storable module

=head3 Arguments

=over 4

=item

A hash reference corresponding to FASTQ index.

=item

An output file path where to store the index.

=back

=head3 Returns

=over 4

=item

The output file path containing index

=back

=cut

sub writeFastqIndexFile{
	my ($self, $file) = @_;
  $logger->info('Writing index ('.(scalar keys %{$self->{index}}).' sequences) in file : '.$file."\n");
  store $self->{index}, $file;
  $logger->info('File '.$file." is now created\n");
}

=head2 retrieveFastqSequence

=head2

=head3 Description

Retrieve FASTQ sequences using a list of ids

=head3 Arguments

=over 4

=item

A sequence id OR an array reference containing the list of sequences id to retrieve.

=back

=head3 Returns

=over 4

=item

A hash reference containing sequences id as keys and sequences as values

$data -> {seq_id} = sequence_corresponding_to_seq_id

=back

=cut

sub retrieveFastqSequence {
	my ($self, $ids) = @_;
	my $data={};
  my $nbSequences = 0;
  if(! ref $ids){$ids = [$ids]}
  $logger->debug('Retrieving sequences of '.scalar(@$ids).' ids from indexed file : '.$self->{file}."\n");
	my $fh = $self->{file_handle};
	foreach my $id (@$ids){
		my $cleanedId = $id;
      if($id =~ /@(\S+)/){$cleanedId = $1}
      $logger->trace('Retrieving informations of id ' . $cleanedId. " from index\n");
      if(exists $self->{index} -> {$cleanedId}){
        $logger->trace('id ' . $cleanedId . ' is present in index (id_begin_position : '. $self->{index} -> {$cleanedId}{'id_begin_position'}. ', id_length : '. $self->{index} -> {$cleanedId}{'id_length'}.")\n");
        seek($fh, $self->{index} -> {$cleanedId}{'id_begin_position'}, 0);
        <$fh>;
        my $sequence = <$fh>;
        $data->{$id} = $sequence;
        $nbSequences ++;
      	$logger->trace('Sequence of id '.$cleanedId.' is : ' . $sequence . "\n")
    	}
    	else{
				$logger->trace('id ' . $cleanedId. " not found in index\n")
			}
	}
  $logger->debug($nbSequences.'/'.scalar(@$ids).' sequences has been retrieved from indexed file ' . $self->{file} . "\n");
	return $data;
}

=head2 retrieveFastqQuality

=head2

=head3 Description

Retrieve FASTQ sequences quality using a list of ids

=head3 Arguments

=over 4

=item

A sequence id OR an array reference containing the list of sequences id to retrieve quality.

=back

=head3 Returns

=over 4

=item

A hash reference containing sequences id as keys and sequences quality as values

$data -> {seq_id} = sequence_quality_corresponding_to_seq_id

=back

=cut

sub retrieveFastqQuality {
	my ($self, $ids) = @_;
	my $data;
  my $nbSequences = 0;
  if(! ref $ids){$ids = [$ids]}
  $logger->debug('Retrieving sequence quality of '.scalar(@$ids).' ids from indexed file : '.$self->{file}."\n");
	my $fh = $self->{file_handle};
	foreach my $id (@$ids){
		my $cleanedId = $id;
		if($id =~ /@(\S+)/){
			$cleanedId = $1;
		}
    $logger->trace('retrieving informations of id ' . $cleanedId. " from index\n");
    if(exists $self->{index} -> {$cleanedId}){
      $logger->trace('id ' . $cleanedId . ' is present in index (id_begin_position : '. $self->{index} -> {$cleanedId}{'id_begin_position'}. ', id_length : '. $self->{index} -> {$cleanedId}{'id_length'}.")\n");
      seek($fh, $self->{index} -> {$cleanedId}{'id_begin_position'}, 0);
      my $quality .= <$fh>.<$fh>.<$fh>;
      $quality = <$fh>;
      $data .= $quality;
      $nbSequences ++;
      $logger->trace('Sequence quality of id '.$cleanedId.' is : ' . $quality. "\n")
    }
    else{
			$logger->trace('id ' . $cleanedId. " not found in index\n");
		}
	}
  $logger->debug($nbSequences.'/'.scalar(@$ids).' sequences qualities has been retrieved from indexed file ' . $self->{file} . "\n");
	return $data;
}

=head2 retrieveFastqBlock

=head2

=head3 Description

Retrieve FASTQ formatted sequences using a list of ids

=head3 Arguments

=over 4

=item

A sequence id OR an array reference containing the list of sequences id to retrieve.

=back

=head3 Returns

=over 4

=item

A scalar containing the sequences corresponding to ids in FASTQ format

=back

=cut

sub retrieveFastqBlock {
	my ($self, $ids) = @_;
	my $data;
  my $nbSequences = 0;

  if(! ref $ids){$ids = [$ids]}

  $logger->debug('Retrieving fastq block of '.scalar(@$ids).' ids from indexed file : '.$self->{file}."\n");
	my $fh = $self->{file_handle};
	foreach my $id (@$ids){
		my $cleanedId = $id;
		if($id =~ /@(\S+)/){
			$cleanedId = $1;
		}
    $logger->trace('Retrieving informations of id ' . $cleanedId. " from index\n");
    if(exists $self->{index} -> {$cleanedId}){
      $logger->trace('id ' . $cleanedId . ' is present in index (id_begin_position : '. $self->{index} -> {$cleanedId}{'id_begin_position'}. ', id_length : '. $self->{index} -> {$cleanedId}{'id_length'}.")\n");
      seek($fh, $self->{index} -> {$cleanedId}{'id_begin_position'}, 0);
    	read($fh, my $block, $self->{index} -> {$cleanedId}{'id_length'});
      $block .= <$fh>.<$fh>.<$fh>;
      $data .= $block;
      $nbSequences++;
      $logger->trace('fastq block of id '.$cleanedId.' is : ' ."\n". $block. "\n");
    }
    else{
			$logger->trace('id ' . $cleanedId. " not found in index\n");
		}
	}
  $logger->debug($nbSequences.'/'.scalar(@$ids).' fastq block has been retrieved from indexed file ' . $self->{file} . "\n");
	return $data;
}

sub printFastqBlock {
  my ($self, $ids, $fh_out) = @_;
  my $nbSequences = 0;
  if(! ref $ids){$ids = [$ids]}
  $logger->debug('Printing fastq block of '.scalar(@$ids).' ids from indexed file : '.$self->{file}."\n");
  my $fh = $self->{file_handle};
  foreach my $id (@$ids){
		my $cleanedId = $id;
		if($id =~ /@(\S+)/){
			$cleanedId = $1;
		}
    $logger->trace('Retrieving informations of id ' . $cleanedId. " from index\n");
    if(exists $self->{index} -> {$cleanedId}){
      $logger->trace('id ' . $cleanedId . ' is present in index (id_begin_position : '. $self->{index} -> {$cleanedId}{'id_begin_position'}. ', id_length : '. $self->{index} -> {$cleanedId}{'id_length'}.")\n");
      seek($fh, $self->{index} -> {$cleanedId}{'id_begin_position'}, 0);
    	read($fh, my $block, $self->{index} -> {$cleanedId}{'id_length'});
      $block .= <$fh>.<$fh>.<$fh>;
      print $fh_out $block;
      $nbSequences++;
      $logger->trace('fastq block of id '.$cleanedId.' is : ' ."\n". $block. "\n");
    }
    else{
			$logger->trace('id ' . $cleanedId. " not found in index\n");
		}
  }
  $logger->debug($nbSequences.'/'.scalar(@$ids).' fastq block has been retrieved from indexed file ' . $self->{file} . "\n");
}

1;
