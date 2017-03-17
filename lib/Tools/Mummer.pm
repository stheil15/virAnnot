#!/usr/bin/perl

package Tools::Mummer;
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename;
use Cwd 'abs_path';
use Tools::Fasta;
use Logger::Logger;

# Si on appelle la classe depuis un terminal

unless (caller) {
  my $referenceFile = "test/reference.fasta";
  my $queryFile = "test/query.fasta";
  my $deltaFile = "test/test.delta";
  my $mummer = Tools::Mummer -> new ({query => $queryFile, reference => $referenceFile, alignments =>  $deltaFile});
  #my %parameters = $mummer->parameters;
  my @alignments = $mummer->getAlignmentsFromDeltaFile($deltaFile);
  $mummer->generateMultiFastaAlignment;
}

=head1 MUMMER RELATED METHODS

=head2

=head2 new

=head2

=head3 Description

Create a new Tools::Mummer

=head3 Arguments

=over 4

=item

A Mummer delta file or a hash of parameters.

Currently accepted keys for the hash of parameters are :

'alignments' => Mummer delta file

'query' => Query FASTA file

'reference' => Reference FASTA file

=back

=head3 Returns

=over 4

=item

A Tools::Mummer object

=back

=cut

sub new {
  my ($class, $parameters) = @_;
  $class = ref($class) || $class;
  my $self = {};
  bless( $self, $class );
  my %defaultParameters = $class->getDefaultParameters;
  if($parameters){
    if(!ref $parameters){
      $parameters = {'alignments' => $parameters}
    }
    %defaultParameters = (%defaultParameters, %$parameters);
    $self->reference($defaultParameters{'reference'});
    $self->query($defaultParameters{'query'});
    $self->delta($defaultParameters{'alignments'});
    #$self->parameters($defaultParameters{'mummer_parameters'});
  }
  return $self;
}

=head2 getDefaultParameters

=head2

=head3 Description

Get default parameters used to create Mummer object

=head3 Arguments

=over 4

=item

None

=back

=head3 Returns

=over 4

=item

A hash of parameters.

Current keys are :

'alignments' => Mummer delta file

'query' => Query FASTA file

'reference' => Reference FASTA file

=back

=cut

sub getDefaultParameters {
  my ($class) = @_;
  my %defaultParameters = (
                  'mummer_binary' 	  => 'nucmer',
                  'reference' 	      => '',
                  'query' 		        => '',
                  'alignments'	 	    => '',
  );
  return %defaultParameters;
}

=head2 query

=head2

=head3 Description

Get/set query

=head3 Arguments

=over 4

=item

Path to FASTA query file (optional)

Will override query path retrieved from delta file if defined and not empty

=back

=head3 Returns

=over 4

=item

Path to FASTA query file

=back

=cut

sub query {
  my ($self, $value) = @_;
  if(defined $value){
    $self->{'mummer_query'} = $value;
  }
  return $self->{'mummer_query'};
}

=head2 reference

=head2

=head3 Description

Get/set reference

=head3 Arguments

=over 4

=item

Path to FASTA reference file (optional)

Will override reference path retrieved from delta file if defined and not empty

=back

=head3 Returns

=over 4

=item

Path to FASTA reference file

=back

=cut

sub reference {
  my ($self, $value) = @_;
  if(defined $value){
    if(-e $value){
      $self->{'mummer_reference'} = $value;
    }
    else{
      $logger->error('Reference sequence file do not exist. ' . $value);
    }
  }
  return $self->{'mummer_reference'};
}

=head2 delta

=head2

=head3 Description

Get/set delta

=head3 Arguments

=over 4

=item

Path to Mummer delta file (optional)

=back

=head3 Returns

=over 4

=item

Path to Mummer delta file

=back

=cut

sub delta {
  my ($self, $value) = @_;
  if(defined $value){
    $self->{'mummer_delta'} = $value;
  }
  return $self->{'mummer_delta'};
}

=head2 binary

=head2

=head3 Description

Get/set binary

=head3 Arguments

=over 4

=item

Valid Mummer binary : nucmer or promer (optional)

=back

=head3 Returns

=over 4

=item

Mummer binary

=back

=cut

sub binary {
  my ($self, $value) = @_;
  if(defined $value){
    if($value eq 'nucmer' || $value eq 'promer'){
      $self->{'mummer_binary'} = $value;
    }
    else{
      print "Error : $value mummer parameters must be nucmer or promer\n";
    }
  }
  return $self->{'mummer_binary'};
}

#sub parameters{
#  my ($self, $value) = @_;
#  if(defined $value){
#    if(ref $value eq 'HASH') {
#      $self->{'mummer_parameters'} = $value;
#    }
#    else{
#      print "Error : mummer parameters must be an hash ref\n";
#    }
#  }
#  return %{$self->{'mummer_parameters'}};
#}

=head2 getSequences

=head2

=head3 Description

Retrieve sequences from reference and query FASTA files

=head3 Arguments

=over 4

=item

None

=back

=head3 Returns

=over 4

=item

A hash containing sequence ids as keys and complete sequences as values

=back

=cut

sub getSequences {
  my ($self) = @_;
  my $idSeq;
  my %sequences;
  foreach my $fastaFile ($self->reference, $self->query){
    open (FASTA, $fastaFile) or die "Error opening FASTA file $fastaFile : $!\n";
    while (my $line = <FASTA>){
      chomp $line;
      if ($line =~ />([^\s]*)/){
        $idSeq = $1;
      }
      else{
        $sequences{$idSeq} .= $line;
      }
    }
    close FASTA;
  }
  return %sequences;
}

=head2 getAlignmentsFromDeltaFile

=head2

=head3 Description

Retrieve sequences alignments from a delta file or a list of delta file

=head3 Arguments

=over 4

=item

Get alignments from delta files

=back

=head3 Returns

=over 4

=item

An array containing a list of alignments

Each alignment is stored as a hashref and contains the following keys :

REFERENCE => reference id

QUERY => query id

ALIGNMENT_LENGTH_REFERENCE => alignment length on reference sequence

ALIGNMENT_LENGTH_QUERY => alignment length on query sequence

BEGIN_ALIGNMENT_REFERENCE => beginning of alignment on reference sequence

END_ALIGNMENT_REFERENCE => end of alignment on reference sequence

BEGIN_ALIGNMENT_QUERY => beginning of alignment on query sequence

END_ALIGNMENT_QUERY => end of alignment on query sequence

MISMATCHES => number of mitmaches in alignment between reference and query

SIMILARY_SCORE_BELOW_1 => number of similarity score below 1 in alignment

STOP_CODONS => number of stop codons in alignment

QUERY_FRAME => Query frame (1 = sense, -1 = anti-sense)

REFERENCE_FRAME => Reference frame (1 = sense, -1 = anti-sense)

=back

=cut

sub getAlignmentsFromDeltaFile {
  my ($self, $deltaFile) = @_;
  my @alignments;
  if(!ref $deltaFile){
    $deltaFile = [$deltaFile];
  }
  foreach my $deltaFile (@$deltaFile){
    my %minPosReferences;
    open (DELTA, $deltaFile) or die "Error opening delta file $deltaFile : $!\n";
    my $alignment;
    my ($reference, $query, $align_length_reference, $align_length_query);
    my ($referenceFile, $queryFile) = split(" ", <DELTA>);
    my ($binary) = split(" ", <DELTA>);
    if(! $self->reference){
      $self->reference($referenceFile);
    }
    if(! $self->query){
      $self->query($queryFile);
    }
    #if(! $self->binary){
      $self->binary(lc($binary));
    #}
    #my %sequences = $self->getSequences;
    my $ref_seq = Tools::Fasta->new(file => $self->reference);
    my $qry_seq = Tools::Fasta->new(file => $self->query);
    while (my $line = <DELTA>) {
      chomp $line;
      if ($line =~ />([^\s]*)\s([^\s]*)\s([^\s]*)\s([^\s]*)/){
        ($reference, $query, $align_length_reference, $align_length_query) = ($1, $2, $3, $4);
        if (! exists $minPosReferences{$reference} ){
          $minPosReferences{$reference} = 0;
        }
      }
      elsif($line =~ /([^\s]*)\s([^\s]*)\s([^\s]*)\s([^\s]*)\s([^\s]*)\s([^\s]*)\s([^\s]*)/){
        my ($beginAlignRef, $endAlignRef, $beginAlignQuery, $endAlignQuery, $numberOfMismatchesInAlign, $numberOfSimilarityScoresBelow1InAlign, $numberOfStopCodonsInAlign) = ($1, $2, $3, $4, $5, $6, $7);
        $alignment = {};
        push (@alignments, $alignment);
        $alignment->{'REFERENCE'} = $reference;
        $alignment->{'QUERY'} = $query;
        $alignment->{'ALIGNMENT_LENGTH_REFERENCE'} = $align_length_reference;
        $alignment->{'ALIGNMENT_LENGTH_QUERY'} = $align_length_query;
        $alignment->{'BEGIN_ALIGNMENT_REFERENCE'} = $1;
        $alignment->{'END_ALIGNMENT_REFERENCE'} = $2;
        $alignment->{'BEGIN_ALIGNMENT_QUERY'} = $3;
        $alignment->{'END_ALIGNMENT_QUERY'} = $4;
        $alignment->{'MISMATCHES'} = $5;
        $alignment->{'SIMILARITY_SCORE_BELOW_1'} = $6;
        $alignment->{'STOP_CODONS'} = $7;
        $alignment->{'QUERY_FRAME'} = 1;
        $alignment->{'REFERENCE_FRAME'} = 1;
        my $ref_seq_hash = $ref_seq->retrieveFastaSequence($alignment->{'REFERENCE'});
        my $qry_seq_hash = $qry_seq->retrieveFastaSequence($alignment->{'QUERY'});
        if ( $alignment->{'END_ALIGNMENT_REFERENCE'} < $alignment->{'BEGIN_ALIGNMENT_REFERENCE'}) {
          ($alignment->{'END_ALIGNMENT_REFERENCE'}, $alignment->{'BEGIN_ALIGNMENT_REFERENCE'}) = (length($ref_seq_hash->{$alignment->{'REFERENCE'}}) - $alignment->{'END_ALIGNMENT_REFERENCE'} + 1, length($ref_seq_hash->{$alignment->{'REFERENCE'}}) - $alignment->{'BEGIN_ALIGNMENT_REFERENCE'} + 1);
          $alignment->{'REFERENCE_FRAME'} = -1;
        }
        if ( $alignment->{'END_ALIGNMENT_QUERY'} < $alignment->{'BEGIN_ALIGNMENT_QUERY'}) {
          ($alignment->{'END_ALIGNMENT_QUERY'}, $alignment->{'BEGIN_ALIGNMENT_QUERY'}) = (length($qry_seq_hash->{$alignment->{'QUERY'}}) - $alignment->{'END_ALIGNMENT_QUERY'} + 1, length($qry_seq_hash->{$alignment->{'QUERY'}}) - $alignment->{'BEGIN_ALIGNMENT_QUERY'} + 1);
          $alignment->{'QUERY_FRAME'} = -1;
        }
        $alignment->{'DELTA_ALIGNMENT'} = [];
      }
      elsif (defined $line && $line ne '0'){
        push (@{$alignment->{'DELTA_ALIGNMENT'}}, $line);
      }
    }
    close DELTA;
  }
  return @alignments;
}

=head2 cleanFrames

=head2

=head3 Description

Select the best frame for each reference and each query mapped on a given reference and return alignments whose query and refererence frames correspond to best frames.

Note that query best frame is dependent to the alignment reference (can be 1 on ref1 and -1 on ref2).

Best frames are determinated by counting the number of alignments corresponding to each couple of query_frame / reference_frame.

If two best frame have the same number of alignments, the select couple of frame for a given reference and query will be the first in appear in the delta file.

=head3 Arguments

=over 4

=item

A reference array containing Mummer alignements.

=back

=head3 Returns

=over 4

=item

An array containing the list of correctly framed alignments.

=back

=cut

sub cleanFrames {

  my ($self, $alignments) = @_;

  my %frames;
  my $firstSeen = 1;

  foreach my $alignment (@$alignments){
    $frames{$alignment->{'REFERENCE'}}{$alignment->{'REFERENCE_FRAME'}}{'QUERIES'}{$alignment->{'QUERY'}}{$alignment->{'QUERY_FRAME'}}{'TOTAL'} += $alignment->{'END_ALIGNMENT_QUERY'} - $alignment->{'BEGIN_ALIGNMENT_QUERY'} + 1;
    $frames{$alignment->{'REFERENCE'}}{$alignment->{'REFERENCE_FRAME'}}{'TOTAL'} += $alignment->{'END_ALIGNMENT_QUERY'} - $alignment->{'BEGIN_ALIGNMENT_QUERY'} + 1;

    if(! $frames{$alignment->{'REFERENCE'}}{$alignment->{'REFERENCE_FRAME'}}{'FIRST_SEEN'}){
      $frames{$alignment->{'REFERENCE'}}{$alignment->{'REFERENCE_FRAME'}}{'FIRST_SEEN'} = $firstSeen;
    }
    if(! $frames{$alignment->{'REFERENCE'}}{$alignment->{'REFERENCE_FRAME'}}{'QUERIES'}{$alignment->{'QUERY'}}{$alignment->{'QUERY_FRAME'}}{'FIRST_SEEN'}){
      $frames{$alignment->{'REFERENCE'}}{$alignment->{'REFERENCE_FRAME'}}{'QUERIES'}{$alignment->{'QUERY'}}{$alignment->{'QUERY_FRAME'}}{'FIRST_SEEN'} = $firstSeen;
    }
    $firstSeen ++;
  }
  my %bestFrames;

  foreach my $reference (keys %frames){

    foreach my $referenceFrame (sort {$frames{$reference}{$b}{'TOTAL'} <=> $frames{$reference}{$a}{'TOTAL'} || $frames{$reference}{$b}{'FIRST_SEEN'} <=> $frames{$reference}{$a}{'FIRST_SEEN'}} keys %{$frames{$reference}}){

      foreach my $query (keys %{$frames{$reference}{$referenceFrame}{'QUERIES'}}){

        foreach my $queryFrame (sort {$frames{$reference}{$referenceFrame}{'QUERIES'}{$query}{$b}{'TOTAL'} <=> $frames{$reference}{$referenceFrame}{'QUERIES'}{$query}{$a}{'TOTAL'} || $frames{$reference}{$referenceFrame}{'QUERIES'}{$query}{$b}{'FIRST_SEEN'} <=> $frames{$reference}{$referenceFrame}{'QUERIES'}{$query}{$a}{'FIRST_SEEN'}} keys %{$frames{$reference}{$referenceFrame}{'QUERIES'}{$query}}){
          $bestFrames{$reference}{$query}{'reference'} = $referenceFrame;
          $bestFrames{$reference}{$query}{'query'} = $queryFrame;
          last;
        }
      }
      last;
    }
  }

  return grep{defined $bestFrames{$_->{'REFERENCE'}}{$_->{'QUERY'}} && $_->{'REFERENCE_FRAME'} == $bestFrames{$_->{'REFERENCE'}}{$_->{'QUERY'}}{'reference'} && $_->{'QUERY_FRAME'} == $bestFrames{$_->{'REFERENCE'}}{$_->{'QUERY'}}{'query'}} @$alignments;
}

=head2 cleanOverlaps

=head2

=head3 Description

Detect totally overlaping alignments and return a list of alignments that don't overlap totally each other.

=head3 Arguments

=over 4

=item

A reference array containing Mummer alignements.

=item

A reference hash containing reference sequences.

=item

A reference hash containing query sequences.

=back

=head3 Returns

=over 4

=item

Array reference containing a list of alignments that don't overlap totally each other.

=back

=cut

sub cleanOverlaps {
  my ($self, $alignments, $ref_seq, $qry_seq) = @_;
  my @cleanedAlignments;
  my %framedSequences;

  foreach my $alignment (sort {$b->{'END_ALIGNMENT_QUERY'} - $b->{'BEGIN_ALIGNMENT_QUERY'} <=> $a->{'END_ALIGNMENT_QUERY'} - $a->{'BEGIN_ALIGNMENT_QUERY'}} @$alignments){
    my $frameKey = $alignment->{'REFERENCE_FRAME'}.$alignment->{'QUERY_FRAME'};
    my $alignKey = $alignment->{'REFERENCE'}.$alignment->{'QUERY'};
    if(! exists $framedSequences{$frameKey}{$alignKey}){
      my $seq = $ref_seq->retrieveFastaSequence($alignment->{'REFERENCE'});
      $framedSequences{$frameKey}{$alignKey}{'REFERENCE'} = $seq->{$alignment->{'REFERENCE'}};
      $seq = $qry_seq->retrieveFastaSequence($alignment->{'QUERY'});
      $framedSequences{$frameKey}{$alignKey}{'QUERY'} = $seq->{$alignment->{'QUERY'}};
    }
    my $resQuery = substr($framedSequences{$frameKey}{$alignKey}{'QUERY'}, $alignment->{'BEGIN_ALIGNMENT_QUERY'} - 1, ($alignment->{'END_ALIGNMENT_QUERY'} - $alignment->{'BEGIN_ALIGNMENT_QUERY'} +1), '-' x ($alignment->{'END_ALIGNMENT_QUERY'} - $alignment->{'BEGIN_ALIGNMENT_QUERY'} +1));
    my $resReference = substr($framedSequences{$frameKey}{$alignKey}{'REFERENCE'}, $alignment->{'BEGIN_ALIGNMENT_REFERENCE'} - 1, ($alignment->{'END_ALIGNMENT_REFERENCE'} - $alignment->{'BEGIN_ALIGNMENT_REFERENCE'} +1), '-' x ($alignment->{'END_ALIGNMENT_REFERENCE'} - $alignment->{'BEGIN_ALIGNMENT_REFERENCE'} +1));
    if($resQuery !~ /^-+$/ && $resReference !~ /^-+$/){
      push(@cleanedAlignments, $alignment);
    }
  }
  return @cleanedAlignments;
}

=head2 cleanPositions

=head2

=head3 Description

Detect non linear position on query and reference for alignments

Exemple 1 :

align 1 => beginReference : 60, endReference : 80, beginQuery : 1, endQuery : 20

align 2 => beginReference : 90, endReference : 100, beginQuery : 30, endQuery : 40

align 3 => beginReference : 10, endReference : 20, beginQuery : 50, endQuery : 60

The alignment 3 will be removed because it implies to back to begining of reference

=head3 Arguments

=over 4

=item

A reference array containing a list of Mummer alignements.

=back

=head3 Returns

=over 4

=item

Array reference containing cleaned alignments.

=back

=cut

sub cleanPositions {

  my ($self, $alignments) = @_;

  my @cleanedAlignments;
  my $lastQuery;
  my $lastReference;
  my ($lastBeginQuery, $lastEndQuery, $lastBeginReference, $lastEndReference);

  foreach my $alignment (sort {$a->{'REFERENCE'} cmp $b->{'REFERENCE'} || $a->{'QUERY'} cmp $b->{'QUERY'} || $a->{'BEGIN_ALIGNMENT_QUERY'} <=> $b->{'BEGIN_ALIGNMENT_QUERY'} || $a->{'BEGIN_ALIGNMENT_REFERENCE'} <=> $b->{'BEGIN_ALIGNMENT_REFERENCE'}} @$alignments){

    if(! defined $lastReference || ! defined $lastQuery || $lastReference ne $alignment->{'REFERENCE'} || $lastQuery ne $alignment->{'QUERY'}){

      $lastReference = $alignment->{'REFERENCE'};
      $lastQuery = $alignment->{'QUERY'};
      ($lastBeginQuery, $lastEndQuery, $lastBeginReference, $lastEndReference) = (0, 0, 0, 0);
    }

    #print $lastReference. ' - ' . $lastQuery . ' : ' .'('.$alignment->{'BEGIN_ALIGNMENT_REFERENCE'} . ', ' . $alignment->{'END_ALIGNMENT_REFERENCE'}.')'. "\n";

    if($alignment->{'BEGIN_ALIGNMENT_QUERY'} > $lastBeginQuery && $alignment->{'END_ALIGNMENT_QUERY'} > $lastEndQuery && $alignment->{'BEGIN_ALIGNMENT_REFERENCE'} > $lastBeginReference && $alignment->{'END_ALIGNMENT_REFERENCE'} > $lastEndReference){

      push(@cleanedAlignments, $alignment);
      ($lastBeginQuery, $lastEndQuery, $lastBeginReference, $lastEndReference) = ($alignment->{'BEGIN_ALIGNMENT_QUERY'}, $alignment->{'END_ALIGNMENT_QUERY'}, $alignment->{'BEGIN_ALIGNMENT_REFERENCE'}, $alignment->{'END_ALIGNMENT_REFERENCE'});
    }
  }

  return @cleanedAlignments;
}

=head2 getGapTableReference

=head2

=head3 Description

Generate a hash containing gaps positions on gapped reference sequence.

=head3 Arguments

=over 4

=item

A reference id

=item

A reference array containing a list of Mummer alignments

=back

=head3 Returns

=over 4

=item

A hash containing gap positions on gapped reference as keys and number of gaps to insert at these position as value.

=back

=cut

sub getGapTableReference {

  my ($self, $reference, $alignments) = @_;
  my %gapTable;
  my $lastQuery = '';
  my $lastQueryBegin = 0;
  my $lastQueryEnd = 0;
  my $lastReference = '';
  my $lastReferenceBegin = 0;
  my $lastReferenceEnd = 0;

  foreach my $alignment ( sort {$a->{'REFERENCE'} cmp $b->{'REFERENCE'} || $a->{'QUERY'} cmp $b->{'QUERY'} || $a->{'BEGIN_ALIGNMENT_REFERENCE'} <=> $b->{'BEGIN_ALIGNMENT_REFERENCE'}} @{$alignments} ){

    if($alignment->{'REFERENCE'} eq $reference){

      if($alignment->{'QUERY'} ne $lastQuery){
        $lastQuery = $alignment->{'QUERY'};
        $lastQueryBegin = 0;
        $lastQueryEnd = 0;
        $lastReferenceBegin = 0;
        $lastReferenceEnd = 0;
      }

      my ($beginReference, $alignmentCorrectionReference, $beginQuery, $alignmentCorrectionQuery) = $self -> gapBeforeAlignment($alignment->{'BEGIN_ALIGNMENT_REFERENCE'}, $lastReferenceEnd, $alignment->{'BEGIN_ALIGNMENT_QUERY'}, $lastQueryEnd);

      if ( $alignmentCorrectionReference && (! exists $gapTable{$lastReferenceEnd+1} || $gapTable{$lastReferenceEnd+1} < $alignmentCorrectionReference )){ $gapTable{$lastReferenceEnd +1 } = $alignmentCorrectionReference }

      my ($gaps) = getGapHashFromDeltaString($alignment->{'DELTA_ALIGNMENT'}, $alignment->{'BEGIN_ALIGNMENT_REFERENCE'}, $alignment->{'BEGIN_ALIGNMENT_QUERY'}, $beginReference, $beginQuery, undef, undef,$self->binary eq 'promer');
      my %gapsToAdd;

      foreach my $gap (keys %$gaps){
          $gapsToAdd{$gap} += $gaps -> {$gap};
      }
      foreach my $position (keys %gapsToAdd){
        if(! exists $gapTable{$position} || $gapTable{$position} < $gapsToAdd{$position}){
          $gapTable{$position} = $gapsToAdd{$position};
        }
      }
      $lastQueryEnd = $alignment->{'END_ALIGNMENT_QUERY'};
      $lastReferenceEnd = $alignment->{'END_ALIGNMENT_REFERENCE'};
      $lastReferenceBegin = $beginReference;
      $lastQueryBegin = $beginQuery;
    }
  }

  return %gapTable;
}

=head2 getGapTableQueries

=head2

=head3 Description

Generate a hash containing gaps positions on gapped query sequence.

=head3 Arguments

=over 4

=item

A reference id

=item

A reference array containing a list of Mummer alignments

=item

A hash reference containing reference sequences

=item

A hash reference containing query sequences

=item

A hash reference containing gap positions on gapped query

=back

=head3 Returns

=over 4

=item

A hash of hash containing gap positions on gapped queries as keys and number of gaps to insert at these position as value.

Exemple :

Query1 => { 10 => 3, 25 => 1 }

Query2 => { 14 => 2, 28 => 1 }

=back

=cut

sub getGapTableQueries {
  my ($self, $reference, $alignments, $ref_seq, $qry_seq, $gapTableReference) = @_;
  my %gapTable;
  my %gapReference;
  my @queryList;
  my $lastQuery = '';
  my $lastQueryBegin = 0;
  my $lastQueryEnd = 0;
  my $lastReference = '';
  my $lastReferenceBegin = 0;
  my $lastReferenceEnd = 0;

  foreach my $alignment ( sort {$a->{'REFERENCE'} cmp $b->{'REFERENCE'} || $a->{'QUERY'} cmp $b->{'QUERY'} || $a->{'BEGIN_ALIGNMENT_REFERENCE'} <=> $b->{'BEGIN_ALIGNMENT_REFERENCE'}} @{$alignments} ){
    if($alignment->{'REFERENCE'} eq $reference){
      if($alignment->{'QUERY'} ne $lastQuery){
        push (@queryList, $alignment->{'QUERY'});
        $lastQuery = $alignment->{'QUERY'};
        $lastQueryBegin = 0;
        $lastQueryEnd = 0;
        $lastReferenceBegin = 0;
        $lastReferenceEnd = 0;
      }
      my ($beginReference, $alignmentCorrectionReference, $beginQuery, $alignmentCorrectionQuery) = $self -> gapBeforeAlignment($alignment->{'BEGIN_ALIGNMENT_REFERENCE'}, $lastReferenceEnd, $alignment->{'BEGIN_ALIGNMENT_QUERY'}, $lastQueryEnd);
      my ($gapsToAddReference, $gapsToAddQuery) = getGapHashFromDeltaString($alignment->{'DELTA_ALIGNMENT'}, $alignment->{'BEGIN_ALIGNMENT_REFERENCE'}, $alignment->{'BEGIN_ALIGNMENT_QUERY'}, $beginReference, $beginQuery, undef, undef, $self->binary eq 'promer');

      foreach my $position (keys %$gapsToAddQuery){
        $gapTable{$alignment->{'QUERY'}}{$position} += $gapsToAddQuery -> {$position};
      }

      foreach my $position (keys %$gapsToAddReference){
        $gapReference{$alignment->{'QUERY'}}{$position} += $gapsToAddReference -> {$position};
      }

      if($alignmentCorrectionQuery){
        $gapTable{$alignment->{'QUERY'}}{$lastQueryEnd+1} += $alignmentCorrectionQuery;
      }

      if($alignmentCorrectionReference){
        $gapReference{$alignment->{'QUERY'}}{$lastReferenceEnd+1} += $alignmentCorrectionReference;
      }

      $lastQueryEnd = $alignment->{'END_ALIGNMENT_QUERY'};
      $lastReferenceEnd = $alignment->{'END_ALIGNMENT_REFERENCE'};
      $lastReferenceBegin = $beginReference;
      $lastQueryBegin = $beginQuery;
    }
  }

  foreach my $query (@queryList){
    my %gapDiffReference = $self->differenceGapTableHash($gapTableReference, $gapReference{$query});
    my %gapsToAdd;
    my $positionReference = 0;
    my $positionQuery = 0;
    my $stopReference = 0;
    my $stopQuery = 0;
    my $lastRef = 0;
    my $queryLength = 0;
    if(defined($ref_seq->{index}->{$query})){
      my $ref_seq_hash = $ref_seq->retrieveFastaSequence($query);
      $queryLength = length($ref_seq_hash->{$query});
    }
    if(defined($qry_seq->{index}->{$query})){
      my $qry_seq_hash = $qry_seq->retrieveFastaSequence($query);
      $queryLength = length($qry_seq_hash->{$query});
    }

    while(keys %gapDiffReference){
      if ($stopQuery > 0) { $stopQuery -- }
      if ($stopReference > 0) { $stopReference -- }
      if(! $stopQuery){
        $positionQuery ++;
        if( $gapTable{$query}{$positionQuery} ){ $stopQuery += $gapTable{$query}{$positionQuery} + 1 }
      }
      if(! $stopReference){
        $positionReference ++;
        if( $gapReference{$query}{$positionReference} ){ $stopReference += $gapReference{$query}{$positionReference} + 1 }
      }
      if($lastRef != $positionReference){
        if( $gapDiffReference{$positionReference} && $positionQuery < $queryLength){ $gapsToAdd{$positionQuery} += $gapDiffReference{$positionReference} }
        delete ($gapDiffReference{$positionReference});
        $lastRef = $positionReference;
      }
    }

    foreach my $position (keys %gapsToAdd){
      $gapTable{$query}{$position} += $gapsToAdd{$position} ;
    }
  }

  return %gapTable;
}

=head2 gapBeforeAlignment

=head2

=head3 Description

Calculate the number of gaps to insert between 2 Mummer alignments in order to correctly align query on the reference

=head3 Arguments

=over 4

=item

Begin position of alignment on the reference

=item

End position of the previous alignment on the reference

=item

Begin position of alignment on the query

=item

End position of the previous alignment on the query

=back

=head3 Returns

=over 4

=item

The new begin position of the alignment on the gapped reference

=item

The number of gap to insert on the reference

=item

The new begin position of the alignment on the gapped query

=item

The number of gaps to insert on the query

=back

=cut

sub gapBeforeAlignment {

  my ($self, $beginReference, $lastReferenceEnd, $beginQuery, $lastQueryEnd) = @_;
  my ($gapReference, $gapQuery) = (0,0);

  if($beginReference <= $lastReferenceEnd){
    $gapReference += $lastReferenceEnd - $beginReference + 1;
    $beginReference += $gapReference;
  }
  if($beginQuery <= $lastQueryEnd){
    $gapQuery += $lastQueryEnd - $beginQuery + 1;
    $beginQuery += $gapQuery;
  }

  my $gapDel = $gapReference - $gapQuery;

  if($gapDel <0){
    $gapReference = 0;
    $gapQuery = -1 * $gapDel;
  }
  else{
    $gapReference = $gapDel;
    $gapQuery = 0;
  }
  if( ($beginReference - $lastReferenceEnd) - ($beginQuery - $lastQueryEnd) > 0){
    $gapQuery += ($beginReference - $lastReferenceEnd) - ($beginQuery - $lastQueryEnd) ;
  }
  if( ($beginQuery - $lastQueryEnd) - ($beginReference - $lastReferenceEnd) > 0){
    $gapReference += ($beginQuery - $lastQueryEnd) - ($beginReference - $lastReferenceEnd);
  }

  return ($beginReference, $gapReference, $beginQuery, $gapQuery);
}

=head2 getGapHashFromDeltaString

=head2

=head3 Description

Retrieve gap positions from a Mummer delta string

=head3 Arguments

=over 4

=item

Mummer delta string

=item

Begin position of alignment on the reference

=item

Begin position of alignment on the query

=item

Minimum position to start retrieving gaps on reference

=item

Minimum position to start retrieving gaps on query

=item

Maximum position to retrieve gaps on reference

=item

Maximum position to retrieve gaps on query

=item

Boolean to indicate if the Mummer aligment is nucleic (undef or 0) or proteic (1)

=back

=head3 Returns

=over 4

=item

A an array containing :

- A hash reference containing the gap positions corresponding the Mummer alignment on reference

- A hash reference containing the gap positions corresponding the Mummer alignment on query

=back

=cut

sub getGapHashFromDeltaString {

  my ($deltaString, $lastNucleotideReference, $lastNucleotideQuery, $minReference, $minQuery, $maxReference, $maxQuery, $promer) = @_;

  if (! $lastNucleotideQuery){ $lastNucleotideQuery = 1 }
  if (! $lastNucleotideReference){ $lastNucleotideReference = 1 }

  my %gapsReference;
  my %gapsQuery;
  my $positionQuery = $lastNucleotideQuery;
  my $positionReference = $lastNucleotideReference;
  my $gapsToAdd = 1;

  if($promer){
    $gapsToAdd = 3;
  }

  foreach my $d (@$deltaString){

    my $distance = $d;

    if(defined $distance && $distance ne ''){
      if($promer){
        $distance = $distance * $gapsToAdd;
      }
      my %gapInfo = ();

      if ($distance < 0){
        $positionReference += (-1 * $distance) - $gapsToAdd;
        $positionQuery += (-1 * $distance);

        if( (! $minReference || $positionReference >= $minReference) && (! $minQuery || $positionQuery >= $minQuery) && (! $maxReference || $positionReference <= $maxReference) && (! $maxQuery || $positionQuery <= $maxQuery)){
          $gapsReference{$positionReference} += $gapsToAdd;
        }
      }
      elsif($distance > 0){
        $positionQuery += $distance - $gapsToAdd;
        $positionReference += $distance;

        if( (! $minReference || $positionReference >= $minReference) && (! $minQuery || $positionQuery >= $minQuery) && (! $maxReference || $positionReference <= $maxReference) && (! $maxQuery || $positionQuery <= $maxQuery)){
          $gapsQuery{$positionQuery} += $gapsToAdd;
        }
      }
    }
  }

  return (\%gapsReference, \%gapsQuery);
}

=head2 differenceGapTableHash

=head2

=head3 Description

Compare two hash of gap positions and retrieve non commons gaps between the two hashes

=head3 Arguments

=over 4

=item

The first hash of gap positions

=item

The second hash of gap positions

=item

Begin position of alignment on the query

=item

Minimum position to start retrieving non common gaps

=item

Maximum position to retrieve non common gaps

=back

=head3 Returns

=over 4

=item

A hash containing all the non common gaps between the two hash of gaps

=back

=cut

sub differenceGapTableHash {

  my ($self, $gapHash1, $gapHash2, $minValue, $maxValue) = @_;
  my %gapTable;

  foreach my $position (keys %$gapHash1){
    if( ( ! defined $minValue || $position >= $minValue) && ( ! defined $maxValue || $position <= $maxValue) ){ $gapTable{$position} = $gapHash1->{$position}; }
  }

  foreach my $position (keys %$gapHash2){
    if( ( ! defined $minValue || $position >= $minValue) && ( ! defined $maxValue || $position <= $maxValue) ){
      if (! exists $gapTable{$position}) { $gapTable{$position} = $gapHash2->{$position} }
      else{ $gapTable{$position} = abs($gapTable{$position} - $gapHash2->{$position}) }
    }
  }

  return %gapTable;
}

=head2 generateMultiFastaAlignment

=head2

=head3 Description

Generate a multi FASTA alignement file for each reference sequence.

=head3 Arguments

=over 4

=item

Output directory where to generate multi FASTA files (optional)

Default is current directory

=back

=head3 Returns

=over 4

=item

An array containing the list of multi FASTA files

=back

=cut

sub generateMultiFastaAlignment {
  my ($self, $outputDirectory) = @_;
  my $lastRef;
  my %alignments;
  my %gapTableReference;
  my %gapTableQueries;
  my @fileList;
  my @deltaAlignments = $self->getAlignmentsFromDeltaFile($self->delta);
  my $ref_seq = Tools::Fasta->new(file => $self->reference);
  my $qry_seq = Tools::Fasta->new(file => $self->query);
  # my %sequences = $self->getSequences;
  #print 'Before filtering : ' . scalar(@deltaAlignments)."\n";

  @deltaAlignments = $self->cleanFrames(\@deltaAlignments);
  #print 'After cleaning frames : ' . scalar(@deltaAlignments)."\n";

  @deltaAlignments = $self->cleanOverlaps(\@deltaAlignments, $ref_seq, $qry_seq);
  #print 'After cleaning overlaping alignments : ' . scalar(@deltaAlignments)."\n";

  @deltaAlignments = $self->cleanPositions(\@deltaAlignments);
  #print 'After cleaning no linear positions : ' . scalar(@deltaAlignments)."\n";

  if(!defined $outputDirectory){
    $outputDirectory = './';
  }
  $outputDirectory = abs_path($outputDirectory);

  if(! -e $outputDirectory){
    make_path($outputDirectory);
  }
  my @order;
  foreach my $alignment (sort {$a->{'REFERENCE'} cmp $b->{'REFERENCE'} || $a->{'QUERY'} cmp $b->{'QUERY'} || $a->{'BEGIN_ALIGNMENT_REFERENCE'} <=> $b->{'BEGIN_ALIGNMENT_REFERENCE'} } @deltaAlignments){
    my $seq;
    if(! exists $alignments{$alignment->{'REFERENCE'}}){
      if(%alignments){
        $self->{alnLengthMax} = 0 ;
        foreach my $id (keys %alignments){
          if(length($alignments{$id}) > $self->{alnLengthMax}){
            $self->{alnLengthMax} = length($alignments{$id});
          }
        }
        open(FILE, '> '.$outputDirectory.'/'.$lastRef.'.fasta');
        push (@fileList, $outputDirectory.'/'.$lastRef.'.fasta');

        foreach my $query (@order){
          print FILE '>' . $query."\n";
          if($self->{alnLengthMax} > length($alignments{$query})){
            my $nb = $self->{alnLengthMax} - length($alignments{$query}) ;
            print FILE $alignments{$query} . '-'x$nb ."\n";
          }
          else{
            print FILE $alignments{$query}."\n";
          }
        }
        close FILE;
      }
      %alignments = ();
      @order = ($alignment->{'REFERENCE'});
      %gapTableReference = $self->getGapTableReference($alignment->{'REFERENCE'}, \@deltaAlignments);
      $seq = $ref_seq->retrieveFastaSequence($alignment->{'REFERENCE'});
      $alignments{$alignment->{'REFERENCE'}} = $seq->{$alignment->{'REFERENCE'}};
      if ($alignment->{'REFERENCE_FRAME'} < 0) {
        $alignments{$alignment->{'REFERENCE'}} = reverse scalar $alignments{$alignment->{'REFERENCE'}};
        $alignments{$alignment->{'REFERENCE'}} =~ tr/ACGT/TGCA/;
      }
      my $nbGapInserted = 0;
      foreach my $gapPosition (sort { $a <=> $b } keys %gapTableReference){
        substr ( $alignments{$alignment->{'REFERENCE'}}, $gapPosition + $nbGapInserted - 1, 0, '-' x $gapTableReference{$gapPosition} );
        $nbGapInserted += $gapTableReference{$gapPosition};
      }
      %gapTableQueries = $self->getGapTableQueries($alignment->{'REFERENCE'}, \@deltaAlignments, $ref_seq, $qry_seq, \%gapTableReference);
      $lastRef = $alignment->{'REFERENCE'};
    }
    if( ! exists $alignments{$alignment->{'QUERY'}} ) {
      push (@order, $alignment->{'QUERY'});
      $seq = $qry_seq->retrieveFastaSequence($alignment->{'QUERY'});
      $alignments{$alignment->{'QUERY'}} = $seq->{$alignment->{'QUERY'}};
      if ($alignment->{'QUERY_FRAME'} < 0) {
        $alignments{$alignment->{'QUERY'}} = reverse scalar $alignments{$alignment->{'QUERY'}};
        $alignments{$alignment->{'QUERY'}} =~ tr/ACGT/TGCA/;
      }
      my $nbGapInserted = 0;
      foreach my $gapPosition (sort { $a <=> $b } keys %{$gapTableQueries{$alignment->{'QUERY'}}}){
        substr ( $alignments{$alignment->{'QUERY'}}, $gapPosition + $nbGapInserted - 1, 0, '-' x $gapTableQueries{$alignment->{'QUERY'}}{$gapPosition} );
        $nbGapInserted += $gapTableQueries{$alignment->{'QUERY'}}{$gapPosition};
      }
    }
  }
  $self->{alnLengthMax} = 0 ;
  foreach my $id (keys %alignments){
    if(length($alignments{$id}) > $self->{alnLengthMax}){
      $self->{alnLengthMax} = length($alignments{$id});
    }
  }
  if($lastRef){
    open(FILE, '> '.$outputDirectory.'/'.$lastRef.'.fasta');
    push (@fileList, $outputDirectory.'/'.$lastRef.'.fasta');
    foreach my $query (@order){
      print FILE '>'.$query."\n";
      if($self->{alnLengthMax} > length($alignments{$query})){
        my $nb = $self->{alnLengthMax} - length($alignments{$query}) ;
        print FILE $alignments{$query} . '-'x$nb ."\n";
      }
      else{
        print FILE $alignments{$query}."\n";
      }
    }
    close FILE;
  }
  return @fileList;
}
