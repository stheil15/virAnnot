package Tools::Blast;

use strict;
use Bio::SeqIO;
use Data::Dumper;
use Logger::Logger;
use Bio::SearchIO;
use Statistics::Descriptive;

my %defaultParameters = (
'file'					=> undef,
'algo'					=> undef,
'format'				=> 'csv',
'max_hsp'				=> undef,
'max_hit'				=> undef,
'max_evalue'			=> undef,
'min_identity'			=> undef,
'min_score'				=> undef,
'min_query_overlap'		=> undef,
'min_hit_overlap'		=> undef,
'min_hsp_length'		=> undef,
'keep_first_hit_only'	=> 0,
'keep_first_hsp_only'	=> 0,
'remove_self_match'	=> 0,
'parse_description' => 0,
);

=head1 INDEXED BLAST RELATED METHODS

=head2

=head2 new

=head2

=head3 Description

Create a new Tools::Blast object

=head3 Arguments

=over 4

=item

A hash of parameters.

Currently accepted keys are :

file => blast file path

algo => blast program used (BLASTN, BLASTP, BLASTX, TBLASTX or TBLASTN)

format => blast format (CSV, M8, M9, TABLE, XML or M7)

max_hsp => maximum number of hsp to keep

max_hit => maximum number of hit to keep

max_evalue => evalue cutoff

min_identity => identity cutoff

min_score => score cutoff

min_query_overlap => query overlap cutoff

min_hit_overlap => hit overlap cutoff

min_hsp_length => hsp length cutoff

keep_first_hit_only => keep only 1 hit by query (0 = false, 1 = true)

keep_first_hsp_only => keep only 1 hsp by query (0 = false, 1 = true)

remove_self_match => remove self match (0 = false, 1 = true)

=back

=head3 Returns

=over 4

=item

A Tools::Blast object

=back

=cut

sub new {
	my ($class, %attrs) = @_;
	my $self = {};
	bless $self;
	my %processedParameters = _merge_parameters($self,$self->getDefaultParameters, \%attrs);
	# print Dumper %processedParameters;
	$self->file($processedParameters{'file'});
	$self->algo($processedParameters{'algo'});
	$self->format($processedParameters{'format'});
	$self->max_hsp($processedParameters{'max_hsp'});
	$self->max_hit($processedParameters{'max_hit'});
	$self->max_evalue($processedParameters{'max_evalue'});
	$self->min_identity($processedParameters{'min_identity'});
	$self->min_score($processedParameters{'min_score'});
	$self->min_query_overlap($processedParameters{'min_query_overlap'});
	$self->min_hit_overlap($processedParameters{'min_hit_overlap'});
	$self->min_hsp_length($processedParameters{'min_hsp_length'});
	$self->keep_first_hit_only($processedParameters{'keep_first_hit_only'});
	$self->keep_first_hsp_only($processedParameters{'keep_first_hsp_only'});
	$self->remove_self_match($processedParameters{'remove_self_match'});
	$self->parse_description($processedParameters{'parse_description'});

	return $self;
}


sub _merge_parameters {
	my ($self,$h_default,$h_param)=@_;
	my %merged;
	foreach my $key (keys(%{$h_default})){
		if(defined($h_param->{$key})){
			$merged{$key} = $h_param->{$key};
		}
		elsif(exists($h_default->{$key})){
			$merged{$key} = $h_default->{$key};
		}
		else{
			$logger->logdie('Parameters ' . $key . ' not handled.');
		}
	}
	return %merged;
}


sub parse_description {
	my ($self, $param) = @_;
	if(defined($param)){
		if($param == 1 || $param == 0){
			$self->{'parse_description'} = $param;
		}
		else{
			$logger->error('Value provided from option parse_description is: ' . $param . '. Must be 0 or 1.');
		}
	}
	else{
		return $self->{'parse_description'};
	}
}


=head2 file

=head2

=head3 Description

Get/set blast file path

=head3 Arguments

=over 4

=item

A blast file path (optional)

=back

=head3 Returns

=over 4

=item

Blast file path

=back

=cut

sub file {
	my ($self, $file) = @_;
	if(defined $file){
		if(! ref $file && -e $file && ! -z $file){
			$self->{'file'} = $file;
		}
		else{
			$logger->warn("Invalid or empty blast file provided. Must be a no empty existing file. No change done");
		}
	}
	return $self->{'file'}
}

=head2 algo

=head2

=head3 Description

Get/set blast algorithm

=head3 Arguments

=over 4

=item

A valid blast program name : BLASTN, BLASTP, BLASTX, TBLASTX or TBLASTN (optional)

=back

=head3 Returns

=over 4

=item

Blast algorithm

=back

=cut

sub algo {
	my ($self, $algo) = @_;
	my %validAlgo = (
	'BLASTN' 	=> 'BLASTN',
	'BLASTP' 	=> 'BLASTP',
	'BLASTX' 	=> 'BLASTX',
	'TBLASTX' 	=> 'TBLASTX',
	'TBLASTN' 	=> 'TBLASTN',
	'DIAMONDX'  => 'DIAMONDX',
	'DIAMONDP'  => 'DIAMONDP'
	);
	if(defined $algo){
		if(! ref $algo && exists $validAlgo{uc($algo)}){
			$self->{'_algo'} = $validAlgo{uc($algo)};
		}
		else{
			$logger->warn("Invalid blast algorithm provided. Must be " . join(', ', sort {$a cmp $b} keys %validAlgo) . " . No change done");
		}
	}
	return $self->{'_algo'}
}

=head2 format

=head2

=head3 Description

Get/set blast file format

=head3 Arguments

=over 4

=item

A valid blast file format : CSV, M8, M9, TABLE, XML M7, blasttable or blastxml (optional)

CSV, M8, M9, TABLE, XML M7 will be converted to blasttable or blastxml format using the following rules :

CSV => blasttable

M8 => blasttable

M9 => blasttable

TABLE => blasttable

XML => blastxml

M7 => blastxml

=back

=head3 Returns

=over 4

=item

Blast file format (blasttable or blastxml)

=back

=cut

sub format {
	my ($self, $format) = @_;
	my %validFormat = (
	'CSV' 	=> 'blasttable',
	'XML' 	=> 'blastxml',
	'blasttable' => 'blasttable',
	'blastxml' => 'blastxml',
	);
	if(defined $format){
		if(! ref $format && exists $validFormat{uc($format)}){
			$self->{'format'} = $validFormat{uc($format)};
		}
		else{
			$logger->warn("Invalid blast format provided. Must be " . join(', ', sort {$a cmp $b} keys %validFormat) . " . No change done");
		}
	}
	elsif(! defined $self->{'format'}){
		$self->{'format'} = 'blasttable';
	}
	return $self->{'format'}
}

=head2 max_hsp

=head2

=head3 Description

Get/set the maximum number of hsp's to keep

=head3 Arguments

=over 4

=item

Maximum number of hsp's to keep (optional)

=back

=head3 Returns

=over 4

=item

Maximum number of hsp's to keep

=back

=cut

sub max_hsp {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff =~ /[0-9]+/ && $cutoff >= 1){
			$self->{'max_hsp'} = $cutoff;
		}
		else{
			$logger->warn("Invalid hsp cutoff provided. Must be a positive integer >= 1. No change done");
		}
	}
	return $self->{'max_hsp'}
}

=head2 max_hit

=head2

=head3 Description

Get/set the maximum number of hits to keep

=head3 Arguments

=over 4

=item

Maximum number of hits to keep (optional)

=back

=head3 Returns

=over 4

=item

Maximum number of hits to keep

=back

=cut

sub max_hit {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff =~ /[0-9]+/ && $cutoff >= 1){
			$self->{'max_hit'} = $cutoff;
		}
		else{
			$logger->warn("Invalid hit cutoff provided. Must be a positive integer >= 1. No change done");
		}
	}
	return $self->{'max_hit'}
}

=head2 max_evalue

=head2

=head3 Description

Get/set evalue cutoff

=head3 Arguments

=over 4

=item

Evalue cutoff (optional)

=back

=head3 Returns

=over 4

=item

Evalue cutoff

=back

=cut

sub max_evalue {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff >= 0){
			$self->{'_max_evalue'} = $cutoff;
		}
		else{
			$logger->warn("Invalid evalue cutoff provided. Must be a positive number. No change done");
		}
	}
	return $self->{'_max_evalue'}
}

=head2 min_identity

=head2

=head3 Description

Get/set identity cutoff

=head3 Arguments

=over 4

=item

Identity cutoff (optional)

=back

=head3 Returns

=over 4

=item

Identity cutoff

=back

=cut

sub min_identity {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff >= 0){
			$self->{'_min_identity'} = $cutoff;
		}
		else{
			$logger->warn("Invalid identity cutoff provided. Must be a positive number. No change done");
		}
	}
	return $self->{'_min_identity'}
}

=head2 min_score

=head2

=head3 Description

Get/set score cutoff

=head3 Arguments

=over 4

=item

Score cutoff (optional)

=back

=head3 Returns

=over 4

=item

Score cutoff

=back

=cut

sub min_score {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff >= 0){
			$self->{'min_score'} = $cutoff;
		}
		else{
			$logger->warn("Invalid score cutoff provided. Must be a positive number. No change done");
		}
	}
	return $self->{'min_score'}
}

=head2 min_query_overlap

=head2

=head3 Description

Get/set query overlap cutoff

=head3 Arguments

=over 4

=item

Query overlap cutoff (optional)

=back

=head3 Returns

=over 4

=item

Query overlap cutoff

=back

=cut

sub min_query_overlap {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff >= 0){
			$self->{'min_query_overlap'} = $cutoff;
		}
		else{
			$logger->warn("Invalid query overlap cutoff provided. Must be a positive number. No change done");
		}
	}
	return $self->{'min_query_overlap'}
}

=head2 min_hit_overlap

=head2

=head3 Description

Get/set hit overlap cutoff

=head3 Arguments

=over 4

=item

Hit overlap cutoff (optional)

=back

=head3 Returns

=over 4

=item

Hit overlap cutoff

=back

=cut

sub min_hit_overlap {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff >= 0){
			$self->{'min_hit_overlap'} = $cutoff;
		}
		else{
			$logger->warn("Invalid hit overlap cutoff provided. Must be a positive number. No change done");
		}
	}
	return $self->{'min_hit_overlap'}
}

=head2 min_hsp_length

=head2

=head3 Description

Get/set hsp length cutoff

=head3 Arguments

=over 4

=item

Hsp length cutoff (optional)

=back

=head3 Returns

=over 4

=item

Hsp length cutoff

=back

=cut

sub min_hsp_length {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && $cutoff =~ /[0-9]+/ && $cutoff >= 1){
			$self->{'_min_hsp_length'} = $cutoff;
		}
		else{
			$logger->warn("Invalid hsp length cutoff provided. Must be a positive integer >= 1. No change done");
		}
	}
	return $self->{'_min_hsp_length'}
}

=head2 keep_first_hit_only

=head2

=head3 Description

Get/set keep_first_hit_only option

=head3 Arguments

=over 4

=item

Boolean to activate or desactivate the keep_first_hit_only option (optional)

=> 0 : keep all hits

=> 1 : keep only 1 hit by query

=back

=head3 Returns

=over 4

=item

keep_first_hit_only state (0 or 1)

=back

=cut

sub keep_first_hit_only {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && ($cutoff == 0 || $cutoff == 1)){
			$self->{'keep_first_hit_only'} = $cutoff;
		}
		else{
			$logger->warn("Invalid value for keep_first_hit_only. Must be 0 or 1. No change done");
		}
	}
	return $self->{'keep_first_hit_only'}
}

=head2 keep_first_hsp_only

=head2

=head3 Description

Get/set keep_first_hsp_only option

=head3 Arguments

=over 4

=item

Boolean to activate or desactivate the keep_first_hsp_only option (optional)

=> 0 : keep all hsp

=> 1 : keep only 1 hsp by hit

=back

=head3 Returns

=over 4

=item

keep_first_hsp_only state (0 or 1)

=back

=cut

sub keep_first_hsp_only {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && ($cutoff == 0 || $cutoff == 1)){
			$self->{'keep_first_hsp_only'} = $cutoff;
		}
		else{
			$logger->warn("Invalid for keep_first_hsp_only. Must be 0 or 1. No change done");
		}
	}
	return $self->{'keep_first_hsp_only'}
}

=head2 remove_self_match

=head2

=head3 Description

Get/set remove_self_match option

=head3 Arguments

=over 4

=item

Boolean to activate or desactivate the remove_self_match option (optional)

=> 0 : Keep self match

=> 1 : Remove self match

=back

=head3 Returns

=over 4

=item

remove_self_match state (0 or 1)

=back

=cut

sub remove_self_match {
	my ($self, $cutoff) = @_;
	if(defined $cutoff){
		if(! ref $cutoff && ($cutoff == 0 || $cutoff == 1)){
			$self->{'remove_self_match'} = $cutoff;
		}
		else{
			$logger->warn("Invalid for remove_self_match. Must be 0 or 1. No change done");
		}
	}
	return $self->{'remove_self_match'}
}

=head2 getDefaultParameters

=head2

=head3 Description

Get Blast default parameters for Tools::Blast

=head3 Arguments

=over 4

=item

None

=back

=head3 Returns

=over 4

=item

A hash containing the default parameters used for Tools::Blast

Current default options = (

'file' => undef,

'algo' => undef,

'format' => 'csv',

'max_hsp' => undef,

'max_hit' => undef,

'max_evalue' => undef,

'min_identity' => undef,

'min_score' => undef,

'min_query_overlap' => undef,

'min_hit_overlap' => undef,

'min_hsp_length' => undef,

'keep_first_hit_only' => 0,

'keep_first_hsp_only' => 0,

'remove_self_match' => 0,

);

=back

=cut

sub getDefaultParameters {
	my ($class) = @_;
	return \%defaultParameters
}

=head2 setDefaultParameters

=head2

=head3 Description

Set Blast default parameters for Tools::Blast

=head3 Arguments

=over 4

=item

A hash reference containing the default parameters to use for Tools::Blast

Currently accepted keys are :

file => blast file path

algo => blast program used (BLASTN, BLASTP, BLASTX, TBLASTX or TBLASTN)

format => blast format (CSV, M8, M9, TABLE, XML or M7)

max_hsp => maximum number of hsp to keep

max_hit => maximum number of hit to keep

max_evalue => evalue cutoff

min_identity => identity cutoff

min_score => score cutoff

min_query_overlap => query overlap cutoff

min_hit_overlap => hit overlap cutoff

min_hsp_length => hsp length cutoff

keep_first_hit_only => keep only 1 hit by query (0 = false, 1 = true)

keep_first_hsp_only => keep only 1 hsp by query (0 = false, 1 = true)

remove_self_match => remove self match (0 = false, 1 = true)

=back

=head3 Returns

=over 4

=item

A hash containing the default parameters used for Tools::Blast

=back

=cut

sub setDefaultParameters {
	my ($class, $parameters) = @_;
	if($parameters){
		if(ref $parameters eq 'HASH'){
			%defaultParameters = %$parameters;
		}
		else{
			$logger->logdie('Parameters passed is not an HASH.');
		}
	}
	else{
		$logger->warn("You have not provided new default parameters. No change done");
	}
	return %defaultParameters
}

=head2 _readInputFile

=head2

=head3 Description

Parse a blast input file

=head3 Arguments

=over 4

=item

A Tools:Taxonomy object used to retrieve taxonomic information assotiated to the hits (optional)

=item

A TSV file containing the number of reads per contig (optional)

Format is :

Contig1 <TAB> ReadsNumer1

Contig2 <TAB> ReadsNumer2

...

=item

A FASTA file containing the query sequences (optional)

=back

=head3 Returns

=over 4

=item

An array of matches (hash reference containing informations on blast hits and hsp's)

=back

=cut

sub _readInputFile {
	my ($self, $taxonomyTools, $readsPerContigFile, $sequencesFile) = @_;
	$logger->debug('Reading blast file: ' . $self->file);
	$self->{_queryCounter} = 0;
	$self->{_searchioObj} = Bio::SearchIO->new(   '-file'   => $self->file,
	'-format' => $self->format);

	my @matches;
	my $readsPerContig;
	my $sequences;

	while($self->{_result} = $self->{_searchioObj}->next_result()){
		#~ Getting the BLAST algo type. Used for the calculation of the coverage.
		my $algo = $self->{_result}->algorithm();
		if(defined $algo && ! $self->{'_algo'}){
			$self->algo($algo);
		}
		$self->{_queryCounter} += 1;
		my $oldQueryId='';
    	while ($self->{_currentHit} = $self->{_result}->next_hit() || ! $self->{_result}->num_hits()){
			my $match = {};

			if($oldQueryId ne $self->{_result}->query_name()){
				$logger->debug('Treating query ' . $self->{_result}->query_name());
				$oldQueryId = $self->{_result}->query_name();
			}
			if($self->{_result}->num_hits() == 0){
				$match->{tax_id} = 0 ;
				$match->{no_hit} = 1 ;
				$match->{taxonomy} = 'unknown';
				$oldQueryId = $self->{_result}->query_name();
			}
			if($self->{'parse_description'} == 1){
				$match->{query_id} = $self->{_result}->query_description();
			}
			else{
				$match->{query_id} = $self->{_result}->query_name();
			}

			if(scalar(split(" ",$match->{query_id}) > 1)){
				$match->{query_id} = $1 if $match->{query_id} =~ /(\S+)\s/;
			}
			$match->{query_length} = $self->{_result}->query_length();
			$match->{nb_hits}      = $self->{_result}->num_hits();
			$match->{'algo'}         = $self->{'_algo'};

			#~ Removing match against itself. Just in case.
			if(! $match->{no_hit} && $self->remove_self_match) {
				$self->{_currentHit}->name eq $match->{query_id} and next;
			}

			#~ Parsing all the hsps. They are stored in @{$match->{hsps}}
			if(! $match->{no_hit}){
				_parseHsp($self,$match);
			}
			#~ Retieving the taxonomy.
			if(! $match->{no_hit}){
				my @list_description = split(">",$self->{_currentHit}->description());
				if(scalar(@list_description) > 1){
					$match->{description} = $list_description[0];
					$match->{description} =~ s/\s$//;
				}
				else{
					$match->{description} = $self->{_currentHit}->description();
				}
				$match->{hit_id} = $self->{_currentHit}->name;
				# $match->{gi} = $self->_parseGi($match);
        		$match->{accession} = $self->{_currentHit}->accession;
				if($taxonomyTools){
					if(! defined($taxonomyTools->{gg})){
						$taxonomyTools->_getNCBITaxonomy($match);
					}
					else{
						$taxonomyTools->_getGreenGenesTaxonomy($match);
					}
				}
				$self->{_hitAboveCutoff}++;
			}else{
        		$match->{'taxonomy'} = 'unknown'
      		}
			push(@matches,$match);
			if($self->keep_first_hit_only || $match->{no_hit}){
				$logger->trace('FirstHit : last.');
				last; # print only the first HSP
			}
		}
	}
	$logger->debug('Done.');
	return \@matches;
}

=head2 _parseGi

=head2

=head3 Description

Parse hit id to extract gi number

=head3 Arguments

=over 4

=item

A hash reference containing hit informations (must contains at least a hit_id key).

=back

=head3 Returns

=over 4

=item

A gi number if the id was successfully parsed or undef in case of error

=back

=cut

sub _parseGi {
	my ($self,$m) = @_;
	if($m->{hit_id} =~ /^gi\|(\d+)\|.*$/){
		return $1;
	}
	elsif($m->{hit_id} =~ /^lcl\|(\d+)$/){
		return $1;
	}
	elsif($m->{hit_id} =~ /^(\d+)$/){
		return $1;
	}
	else{
		# $logger->error('gi not found for ' . $m->{hit_id});
	}
}


sub _parseAcc {
	my ($self,$m) = @_;
	if($m->{hit_id} =~ /^gi\|\d+\|\S+\|(\S+)\|$/){
		return $1;
	}
	elsif($m->{hit_id} =~ /^lcl\|(\d+)$/){
		return $1;
	}
	elsif($m->{hit_id} =~ /^(\d+)$/){
		return $1;
	}
	else{
		# $logger->error('gi not found for ' . $m->{hit_id});
	}
}

=head2 _parseHsp

=head2

=head3 Description

Parse hit hsp's and add valid hsp's to hit

=head3 Arguments

=over 4

=item

A hash reference containing hit informations (must contains at least a hit_id key).

=back

=head3 Returns

=over 4

=item

1 if the hit contains at least one valid hsp that successfully passed the filters (evalue cutoff, identity cutoff ...)

=item

0 if the hit don't have valid hsp

=back

=cut

sub _parseHsp {
	my ($self,$match) = @_;
  $logger->debug('Parsing HSPs for ' . $self->{_result}->query_name());
	my @hsps = $self->{_currentHit}->hsps();
	$match->{hsps} = [];
	my $aboveCutoff = 0;

	# all HSPs in different strand than 1st HSPs will be discarded.
	$self->{_strand} = $hsps[0]->strand('query');

	#~ Initializing cumulative variables.
	for ('query','hit'){
		@{$self->{_alignCoords}->{$_}} = ();
		@{$self->{_alignCoords}->{$_}} = ();
		$self->{_cumulHspLength}->{$_} = 0;
		$self->{_cumulNumIdentical}->{$_} = 0;
		$self->{_cumulEvalue} = 0;
	}

	#~ Running through the hsps array.
	foreach ( @hsps ) {
		$self->{_hsp} = $_;
		#~ Test if an hsp is not on the same strand than the first one.
		$self->{_hsp}->strand('query') ne $self->{_strand} and next;

		my $hsp;
		$hsp->{hsp_length} = $self->{_hsp}->length('total');
		$hsp->{evalue} = $self->{_hsp}->evalue;
		$hsp->{num_identical} = $self->{_hsp}->num_identical();
		$hsp->{percentIdentity} = sprintf("%.1f", $self->{_hsp}->percent_identity());
		$hsp->{score} = $self->{_hsp}->bits;
		$hsp->{strand} = ($self->{_hsp}->strand('hit') == 0 ? $self->{_hsp}->strand('query') : $self->{_hsp}->strand('hit'));
		$hsp->{startQ} = $self->{_hsp}->start('query');
		$hsp->{endQ}   = $self->{_hsp}->end('query');
		$hsp->{startH} = $self->{_hsp}->start('hit');
		$hsp->{endH}   = $self->{_hsp}->end('hit');
		$hsp->{gaps}   = $self->{_hsp}->gaps('query');
		$hsp->{strandH} = $self->{_hsp}->strand('hit');
		$hsp->{strandQ} = $self->{_hsp}->strand('query');
		$self->getOverlapValues($hsp,'hsp',$match->{query_length});
		$match->{gaps} += $hsp->{gaps};
		my @identical = $self->{_hsp}->matches();
		$hsp->{mismatches} = ( $hsp->{hsp_length} - ($identical[0]) );
		$match->{mismatches} += $hsp->{mismatches};
		if(defined($self->{_max_evalue})){
			if( $hsp->{evalue} >= $self->{_max_evalue} ){
				$logger->trace($match->{query_id} . ' | ' . $self->{_currentHit}->name . ' under evalue threshold ' . $hsp->{evalue} . ' >= ' . $self->{_max_evalue});
				next;
			}
		}
		#~ All threshold tests for hsps.
		if(defined($self->{_min_identity})){
			if($hsp->{percentIdentity} < $self->{_min_identity}){
				$logger->trace($match->{query_id} . ' | ' . $self->{_currentHit}->name . ' under identity threshold ' . $hsp->{percentIdentity} . ' < ' . $self->{_min_identity});
				next;
			}
		}

		if($self->{_min_score}){
			if($hsp->{score} < $self->{_min_score}){
				$logger->trace($match->{query_id} . ' | ' . $self->{_currentHit}->name . ' under score threshold ' . $hsp->{score} . ' < ' . $self->{_min_score});
				next;
			}
		}
		if($self->format eq 'blasttable'){
			$hsp->{_numIdentical} = sprintf("%.0f", ($hsp->{hsp_length}-($hsp->{mismatches}-$hsp->{gaps})));
		}
		else{
			$hsp->{_numIdentical} = $self->{_hsp}->num_identical;
		}

		foreach (@{$self->{_alignCoords}->{'hit'}}){
			if ($hsp->{startH} >= $_->[0] and $hsp->{endH} <= $_->[1]) {
				next;
			}
		}
		push(@{$self->{_alignCoords}->{'hit'}}, [$hsp->{startH},$hsp->{endH}]);

		foreach (@{$self->{_alignCoords}->{'query'}}){
			if ($hsp->{startQ} >= $_->[0] and $hsp->{endQ} <= $_->[1]) {
				next;
			}
		}
		push(@{$self->{_alignCoords}->{'query'}}, [$hsp->{startQ},$hsp->{endQ}]);
		$self->{_cumulNumIdentical}->{'query'} += $hsp->{_numIdentical};
		$self->{_cumulHspLength}->{'query'}    += $hsp->{hsp_length};
		$self->{_cumulHspLength}->{'hit'}    += $hsp->{hsp_length};
		$self->{_cumulEvalue} += $hsp->{evalue};
		$match->{score} += $hsp->{score};
		$aboveCutoff+=1 ;
		push(@{$match->{hsps}},$hsp);
		last if ($self->{_keep_first_hsp_only}); # print only the first HSP
	} # END OF foreach hsp

	#~ Test if the hsp's array is not empty.
	if(scalar( @{$match->{hsps}} ) == 0){
		$match->{no_hit} = 1 ;
		return 0;
	}
	else{
		$match->{nb_hsps} = scalar( @{$match->{hsps}} );
	}
	$match->{percentIdentity} = sprintf("%.1f", $self->{_cumulNumIdentical}->{query}/$self->{_cumulHspLength}->{query}*100);
	$match->{evalue} = $self->{_cumulEvalue} / scalar(@{$match->{hsps}});
	$match->{hsp_length} = $self->{_cumulHspLength}->{query} ;
	# test cutoff for length and identity [ONLY FOR QUERY values !!!]
	if(defined($self->{_min_hsp_length})){
    if($match->{hsp_length} < $self->{_min_hsp_length}){
      $match->{no_hit} = 1 ;
      return 0;
    }
  }
	if(defined($self->{_min_identity})){
    if($match->{percentIdentity} < $self->{_min_identity}){
      $match->{no_hit} = 1 ;
      return 0;
    }
  }

	($match->{queryOverlap}, $match->{hitOverlap}) = ("NA.", "NA.");
	unless($self->format eq 'blasttable'){ # No query_length and hit_length available when inFormat is 'm8'
		$self->getOverlapValues($match, 'hit',$match->{query_length});
		if(defined($self->{'min_query_overlap'})){
      if($match->{queryOverlap} < $self->{'min_query_overlap'}){
        $match->{no_hit} = 1 ;
        $logger->debug('excluded min_query_overlap: ' . $match->{queryOverlap});
        return 0;
      }
    }
		if(defined($self->{'min_hit_overlap'})){
      if($match->{hitOverlap} < $self->{'min_hit_overlap'}){
        $match->{no_hit} = 1 ;
        $logger->debug('excluded min_hit_overlap: ' . $match->{queryOverlap});
        return 0;
      }
    }
	}
	$self->getMinMaxCoordsOfHsp($match);
	return 1;
} # END OF sub parseHsp

=head2 getMinMaxCoordsOfHsp

=head2

=head3 Description

Set start and end position of the hsp

=head3 Arguments

=over 4

=item

A hash reference containing hsp informations.

=back

=head3 Returns

=over 4

=item

1

=back

=cut

sub getMinMaxCoordsOfHsp {
	my ($self,$match) = @_;
	for my $qh ('query', 'hit'){
		my $statObjStart = Statistics::Descriptive::Full->new() ;
		my $statObjEnd = Statistics::Descriptive::Full->new() ;
		foreach (@{$self->{_alignCoords}->{$qh}}) {
			$statObjStart->add_data($_->[0]);
			$statObjEnd->add_data($_->[1]);
		}
		if($qh eq 'query'){
			$match->{startQ} = $statObjStart->min();
			$match->{endQ}   = $statObjEnd->max();
		}
		else{
			$match->{startH} = $statObjStart->min();
			$match->{endH}   = $statObjEnd->max();
		}
	}
	return 1;
} # END OF getMinMaxCoordsOfHsp

=head2 getOverlapValues

=head2

=head3 Description

Set hsp or hit overlap values for hit and query

=head3 Arguments

=over 4

=item

A hash reference containing hsp or hit informations.

=item

A scalar to indicate if the alignment is a hit or an hsp :

hsp = > the alignment is an hsp

other string => the alignment is a hit

=item

The query length (optional)

=back

=head3 Returns

=over 4

=item

1

=back

=cut

sub getOverlapValues {
	my ($self,$match,$t,$qLength) = @_;
	my ($qAlignLen, $hAlignLen)=(0, 0);
	if($t eq 'hsp'){
		$qAlignLen = $self->{_hsp}->length;
		$hAlignLen = $self->{_hsp}->length;
	}
	else{
		$qAlignLen = $self->{_cumulHspLength}->{query};
		$hAlignLen = $self->{_cumulHspLength}->{hit};
	}

	if ($self->{'_algo'} eq 'BLASTX') {
		if ($qLength){ $match->{queryOverlap} = sprintf "%.0f", (($qAlignLen*3/$qLength)*100) }
		if ($self->{_currentHit}->length){ $match->{hitOverlap} = sprintf "%.0f", (($hAlignLen/$self->{_currentHit}->length)*100) }
	}
	elsif ($self->{'_algo'} eq 'TBLASTN') {
		if ($qLength){ $match->{queryOverlap} = sprintf "%.0f", (($qAlignLen/$qLength)*100) }
		if ($self->{_currentHit}->length){ $match->{hitOverlap} = sprintf "%.0f", (($hAlignLen*3/$self->{_currentHit}->length)*100) }
	}
	elsif ($self->{'_algo'} eq 'TBLASTX') {
		if ($qLength){ $match->{queryOverlap} = sprintf "%.0f", (($qAlignLen*3/$qLength)*100) }
		if ($self->{_currentHit}->length){ $match->{hitOverlap} = sprintf "%.0f", (($hAlignLen*3/$self->{_currentHit}->length)*100) }
	}
	else {
		if ($qLength){ $match->{queryOverlap} = sprintf "%.0f", (($qAlignLen/$qLength)*100) }
		if ($self->{_currentHit}->length){ $match->{hitOverlap} = sprintf "%.0f", (($hAlignLen/$self->{_currentHit}->length)*100) }
	}

	$match->{queryOverlap} = 0 if(! $match->{queryOverlap});
	$match->{queryOverlap} = 100 if($match->{queryOverlap} > 100);
	$match->{hitOverlap} = 0 if(! $match->{hitOverlap});
	$match->{hitOverlap} = 100 if($match->{hitOverlap} > 100);

	return 1;
}

=head2 getBlastStats

=head2

=head3 Description

Generate taxonomic statistics on blast results

=head3 Arguments

=over 4

=item

An array reference containing a list of blast hits.

=back

=head3 Returns

=over 4

=item

Hash reference containing blast taxonomic statistics :

no_hit => number of no hit sequences,

no_viruses => number of annotated sequences with no viral taxid,

viruses => number of viral sequences,

dsRNA_viruses => number of dsRNA viral sequences,

ssRNA_viruses => number of ssRNA viral sequences,

dsDNA_viruses => number of dsDNA viral sequences,

ssDNA_viruses => number of ssDNA viral sequences,

undeterminated_viruses => 0,

viroids => number of viroids sequences,

bacteria => number of bacteria sequences,

archaea => number of archaea sequences,

cellular_organisms => number of cellular organisms sequences,

eukaryota => number of eukaryota sequences,

no_taxonomy => number of annotated sequences with no taxid,

taxonomy => number of annotated sequences with taxid,

other => number of annotated sequences with taxid corresponding to "other sequences" or "unclassified sequences",

annotated => number of annotaded sequences,

taxonomy => number of annotaded sequences with taxid

=back

=cut

sub getBlastStats {

	my ($class, $matches) = @_;

	my %blastStats = (no_hit => 0,
	no_viruses => 0,
	viruses => 0,
	dsRNA_viruses => 0,
	ssRNA_viruses => 0,
	dsDNA_viruses => 0,
	ssDNA_viruses => 0,
	undeterminated_viruses => 0,
	viroids => 0,
	bacteria => 0,
	archaea => 0,
	cellular_organisms => 0,
	eukaryota => 0,
	no_taxonomy => 0,
	taxonomy => 0,
	other => 0,
	annotated => 0,
	taxonomy => 0
	);

	foreach my $match (@{$matches}){
		if(! defined($match->{'no_hit'})){
			$blastStats{'annotated'} ++;
			if($match->{'taxonomy'} && $match->{'organism'} && $match->{'organism'} ne 'unknown' ){
				$blastStats{'taxonomy'} ++;
				if($match->{'taxonomy'} =~ /Viroids/){
					$blastStats{'viroids'} ++;
				}
				elsif($match->{'taxonomy'} =~ /other sequences/ || $match->{'taxonomy'} =~ /unclassified sequences/){
					$blastStats{'other'} ++;
				}
				elsif($match->{'taxonomy'} =~ /cellular organisms/){
					$blastStats{'cellular_organisms'} ++;
					if($match->{'taxonomy'} =~ /Bacteria/){
						$blastStats{'bacteria'} ++;
					}
					elsif($match->{'taxonomy'} =~ /Archaea/){
						$blastStats{'archaea'} ++;
					}
					elsif($match->{'taxonomy'} =~ /Eukaryota/){
						$blastStats{'eukaryota'} ++;
					}
					else{
						$blastStats{'undeterminated_cellular_organism'} ++;
					}
				}
				elsif($match->{'taxonomy'} =~ /Virus/){
					$blastStats{'viruses'} ++;
					if($match->{taxonomy} =~ /dsRNA viruses/){
						$blastStats{'dsRNA_viruses'} ++;
					}
					elsif($match->{taxonomy} =~ /ssRNA viruses/){
						$blastStats{'ssRNA_viruses'} ++;
					}
					elsif($match->{taxonomy} =~ /Retro-transcribing viruses|dsDNA viruses/){
						$blastStats{'dsDNA_viruses'} ++;
					}
					elsif($match->{taxonomy} =~ /ssDNA viruses/){
						$blastStats{'ssDNA_viruses'} ++;
					}
					else{
						$blastStats{'undeterminated_viruses'} ++;
					}
				}
			}
		}
		else{
			$blastStats{'no_hit'} ++;
		}
	}
	$blastStats{'no_viruses'} = $blastStats{'taxonomy'} - $blastStats{'viruses'};
	$blastStats{'no_taxonomy'} = $blastStats{'annotated'} - $blastStats{'taxonomy'};

	return \%blastStats;
}

1;
