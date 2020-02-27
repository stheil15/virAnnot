package Tools::Stat;

use strict;
use warnings;

use POSIX qw(ceil floor);

use base 'Exporter';
our @EXPORT = qw(getMinMax getNumericSortedValues getSum getMean getMedian getStandardDeviation getNandLStats);

=head1 STATISTIC RELATED METHODS 

=head2

=head2 getNumericSortedValues

=head2

=head3 Arguments

=over 4

=item

An array of text or numeric values.

=item

The numeric conversion to perform =

- length : the strings will be replaced by their length

- ascii : the strings will be splitted into characters and each character will be replaced by its ASCII code

=back

=head3 Returns

=over 4

=item

An array of numeric values sorted by numeric crescent order

=back

=cut

sub getNumericSortedValues{
    
        my ($valuesToConvert, $numericConversion) = @_;
        
        my @values = $valuesToConvert;
        
        if(! ref $valuesToConvert){ @values = ($valuesToConvert) }
        
        if ($numericConversion){
        
            if ($numericConversion eq 'length') { @values = map { length $_ } @$valuesToConvert }
            elsif ($numericConversion eq 'ascii') { @values = map { map { ord $_ } split('', $_) } @$valuesToConvert }
        }

        return sort { $a <=> $b } @values;
}

=head2 getMinMax

=head2

=head3 Arguments

=over 4

=item

An array of numeric sorted values.

=back

=head3 Returns

=over 4

=item

An array containing the minimum and the maximum value for the input values

=back

=cut

sub getMinMax{

        my ($numericSortedValues) = @_;
        
        if(! ref $numericSortedValues){ $numericSortedValues = [$numericSortedValues] }
        
        my ($minValue, $maxValue) = @$numericSortedValues[0, -1];

        if( $minValue > $maxValue ){ 
            
            ($minValue, $maxValue) = ($maxValue, $minValue)
            
        }        

        return ($minValue, $maxValue);
}

=head2 getSum

=head2

=head3 Arguments

=over 4

=item

An array of numeric values.

=back

=head3 Returns

=over 4

=item

The sum of the input values

=back

=cut

sub getSum{

        my ($numericValues) = @_;
        
        if(! ref $numericValues){ $numericValues = [$numericValues] }

        my $sum = 0;
        
        foreach my $value (@$numericValues){

            $sum += $value;                
        }

        return $sum;
}

=head2 getMean

=head2

=head3 Arguments

=over 4

=item

An array of numeric values.

=item

The sum of the input values (optional)

If the sum is provided, the fonction won't spend time recalculating it

=back

=head3 Returns

=over 4

=item

The mean of the input values

=back

=cut

sub getMean{

        my ($numericValues, $sum) = @_;

        if(! ref $numericValues){ $numericValues = [$numericValues] }

        my $mean = 0;
        
        if(@$numericValues){
            
            if (! defined $sum) { $sum = getSum($numericValues) }
            $mean = $sum / scalar @$numericValues;
        }

        return $mean;
}

=head2 getMedian

=head2

=head3 Arguments

=over 4

=item

An array of numeric values.

=item

The cutoff in percent (optional)

The default cutoff is 50% and correspond to the median

A cutoff of 25% will return the first quartile

A cutoff of 75% will return the third quartile

=back

=head3 Returns

=over 4

=item

By default the median of input values

In a more general case, the value where there is [cutoff]% of the elements that are lower than this value

=back

=cut

sub getMedian{

        my ($numericSortedValues, $cutoff) = @_;

        my $values = $numericSortedValues;
        
        if(! ref $numericSortedValues){ $values = [$numericSortedValues] }
        if(! defined $cutoff){ $cutoff = 50 }

        my $median;

        if(@$values){

            my $medianPosition;

            if ($values->[0] > $values->[@$values - 1]) { $values = [reverse @$values] }

            if ($cutoff == 50) {

                $medianPosition = (scalar @$values + 1) * $cutoff / 100;
                if  ( $medianPosition == int $medianPosition ){ $median = $values->[$medianPosition - 1] }
                else{$median = ($values->[floor($medianPosition) - 1] +  $values->[ceil($medianPosition) - 1])/2}
            }

            else{

                $medianPosition = (scalar @$values) * $cutoff / 100;
                if ($medianPosition == int $medianPosition ){ $median = $values->[$medianPosition - 1]}
                else {$median = $values->[ceil($medianPosition) - 1]}
            }

        }

        return $median;
}

=head2 getStandardDeviation

=head2

=head3 Arguments

=over 4

=item

An array of numeric values.

=item

The mean of the input values (optional)

If the mean is provided, the fonction won't spend time recalculating it

=back

=head3 Returns

=over 4

=item

The standard deviation of the input values

=back

=cut

sub getStandardDeviation{
	
        my ($numericValues, $mean) = @_;

        if(! ref $numericValues){ $numericValues = [$numericValues] }

        my $standardDeviation = 0;
        
        if(@$numericValues){

            if (! defined $mean) { $mean = getMean($numericValues) }

            my $numberOfValues = 0;

            foreach my $value (@$numericValues){

                $standardDeviation += (($mean - $value) ** 2);
                $numberOfValues ++;		
            }
	
            $standardDeviation = ($standardDeviation/$numberOfValues) ** 0.5;
        }
            
        return $standardDeviation;	
}

=head2 getNandLStats

=head2

=head3 Arguments

=over 4

=item

An array of numeric sorted values.

=item

The cutoff in percent (optional)

If no cutoff is provided, it will be set to 50

=item

The sum of the input values (optional)

If the sum is provided, the fonction won't spend time recalculating it

=back

=head3 Returns

=over 4

=item

An array containing to the N and L values corresponding to the cutoff.

Assuming that the sequence are sorted by decrescent order :

N = The length of the last sequence used to cover [cutoff] percent of the total sequences length

L = The number of sequences needed to cover [cutoff] percent of the total sequences length

=back

=cut

sub getNandLStats{
	
        my ($numericSortedValues, $cutoff, $sum) = @_;

        my $values = $numericSortedValues;
        
        if(! ref $numericSortedValues){ @$values = ($numericSortedValues) }
        if(! defined $cutoff){ $cutoff = 50 }

        my ($N, $L) = (0, 0);

        if(@$values){

            if ($values->[0] < $values->[@$values -1]) { $values = [reverse @$values] }

            if (! defined $sum) { $sum = getSum($numericSortedValues) }

            if ($sum){
                
                my $totalNucleotides = 0;

                foreach my $sequenceLength (@$values){

                    $totalNucleotides += $sequenceLength;
                    $N = $sequenceLength;
                    $L ++ ;
        
                    last if ( ($totalNucleotides*100) / $sum >= $cutoff );				
                }
            }
        }
    
        return ($N, $L);	
}

1;
