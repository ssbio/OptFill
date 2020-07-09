#This function will remove redundant entires in an array, and return a reference to a non redundant array.
#inputs:
#	$_[0]=$@arrayRef = a reference to an array suspected of having redundancies
#outputs:
#	\@nr = reference to an array without those redundancies

use strict;

sub RemoveRedundancies {

	my $arrayRef = $_[0];
	my @array = @$arrayRef; #input array
	my %arrayKey = ( ); #allows us to check if a value exists in the array already. This is the inverse function of the non-redundant
	my @nr = ( ); #this is the non-redundant array

	for (my $i=0; $i<=$#array; $i++) {
	
		if (! exists $arrayKey{$array[$i]}) { #if the value does not have a key in the array
			
			$arrayKey{$array[$i]} = $#nr + 1; #add the value and the key to the hash, so it is known that this value is now in the array
			push @nr, $array[$i];
		
		} #otherwise do nothing, already in the array
	
	}
	
	return \@nr;

}

1;