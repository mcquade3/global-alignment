#!/usr/local/bin/perl
# Mike McQuade
# Global-alignment.pl
# Finds the highest scoring alignment between
# two strings using a scoring matrix.

# Define the packages to use
use strict;
use warnings;
use List::Util qw(max);
use List::MoreUtils qw(firstidx);

# Initialize variables
my ($indelPenalty,$firstString,$secondString,@grid,@blosum62);

# Define the variables for indel penalty and blosum62
$indelPenalty = -5;

@blosum62 = (
	[4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2],
	[0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
	[-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3],
	[-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2],
	[-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3],
	[0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3],
	[-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2],
	[-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1],
	[-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2],
	[-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1],
	[-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1],
	[-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2],
	[-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3],
	[-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1],
	[-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2],
	[1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2],
	[0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2],
	[0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1],
	[-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2],
	[-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7]
);

# Open the file to read
open(my $fh,"<ba5e.txt") or die $!;

# Read in the values from the file
$firstString = <$fh>;
chomp($firstString);
$secondString = <$fh>;
chomp($secondString);

# Populate the first row of the grid with as many zeros
# as the length of the second string.
my @tempArr = (0);
for (my $i = 0; $i <= length($secondString); $i++) {
	push @tempArr, $tempArr[$i] - $indelPenalty;
}
push @grid,[@tempArr];

# Populate the first column of the grid with as many zeros
# as the length of the first string.
for (my $i = 1; $i <= length($firstString); $i++) {
	push @grid,[$indelPenalty*$i];
}

# Calculate the rest of the grid
for (my $i = 1; $i <= length($firstString); $i++) {
	for (my $j = 1; $j <= length($secondString); $j++) {
		my $firstChar = substr($firstString,$i-1,1);
		my $secondChar = substr($secondString,$j-1,1);
		$grid[$i][$j] = max(
								$grid[$i-1][$j] + $indelPenalty,
								$grid[$i][$j-1] + $indelPenalty,
								$grid[$i-1][$j-1] + matchVal($firstChar,$secondChar)
							);
	}
}

# Call the output function with the lengths of the
# given strings.
outputGlobalAlign();

# Close the file
close($fh) || die "Couldn't close file properly";


# Returns the penalty for a specific match of characters
sub matchVal {
	my $firstChar = $_[0];
	my $secondChar = $_[1];
	
	my @index = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
	my $firstIndex = firstidx { $_ eq $firstChar } @index;
	my $secondIndex = firstidx { $_ eq $secondChar } @index;

	return $blosum62[$firstIndex][$secondIndex];
}

# Print out the alignment of the two given strings
sub outputGlobalAlign {
	# Define local variables
	my $alignmentA = "";
	my $alignmentB = "";
	my $i = length($firstString);
	my $j = length($secondString);

	# Backtrack the path defined by the grid
	while ($i > 0 || $j > 0) {
		# If there is a match or mismatch, concatenate the last letter of each string
		# with its respective alignment.
		if ($i > 0 && $j > 0 && $grid[$i][$j] == $grid[$i-1][$j-1] + matchVal(substr($firstString,$i-1,1),substr($secondString,$j-1,1))) {
			$alignmentA = substr($firstString,$i-1,1).$alignmentA;
			$alignmentB = substr($secondString,$j-1,1).$alignmentB;
			$i--;
			$j--;
		# If there is a gap, one alignment string receives the last letter
		# of one of the strings while the other receives a dash.
		} elsif ($i > 0 && ($grid[$i][$j] == $grid[$i-1][$j] + $indelPenalty)) {
			$alignmentA = substr($firstString,$i-1,1).$alignmentA;
		    $alignmentB = "-".$alignmentB;
		    $i--;
		} else {
			$alignmentA = "-".$alignmentA;
		    $alignmentB = substr($secondString,$j-1,1).$alignmentB;
		    $j--;
		}
	}
	# Output the score and the alignments for the two strings
	print $grid[length($firstString)][length($secondString)]."\n";
	print $alignmentA."\n";
	print $alignmentB."\n";
}