#!/usr/bin/perl

use strict;
use warnings;

open (my $fh1, "<", $ARGV[0]) or die "Can't open fasta file: $ARGV[0]\n"; #File with 1 sequence in each line. all sequences should be in equal length   
open (my $fh2, ">", $ARGV[1]) or die "Can't open fasta file: $ARGV[1]\n"; #txt file of a matrix with the following dimensions:  the number of possible characters (e.g. a,c,g,t) X number of position in the with   
my $first_pos = $ARGV[2] ||= 1;

my @AoH; #array of hashes.
my $line_counter = 0; #count the number of sequences for the calculation of proportion
my $line = <$fh1>;
chomp $line;
my $length = length($line); #calculates the length of the sequences form the first line
seek($fh1,0,0); #goes back to the beginning of the file
while ($line = <$fh1>) { #reads through the sequences' file	
	chomp $line;
	$line_counter++;
	for my $i (0..$length-1) {
		my $nuc = uc(substr($line,0,1,"")); #assigns the first letter to a scalar and removes it from the sequence
		$AoH[$i]{$nuc}++;
	}
}
foreach my $nuc ("A","C","G","T") {
	for my $i (0..$length-1) {
		unless (defined ($AoH[$i]{$nuc})) {
			$AoH[$i]{$nuc} = 0; #initilize the AoH
		}
	}
}
# creating a uniq and sorted array with all the possible letters (e.g. nuclotides, amino acids):	
my @words=();
for my $i (0..$length-1) {
	push(@words,keys %{ $AoH[$i] });
}
my @unique;
my %seen;

foreach my $value (@words) {
  if (! $seen{$value}) {
    push @unique, $value;
    $seen{$value} = 1;
  }
}

@unique = sort(@unique);

# till here

print $fh2 "\t",join("\t", (@unique)),"\n"; #printing output matrix to file
for my $i (0..$length-1) {
	print $fh2 $first_pos+$i; 
	foreach my $keys (@unique) {
		if (not(defined($AoH[$i]{$keys}))) {
			$AoH[$i]{$keys}=0;
		}
		print $fh2 "\t".$AoH[$i]{$keys}/$line_counter;
	}
	print $fh2 "\n";
}

close $fh1;
close $fh2;
