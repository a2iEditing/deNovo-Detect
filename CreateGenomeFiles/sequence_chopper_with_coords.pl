#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper; 
# print Dumper %hash;

open (my $fh1, "<", $ARGV[0]) or die "Can't open $ARGV[0]\n"; #input fasta file with name consists of gene name, strand, chromosome, start-end seperated by space (" ")
my $window = $ARGV[1]//=152; #window width for sub sequences
my $jump = $ARGV[2]//=76; #window jumps. if < $window the windows will overlap one another 
open (my $fh2, ">", $ARGV[3]) or die "Can't open $ARGV[1]\n"; #output fasta file with sequences of sub sequences.

my @data; #array for lines splitted data
my $seq; #scalar for fasta sequence 
my %hash; #hash of arrays of arrays (HoAoA). key - element name. value - array which stores arrays for each of the parts that make the element. the subarrays includes the following fields: chr, start (0-based), end, sequence of sub-element
my $key;

while (my $line = <$fh1>) {
	chomp $line;
	if ($line =~ m/^>/) {
		my ($name) = $line =~ m/^>(.+)/;
		@data = split(" ",$name); #splits the fasta name to an array
		$key = $data[0]; #key is the sequence name
		$hash{$key}[0][0] = $data[1]; #sequence strand
		$hash{$key}[0][1] = $data[2]; #sequence chromosome
		for my $i (3..$#data) {
			my @intervals = split("-",$data[$i]); #splits start and end pairs for each subelement
			if ($hash{$key}[0][0] eq "-") {
				@intervals = reverse(@intervals); #reverse start and end coordinates if the sequence is from the negative strand
			}
			# $intervals[2] = $intervals[1] - $intervals[0]; #calculates the length of the intrval
			if ($i == 3) { 
				# $intervals[3] = 0; 
				$intervals[2] = 0; 
			} else {
				# $intervals[3] = abs($intervals[0] - $hash{$key}[$i-2][0]); #calculates the length of the gap between the current interval and the one before it 
				$intervals[2] = abs($intervals[0] - $hash{$key}[$i-3][1]); #calculates the length of the gap between the current interval and the one before it 
			}
			push(@{ $hash{$key} }, \@intervals ); #push the interval info arrays (as reference) to the appropriate key in the hash
		}
	} else {
		$hash{$key}[0][2].=$line #adds up the sequence in a specific array element
	}
}
# print Dumper %hash;
# exit; 

SEQ: foreach $key (keys(%hash)) {
	my $seq = $hash{$key}[0][2]; #gets the sequence from the hash
	my $seq_length = length($seq); #calculates the length of the sequence
	next SEQ if ($seq_length < $window); #skips the sequence if it is too short
	my $strand_flag=1; #sets flag to 1 to indicate positive strand as default
	if ($hash{$key}[0][0] eq "-") {
		$strand_flag=-1; #change flag to -1 if on the negative strand
	}
	print ">>>",$key,"\n";
	# print ">>>",$strand_flag,"\t",$hash{$key}[0][0],"\n";
	my $cur_start_coord = $hash{$key}[1][0] - ($jump*$strand_flag); #uses the first start coordinate as the start coordinate of the first subsequence
	# print ">>>",$cur_start_coord,"\n";
	my $start_index = 1; #the index in AoA to look at if the start coordinate is in
	my $end_index = 1; #the index in AoA to look at if the end coordinate is in
	my $cur_end_coord = $cur_start_coord+(1*$strand_flag); #initialize the end coordinate to the lowest value it can have whatsoever 
	my $cur_pos; #declared here for being available for next block
	SUB_SEQS: for ($cur_pos = 0; $cur_pos + $window <= $seq_length; $cur_pos += $jump) { # iterates over the sequence to extract sub sequences in length of $window and located in $jump from one another
		# $cur_start_coord += $cur_pos;
		# print "<<",$cur_pos,"\n";
		$cur_start_coord += ($jump*$strand_flag);
		# print ">>>",$cur_start_coord,"\n";
		# print ">>",$cur_start_coord,"\t",$hash{$key}[$start_index][1],"\n";
		# print "<<",$cur_pos,"\n";
		START_SEARCH: while (not (($cur_start_coord*$strand_flag >= $hash{$key}[$start_index][0]*$strand_flag) and ($cur_start_coord*$strand_flag < $hash{$key}[$start_index][1]*$strand_flag))) { #iterates if the start coordinate is not falling within the current interval
		# START_SEARCH: while (not $cur_start_coord < $hash{$key}[$start_index][1])) { #iterates if the start coordinate is not falling within the current interval
			$start_index++; #move to the next interval
			# $cur_start_coord+=$hash{$key}[$start_index][3]; #adds the gap length to the start coordinate 
			$cur_start_coord+=($hash{$key}[$start_index][2]*$strand_flag); #adds the gap length to the start coordinate 
			# print ">>>",$cur_start_coord,"\t",$start_index,"\n";
			# print ">>>",$cur_start_coord,"\t",$hash{$key}[$start_index][1],"\n";
		}
		my $gaps = 0; #sums the gaps between start coordinate and end coordinate
		if ($end_index <= $start_index) {
			$end_index = $start_index; #assign the index of the interval in which the end coord is located to be at least as small as the index of the strar coordinate
			$cur_end_coord = $cur_start_coord + ($window*$strand_flag); #in general the end is the start coordinate + window length unless it falls within a gap
		} else {
			for my $i (($start_index+1)..$end_index) {
				$gaps+=($hash{$key}[$i][2]);
			}
			$cur_end_coord = $cur_start_coord + (($window + $gaps)*$strand_flag); #in general the end is the start coordinate + window length unless it falls within a gap
		}
		# print "<<",$start_index,"\t",$end_index,"\n";
		END_SEARCH: while (not (($cur_end_coord*$strand_flag > $hash{$key}[$end_index][0]*$strand_flag) and ($cur_end_coord*$strand_flag <= $hash{$key}[$end_index][1]*$strand_flag))) { #iterates if the end coordinate is not falling within the current interval
		# END_SEARCH: while ($cur_end_coord <= $hash{$key}[$end_index][1]) { #iterates if the end coordinate is not falling within the current interval
			$end_index++; #move to the next interval
			# $cur_end_coord+=$hash{$key}[$end_index][3]; #adds the gap length to the end coordinate 
			$cur_end_coord+=($hash{$key}[$end_index][2]*$strand_flag); #adds the gap length to the end coordinate 
			# print ">>>",$cur_start_coord,"\t",$cur_end_coord,"\t",$end_index,"\n";
		}
		my $sub_seq = substr($seq,$cur_pos,$window); #extracts the subsequence from the original sequence
		my $coord_string=""; #includes the the different intervals that make the subsequence
		my @coord_sub = (); #Auxilary array of arrays
		if ($start_index == $end_index) { 
			# $coord_string = $cur_start_coord ."-". $cur_end_coord; #start and end are in the same interval
			push ( @coord_sub, [ $cur_start_coord, $cur_end_coord ]); #start and end are in the same interval
		} else { 
			# $coord_string = $cur_start_coord ."-". $hash{$key}[$start_index][1]; #from the start coordinate to the end of the "start interval"
			push ( @coord_sub, [ $cur_start_coord, $hash{$key}[$start_index][1] ]); #start and end are in the same interval #from the start coordinate to the end of the "start interval"
			if ($end_index - $start_index > 1) { #start and end are not in adjacent interval 
				for my $i (($start_index+1)..($end_index-1)) {
					# $coord_string .= "_". $hash{$key}[$i][0] ."-". $hash{$key}[$i][1]; #adds the intermediate intervals as a whole
					push ( @coord_sub, [ $hash{$key}[$i][0],$hash{$key}[$i][1] ]);  #adds the intermediate intervals as a whole
				}
			}
			# $coord_string .= "_". $hash{$key}[$end_index][0] ."-". $cur_end_coord; #from the start of the "end interval" to the end coordinate
			push ( @coord_sub, [ $hash{$key}[$end_index][0],$cur_end_coord ]); #from the start of the "end interval" to the end coordinate
		}
		for my $i (0..$#coord_sub) { 
			if ($hash{$key}[0][0] eq "-") {
				$coord_string .= "_" . join("-",reverse(@{ $coord_sub[$i]}[0..1]));
			} else {
				$coord_string .= "_" . join("-",@{ $coord_sub[$i]}[0..1]); 
			}
		}
		print $fh2 ">",join("_",($key,@{ $hash{$key}[0]}[0..1])),$coord_string,"\n"; #print the sequence name, strand, chromosome and associated intervals 
		PRINT_SEQ: while ($sub_seq ne "") {
			my $temp_seq = substr($sub_seq,0,60,"");
			print $fh2 $temp_seq,"\n"; #pulls out sub-strings of the sequence in order to print them in fasta format (aka - line length < 80) 
		}
	}
	# print "<<",$cur_pos,"\t",$window,"\t",$jump,"\t",$seq_length,"\n";
	if ($cur_pos + $window - $jump < $seq_length) { #deals with the last window that starts from the end of the sequence to avoid overhangs in the running window. jump is subtracted to compensate on the last $jump addition to SUB_SEQ loop
		$cur_end_coord = $hash{$key}[$#{ $hash{$key} }][1]; #uses the last end coordinate as the end coordinate of the last subsequence
		$end_index = $#{ $hash{$key} }; #the index in AoA to look at
		$start_index = $end_index; #starts to look for the interval containing the start coordinate form the last interval
		$cur_start_coord = $cur_end_coord - ($window*$strand_flag); #in general when looking from the back of the sequence the start is the end coordinate - window length unless it falls within a gap
		# START_BACK_SEARCH: while ($cur_start_coord < $hash{$key}[$start_index][0]) { #iterates if the start coordinate is not falling within the current interval
		START_BACK_SEARCH: while (not (($cur_start_coord*$strand_flag >= $hash{$key}[$start_index][0]*$strand_flag) and ($cur_start_coord*$strand_flag < $hash{$key}[$start_index][1]*$strand_flag))) { #iterates if the start coordinate is not falling within the current interval
			### $cur_start_coord+=$hash{$key}[$start_index][3]; #adds the gap length to the start coordinate 
			$cur_start_coord -= ($hash{$key}[$start_index][2]*$strand_flag); #adds the gap length to the start coordinate 
			$start_index--; #move to the next interval
		}
		my $sub_seq = substr($seq,-$window); #extracts the last subsequence from the end of original sequence
		my $coord_string=""; #includes the the different intervals that make the subsequence
		my @coord_sub = (); #Auxilary array of arrays
		if ($start_index == $end_index) { 
			# $coord_string = $cur_start_coord ."-". $cur_end_coord; #start and end are in the same interval
			push ( @coord_sub, [ $cur_start_coord, $cur_end_coord ]); #start and end are in the same interval
		} else { 
			# $coord_string = $cur_start_coord ."-". $hash{$key}[$start_index][1]; #from the start coordinate to the end of the "start interval"
			push ( @coord_sub, [ $cur_start_coord, $hash{$key}[$start_index][1] ]); #start and end are in the same interval #from the start coordinate to the end of the "start interval"
			if ($end_index - $start_index > 1) { #start and end are not in adjacent interval 
				for my $i (($start_index+1)..($end_index-1)) {
					# $coord_string .= "_". $hash{$key}[$i][0] ."-". $hash{$key}[$i][1]; #adds the intermediate intervals as a whole
					push ( @coord_sub, [ $hash{$key}[$i][0],$hash{$key}[$i][1] ]);  #adds the intermediate intervals as a whole
				}
			}
			# $coord_string .= "_". $hash{$key}[$end_index][0] ."-". $cur_end_coord; #from the start of the "end interval" to the end coordinate
			push ( @coord_sub, [ $hash{$key}[$end_index][0],$cur_end_coord ]); #from the start of the "end interval" to the end coordinate
		}
		for my $i (0..$#coord_sub) { 
			if ($hash{$key}[0][0] eq "-") {
				$coord_string .= "_" . join("-",reverse(@{ $coord_sub[$i]}[0..1]));
			} else {
				$coord_string .= "_" . join("-",@{ $coord_sub[$i]}[0..1]); 
			}
		}
		print $fh2 ">",join("_",($key,@{ $hash{$key}[0]}[0..1])),$coord_string,"\n"; #print the sequence name, strand, chromosome and associated intervals 
		PRINT_SEQ: while ($sub_seq ne "") {
			my $temp_seq = substr($sub_seq,0,60,"");
			print $fh2 $temp_seq,"\n"; #pulls out sub-strings of the sequence in order to print them in fasta format (aka - line length < 80) 
		} 
		PRINT_LAST_SEQ: while ($sub_seq ne "") {
			print $fh2 substr($sub_seq,0,60,""),"\n"; #pulls out sub-strings of the sequence in order to print them in fasta format (aka - line length < 80) 
		}
	}
}

close $fh1;
close $fh2;
