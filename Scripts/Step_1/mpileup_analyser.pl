#!/usr/bin/perl 
use strict;
use warnings; 
# this script gets mpileup output file and returns a bed-like file with counting of reads with each nucleotide\ with indels=

open (my $fh1, "<", $ARGV[0]) or die "Can't open input file: $ARGV[0]\n"; #input mpileup output file
open (my $fh2, ">", $ARGV[1]) or die "Can't create output file: $ARGV[1]\n";# bed-like output file
my $min_q=$ARGV[2]//=30; #minimum base quality score to include in results 
my $q_offset=$ARGV[3]//=33; #offset for ascii to base quality score conversion 

my @data;
# print $fh2 join("\t",("chr","start","end","ref","total_cov","a_plus","c_plus","g_plus","t_plus","n_plus","total_with_q_plus","read_start_plus","read_end_plus","del_plus","inser_plus","ref_skip_plus","low_q_plus","a","c","g","t","n","total_with_q","read_start","read_end","del","inser","ref_skip","low_q","del_placehold")),"\n"; #print header 
while (my $line=<$fh1>) {
	my %hash; #used to count the different characters\ nucleotides that found in the $seq string and represent reads covering the currnet position. the variable "$char" is used as a hash key and its content defines the categorization. upper case "$char" value represents read that was mapped to the positive strand of the reference sequence whereas lower case represents read that was mapped to the negative strand of the reference sequence. The different categories are: A, C,G,T or N nucleotide, TOTALQ: total nucleotides with quality >= $min_q, START: read start (^.), END: Read end ($), SKIP: reference skipping ([<>]) - N cigar, IN: insertion (\+[0-9]+[AGCTNacgtn]+), DELPH: deletion place holder (*), DEL: deletion in the upstream base(s) (\-[0-9]+[AGCTNacgtn]+) and LOW: low quality base - lower then $min_q
	chomp $line;
	@data = split("\t",$line); #splits each line to the fields
	my $chr = $data[0]; #position chromosome
	my $start = $data[1] - 1; #position start coordinate in 0-based format
	my $end = $data[1]; #position end coordinate
	my $ref = $data[2]; #reference nucleotide as annotated in the mpileup file. Equals "N" if reference nucleotide is unknown
	my $total = $data[3]; #total number of reads that covering the position
	my $seq = $data[4]; #string that contains the sequences of reads covering this position 
	my $quality = $data[5]; #string that contains the bases quality scores at the position
	my $prev_char = ""; #stores the previous $char value
	my @keys = ("a","c","g","t","n","totalq","start","end","del","in","skip","low"); #list of keys to initialize in the hash except DELPH (deletion place holder) which does not have strand-specific representation
	foreach my $key (@keys) { #initializes counters for each category - both for the negative strand (lower case) and the positive strand (upper case) 
		$hash{$key}=0;
		$hash{uc($key)}=0;
	}
	$hash{"DELPH"}=0; #initialize DELPH key's value - DELPH stands for deletion place holder which does not have strand-specific representation 
	while ($seq ne "") {
		$seq =~ s/^(.)//; #trims the first character from the $seq string and assign it to $1 special variable
		my $char = $1; #save the character for further analysis
		# print $seq,"\n";
		if ($char eq ".") {
			$char = uc($ref); #if reference nucleotide is known replace . (i.e match to the positive strand of the reference) with the nucleotide itself
		} elsif ($char eq ",") {
			$char = lc($ref); #if reference nucleotide is known replace , (i.e match to the negative strand of the reference) with the nucleotide itself
		} 
		if ($char eq "^") {
			$seq =~ s/^.//;  #trims the read quality score character that comes after "^"
			my ($next_char) = $seq =~ /(^[ACGTN\.acgtn\,])/; #extract the next character which is assumed to be a nucleotide
			if ($next_char =~ /([ACGTN\.])/) {
				$char = "START"; #start of read that was mapped to the positive strand of the reference
			} else {
				$char = "start"; #start of read that was mapped to the negative strand of the reference
			}
		} elsif ($char eq "\$") {
			# print ">>",$prev_char,"\n";
			if ($prev_char =~ /(^[ACGTN]$)/) { #the character preceding "$" is assumed to be a nucleotide
				$char = "END"; #end of read that was mapped to the positive strand of the reference
			} else {
				$char = "end"; #end of read that was mapped to the negative strand of the reference
			}
		} elsif ($char eq "-") {
			$seq =~ s/^([0-9]+)//;  #trims the deletion length "-" and saves length in $1
			$seq =~ s/^([ACGTNacgtn]){$1}//;  #trims the and sequence after deletion length and saves the last nucleotide in $1
			# $seq =~ s/^([0-9]+)([ACGTNacgtn])+{$1}//;  #trims the deletion length and sequence after "-" and saves the first nucleotide
			if ($1 =~ /(^[ACGTN]$)/) { #the character following deletion length is assumed to be a nucleotide
				$char = "DEL"; # the deletion was mapped to the positive strand of the reference
			} else {
				$char = "del"; #the deletion was mapped to the negative strand of the reference
			}
		} elsif ($char eq "+") {
			$seq =~ s/^([0-9]+)//;  #trims the insertion length "-" and saves length in $1
			$seq =~ s/^([ACGTNacgtn]){$1}//;  #trims the and sequence after insertion length and saves the last nucleotide in $1
			# $seq =~ s/^([0-9]+)([ACGTNacgtn])+{$1}//;  #trims the insertion length and sequence after "-" and saves the first nucleotide
			if ($1 =~ /(^[ACGTN]$)/) { #the character following insertion is assumed to be a nucleotide
				$char = "IN"; # the insertion was mapped to the positive strand of the reference
			} else {
				$char = "in"; #the insertion was mapped to the negative strand of the reference
			}
		} elsif ($char =~ /(^[ACGTNacgtn<>*]$)/) { #character is a nucleotide,reference skip or deletion (in case of reference skip or deletion place holder mpileup keeps the quality of the next base after the N or D cigar, i.e in the next M cigar) 
			$quality =~ s/^(.)//; #trims the first base quality score and assign it to $1 special variable
			my $score = ord($1) - $q_offset; #calculates and stores the base quality phred score
			if ($char eq ">") { 
				$char = "SKIP"; #reference skipping (N cigar) of read that was mapped to the positive strand of the reference 
			} elsif ($char eq "<") {
				$char = "skip"; #reference skipping (N cigar) of read that was mapped to the minus strand of the reference 
			} elsif ($char eq "*") {
				$char = "DELPH"; #deletion place holder 
			} elsif ($score < $min_q) { #mark low quality bases that do not meet the quality criterion that was specified 
				if ($char =~ /(^[ACGTN]$)/) { #the low quality base is located on a read that was mapped to the positive strand of the reference 
					$char = "LOW";
				} else {
					$char = "low"; #the low quality base is located on a read that was mapped to the negative strand of the reference 
				}
			} else {
				; #do nothing. $char is a nucleotide with sufficient quality and its sequence and letter case (upper\lower) are used for counting
			}
		}
		$hash{$char}++; #increase the appropriate counter in %hash
		$prev_char = $char; #keeps the $char for the next iteration
	}
	foreach my $key ("a","c","g","t") { #calculates the totals for reads with quality >= $min_q
		$hash{"totalq"} += $hash{$key};
		$hash{uc("totalq")} += $hash{uc($key)};
	}
	# use Data::Dumper; 
	# print Dumper %hash;
	print $fh2 join("\t",($chr,$start,$end,$ref,$total)); #print position coordinates to the output 
	foreach my $key (@keys) { #print the content of the counters for the positive strand
		print $fh2 "\t",$hash{uc($key)};
	}
	foreach my $key (@keys) { #print the content of the counters for both strands
		print $fh2 "\t",$hash{uc($key)}+$hash{$key};
	}
	print $fh2 "\t",$hash{"DELPH"}; #print DELPH key's value - DELPH stands for deletion place holder which does not have strand-specific representation
	print $fh2 "\n"; #move to next line in output
}

close $fh1;
close $fh2;
