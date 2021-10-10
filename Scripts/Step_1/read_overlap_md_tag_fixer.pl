#!/usr/bin/perl 
use strict;
use warnings; 


#THIS SCRITP GETS SAM FILE IN WHICH OVERLAPPING BASES IN ONE MATE WERE SOFT-CLIPPED AND ORIGINAL CIGAR WAS KEPT IN A SPECIFIC BAM TAG.
#THE SCRIPT LOOKS FOR THE MD TAG AND FIXS IT TO MATCH TO THE NEW SOFT-CLIPPED CIGAR.
#OUTPUTS SAM FILE TO STDOUT.



open (my $fh1, "<", $ARGV[0]) or die "Can't open data file: $ARGV[0]\n"; #SAM file
my $cigar_tag_name=$ARGV[1]//="CG"; #name of tag in-which the original CIGAR is stored. The tag is assumed to be found only in reads with overlapping bps.  
my $md_tag_name=$ARGV[2]//=""; #name of tag in-which the original MD tag would be is stored. If equals "" the original MD would not be saved.

while (my $line=<$fh1>) {
	if ($line =~ m/^@/) {
		print $line;
	} else {
		chomp $line;
		my @data=split("\t",$line);
		# next if ($data[5] =~ /^[0-9]+S$/); #skips reads that were entirely soft-clipped 
		if ($data[5] =~ /^[0-9]+S$/) { #skips reads that were entirely soft-clipped 
			print $line,"\n";
			next;
		}
		my $i = 11; #setted to the index of the first cell in the array that contains a TAG
		my $old_cigar=""; 
		CIGAR_FIND: until (not(defined($data[$i]))) { #iterating on the TAGs columns to look for the CIGAR tag
			last CIGAR_FIND if (($old_cigar) = $data[$i] =~ /\Q$cigar_tag_name\E:Z:(([0-9]+[MIDNSHPX=])+)/);
			$i++;
		}
		if ($i > $#data) { #if no CIGAR tag was found
			print $line,"\n";
		} else {
			my $md="";
			$i = 11; #setted to the index of the first cell in the array that contains a TAG
			MD_FIND: until (not(defined($data[$i]))) { #iterating on the TAGs columns to look for the MD tag
				last MD_FIND if (($md) = $data[$i] =~ /MD:Z:([0-9ACGTN\^]+)/);
				$i++;
			}
			my $old_md = $md; #store the original MD tag
			my (@i_number_new_arr) = $data[5] =~ m/([0-9]+)I/g; #retrieving the number of bases in insertion
			my (@i_number_old_arr) = $old_cigar =~ m/([0-9]+)I/g; #retrieving the number of bases in insertion in the original alignment
			my ($i_number_new,$i_number_old,$i_number)=(0,0,0);
			if (@i_number_new_arr) { 
				# print join("\t",@i_number_new_arr,$#i_number_new_arr),"\n";
				foreach my $val (@i_number_new_arr) {
					$i_number_new+=$val;
				}
			}
			if (@i_number_old_arr) { #summing the number of of bases in insertion in the original alignment
				# print join("\t",@i_number_old_arr,$#i_number_old_arr),"\n";
				foreach my $val (@i_number_old_arr) {
					$i_number_old+=$val;
				}
			}
			$i_number = $i_number_old-$i_number_new; #calculating the number of bases in insertion that should be added to MD
			my ($s_number_new) = $data[5] =~ /^([0-9]+)S/; #retrieving the number of left soft-clipped bases
			my ($s_number_old) = $old_cigar =~ /^([0-9]+)S/; #retrieving the number of left soft-clipped bases in the original alignment 
			# if ((defined($s_number_new) and not defined($s_number_old)) or ($s_number_new > $s_number_old)) { #checks if the read was clipped from its left end
			if (defined($s_number_new) and (not defined($s_number_old) or ($s_number_new > $s_number_old))) { #checks if the read was clipped from its left end
				if (defined($s_number_old)) {
					$s_number_new -= $s_number_old; #remove the originally soft-clipped bases from the new soft-clipping value
				}
				$s_number_new -= $i_number; #for compensating for inserted reads that are not represented in MD tag 
				while ($s_number_new > 0) { #fixing MD tag to match the new CIGAR
					# print ">",$md,"\t",$i,"\n";
					my ($md_element) = $md =~ /^([0-9]+|[ACGTN]|\^[ACGTN]+)/; #extracting the leftmost element in MD tag
					# print $md_element,"\t",$md,"\n";
					$md =~ s/^([0-9]+|[ACGTN]|\^[ACGTN]+)//;  # removing the leftmost element in MD tag
					if ($md_element =~ m/^[0-9]+/) {
						if ($s_number_new >= $md_element) {
							$s_number_new -= $md_element;
						} else {
							$md_element -= $s_number_new; #calculates the "remain" of MD "match" element
							$s_number_new = 0; #ends the WHILE
							$md = $md_element.$md; #pasting the "remain" of MD "match" element in the leftmost end of MD tag 
						}
					} elsif ($md_element =~ m/^[ACGTN]/) {
						$s_number_new--; #a mismatch in the former alignment
					} #else: $md_element=~/\^[ACGTN]+/, i.e. deletion from the read. therefore no action should be taken but removing the deletion annotation from the MD tag (as already done above).
				}
			} else { #the read was clipped from its right end
				($s_number_new) = $data[5] =~ /([0-9]+)S$/; #retrieving the number of right soft-clipped bases
				($s_number_old) = $old_cigar =~ /([0-9]+)S$/; #retrieving the number of right soft-clipped bases in the original alignment 
				if (defined($s_number_old)) {
					$s_number_new -= $s_number_old; #remove the originally soft-clipped bases from the new soft-clipping value
				}
				$s_number_new -= $i_number; #for compensating for inserted reads that are not represented in MD tag 
				while ($s_number_new > 0) { #fixing MD tag to match the new CIGAR
					my ($md_element) = $md =~ /([0-9]+|[ACGTN]|\^[ACGTN]+)$/; #extracting the rightmost element in MD tag
					# print $md_element,"\n";
					$md =~ s/([0-9]+|[ACGTN]|\^[ACGTN]+)$//;  # removing the rightmost element in MD tag
					if ($md_element =~ m/[0-9]+$/) {
						if ($s_number_new >= $md_element) {
							$s_number_new -= $md_element;
						} else {
							$md_element -= $s_number_new; #calculates the "remain" of MD "match" element
							$s_number_new = 0; #ends the WHILE
							$md = $md.$md_element; #pasting the "remain" of MD match element in the rightmost end of MD tag
						}
					} elsif ($md_element =~ m/[ACGTN]$/) {
						$s_number_new--; #a mismatch in the former alignment
					} #else: $md_element=~/\^[ACGTN]+$/, i.e. deletion from the read. therefore no action should be taken but removing the deletion annotation from the MD tag (as already done above).
				}			
			}
			my $last_md_element;
			if ((($last_md_element) = $md =~ /([0-9]+)$/) and ($last_md_element==0)) {
				$md =~ s/0$//; #remove trailing zero if exists 
			} elsif ($last_md_element = $md =~ /\^[ACGTN]+$/) {
				$md =~ s/\^[ACGTN]+$//; #remove trailing deletion if exists 
			}
			if (my $first_md_element = $md =~ /^\^[ACGTN]+/) {
				$md =~ s/^\^[ACGTN]+//; #remove leading deletion if exists
			}
			if ($md_tag_name ne "") { #storing the original MD tag in a new tag
				push(@data,"${md_tag_name}:Z:${old_md}");
			}
			$data[$i] = "MD:Z:${md}"; #updating the MD tag
			print join("\t",@data),"\n";
		}
	}
}

close $fh1;
