#!/usr/bin/perl 
use strict;
use warnings; 

#This script gets sam file and output bed (bed6) file with all mismatches in the read level. Several filters are available: base quality, mapping quality, avoiding reads ends, avoiding bps around splice sites and ignoring duplicated reads.
#The script will output the mismatches relative to the positive DNA strand unless strand-specific data is given. In this case the read which include the RNA sequence should be specified.
#The bed output is as follows (tab delimited): chromosome, start coordinate (0-based),end coordinate(1-based), bed name, bed score: base quality, strand: mismatch strand.
#The name field (4th field in the output) is designed (semi-column ";" delimited): read name, read orientation (F or R), read "end" (1 or 2), mismatch type, mismatch position along the read sequence, mismatch position along the alignment,alignment length, read flag.
#INDELs are not reported.
#MD tag must be compatible with CIGAR string (e.g. number of matched/deleted/inserted bps etc). Discrepancy might occur after applying tools that change the CIGAR on the input BAM/SAM file.   

sub uniq { #adopted from http://stackoverflow.com/questions/7651/how-do-i-remove-duplicate-items-from-an-array-in-perl
    my %seen;
    grep !$seen{$_}++, @_;
}

sub compSeq {
#Function gets a DNA sequence and returns its compliment in the same directionality. Function is compatible
#with the IUPAC code and therefore recognize the following letters: ACTGRYSWKMBDHVN.
#Non standard characters will not be handled and will remain intact in the function output. 
	my ($seq) = @_;
	chomp($seq);
	my $comp_seq = "";
	for (my $i = 0; $i < length($seq); $i++) { #Loop is scanning the sequence starting from its end
		my $nuc = substr($seq,$i,1); #Reads letters one by one and translate them
		$nuc =~ tr/ACTGRYSWKMBDHVNactgryswkmbdhvn/TGACYRSWMKVHDBNtgacyrswmkvhdbn/; 
		$comp_seq .= $nuc; #Complimentary letters are chained to a reverse complement sequence
	}
	return $comp_seq;
}

open (my $fh1, "<", $ARGV[0]) or die "Can't open data file: $ARGV[0]\n"; #(STAR) sam file with MD tags.
open (my $fh2, ">", $ARGV[1]) or die "Can't open output file: $ARGV[1]\n"; #bam6 file.
open (my $fh3, ">", $ARGV[2]) or die "Can't open log file: $ARGV[2]\n"; #Log file.
open (my $fh4, ">", $ARGV[3]) or die "Can't open remove reads list file: $ARGV[3]\n"; #list of name of removed reads.
open (my $fh5, ">", $ARGV[4]) or die "Can't open remove output SAM file: $ARGV[4]\n"; #output SAM file with all the reads that were not excluded. This is NOT a standard SAM file since it may contain paired reads for whom one mate was included in the file and the other mate was excluded due to filters. 
my $md_col=$ARGV[5]//=12; #column number in-which the MD flag is found (1-based). if set to 0 the script will look for the correct column automatically.
my $exclude_tags=$ARGV[6]//=""; #BAM tags to exclude separated by commas (",")
my $min_qual=$ARGV[7]//=30; #min base quality for calling a mismatch.
my $offset=$ARGV[8]//=33; #quality score offset (def: ASCII-33).
my $avoid_bp_5p=$ARGV[9]//=5; #avoid calling mismatches if occur within N bp from 5' end.  
my $avoid_bp_3p=$ARGV[10]//=5; #avoid calling mismatches if occur within M bp from 3' end.
my $avoid_bp_ss=$ARGV[11]//=4; #avoid calling mismatches if occur within L bp around splice site.
my $avoid_bp_5p_aln=$ARGV[12]//=5; #avoid calling mismatches if occur within O bp from 5' end OF THE ALIGNMENT.  
my $avoid_bp_3p_aln=$ARGV[13]//=5; #avoid calling mismatches if occur within P bp from 3' end OF THE ALIGNMENT.
my $avoid_homopol=$ARGV[14]//=5; #specify the number of base-pairs to inspect from each side of the mismatch in order to avoid homopolymers (e.g. AAAAA). use 0 to not exclude homopolymers. 
my $mm_types_str=$ARGV[15]//="ALL"; #specify which types of mismatches to output. if set to "ALL" no type check will be preformed. mismatch types should be given as pair of capitalized letters, separated by commas with no spaces, i.e. Refnuc1Altnuc1,Refnuc2Altnuc2... (e.g. AG,TC,CT). Mind that if input data is not strand-specific both forward and reverse mismatch types should be noted (e.g. for recalling AG mismatches the value must be set to "AT,TC")
my $max_mms=$ARGV[16]//=0; #ignore reads that contains number of mismatches higher that this value. use 0 to avoid filtering for maximal number of mismatches.
my $max_mm_types=$ARGV[17]//=0; #ignore reads that contains number of TYPES of mismatches higher that this value. use 0 to avoid filtering for maximal number of mismatches' type. use 1 to allow only 1 mismatch TYPE on the same read (several mutations from the same type will pass).
my $min_mms_for_max_mm_types=$ARGV[18]//=2; #specify the minimal number of mismatches the read should contain in order to preform type-based filtering (as specified by $max_mm_types). using 3, for instance, will allow reads with 2 different mismatch TYPES to pass.
my $when_apply_max_mms=$ARGV[19]//=1; #0- apply max_mms/max_mm_types filters BEFORE mismatch-level filters (such as ends-avoidance, base quality etc.). 1- (or any nonzero value) apply  max_mms/max_mm_types filters AFTER filtering out undesired mismatches  
my $low_qual_reads_qual=$ARGV[20]//=20;  #ignore reads that do not have at least certain percent of bases (as defined by $low_qual_reads_pcent) in this base quality. use 0 to avoid this filter
my $low_qual_reads_pcent=$ARGV[21]//=0.75; #ignore reads that do not have at least this percent of bases in  a certain base quality(as defined by $low_qual_reads_qual). use $low_qual_reads_qual=0 to avoid this filter.
# my $ignore_dup=$ARGV[9]//=1; #skip optical/PCR duplicate in analysis. -1=duplicated reads would be included. 1 (or any other value)= skip duplicated reads
my $ignore_reads=$ARGV[22]//=1804; #skip reads with any of the specified flags set optical. Default 1804 (HEX 0x70C) excludes record that include: read unmapped, mate unmapped, not primary alignment, read fails platform/vendor quality checks, read is PCR or optical duplicate. 0- do not exclude any record 
my $min_map_qual=$ARGV[23]//=0; #mapping quality (0-255). STAR aligner will output 255 for any uniquely-mapped read. 
my $strand_sp=$ARGV[24]//=3; #For strand-specific sequencing. Determines which read (1 or 2) is mapped to the coding strand. Default: 3- Data is not strand specific: mismatches will be outputted relative to the DNA positive strand.  

my @mm_types; #an array to store the mismatches types to check for
if ($mm_types_str ne "ALL") {
	@mm_types=split(",",$mm_types_str);
	foreach my $cur_type (@mm_types) {
		die "Wrong mismatch types input" if ($cur_type !~ m/[ACGT][ACGT]/);
	}
}
my $avoid_homopol_minus_one = $avoid_homopol-1; #for the regular expression
my @exclude_tags_list=split(",",$exclude_tags); #split TAGS names to exclude
LINE: while (my $line=<$fh1>) { #reading sam file's lines
	chomp $line;
	if ($line =~ m/^@/) { #skips sam header lines, if exist. 
		print $fh5 $line,"\n"; #write the sam header line to the sam output
		next LINE;
	}
	my @data=split("\t",$line); #splitting fields in each line to an array 
	my ($name,$flag,$chr,$coord,$map_qual,$cigar,$seq,$qual)=($data[0],$data[1],$data[2],$data[3],$data[4],$data[5],$data[9],$data[10]); #retrieving read data: read name, flag, aligned chromosome, aligned start coordinate (1-based), CIGAR string, read sequence, read base quality.
	if ($flag & $ignore_reads){ #filter undesired reads as specified in the flags set in $ignore_reads. 
		print $fh4 $name,"\t",$flag,"\n"; #write excluded read name and flag to a file
		# print $fh5 $line,"\n"; #write the sam record to the sam output
		next LINE;
	}
	my $md_tag; #declaring an variable to store the MD tag
	if ($md_col > 0) { #retrieve MD tag: 
		$md_tag = $data[$md_col-1]; #FROM a pre-defined column OR:
	} else { #FROM a column which found to contain MD tag
		my $i = 11; #setted to the index of the first cell in the array that contains a TAG
		MD_FIND: until (not(defined($data[$i]))) { #iterating on the TAGs columns to look for the MD tag
			$md_tag = $data[$i];
			last MD_FIND if ($md_tag =~ m/MD:Z:[0-9ACGTN\^]+/);
			undef($md_tag);
			$i++;
		}
		if (not(defined $md_tag)) { #no MD tag was found
			print $fh3 "Read ${name} do not contains a supported MD TAG. Moving to next read.\n";
			print $fh5 $line,"\n"; #write the sam record to the sam output
			next LINE; 
		}
	}
	EXCLUDE_TAG: foreach my $tag_name (@exclude_tags_list) {
		my $i = 11; #setted to the index of the first cell in the array that contains a TAG
		TAG_FIND: until (not(defined($data[$i]))) { #iterating on the TAGs columns to look for the MD tag
			my $cur_tag = $data[$i];
			if ($cur_tag =~ m/\Q$tag_name\E/) {
				print $fh4 $name,"\t",$flag,"\n"; #write excluded read name and flag to a file
				next LINE;
			}
			$i++;
		}
	}
	#FILTERING (read-level):
	if ($map_qual < $min_map_qual) { #filter out flags that do not meet the mapping quality criterion  
		print $fh4 $name,"\t",$flag,"\n"; #write excluded read name and flag to a file
		next LINE;
	}
	# next if (($ignore_dup != -1) and ($flag & 1024)); #skips duplicated read if asked for. 
	# next LINE if ($flag & $ignore_reads); #filter undesired reads as specified in the flags set in $ignore_reads. 
	# print "~",$cigar,"\n";
	my $read_len = length($seq); #calculating read length
	my $qual_for_low_qual_reads = $qual; #gets the base quality strings to filter low quality reads
	my $high_qual_bases = 0; #counter for the number of bases with at least the requested quality 
	LOW_QUAL: while ($qual_for_low_qual_reads ne "") {
		my $base_qual = ord(substr($qual_for_low_qual_reads,0,1,""))-$offset; #retrieve the base quality from the quality string
		# print $base_qual,"\n";
		$high_qual_bases++ if ($base_qual >= $low_qual_reads_qual); #increase the counter value if the current base quality is equals or higher to the requested quality
	}
	if (($high_qual_bases/$read_len) < $low_qual_reads_pcent) { #exclude the read if it does not have bases in the sufficient quality in the requested fraction of read length 
		print $fh4 $name,"\t",$flag,"\n"; #write excluded read name and flag to a file
		next LINE;
	}
	# my $md1 = substr($md_tag,rindex($md_tag,":")+1); #extract the MD string.
	my ($md) = $md_tag =~ /MD:Z:([0-9ACGTN\^]+)/; #extract the MD string.
	# print ">",$md,"\n";
	my (@m_cigar) = $cigar =~ /(\d+)M/g;  #store the number of matched bases in array
	# print join("\t",@m_cigar) if ($cigar eq "76S"); #DEBUG
	unless (@m_cigar) { #move to next line if no matching base-pairs are found. This situation can occur if the SAM file was treated with external tools after the alignment (e.g. removing overlapping sequences in mates pair).
		print $fh4 $name,"\t",$flag,"\n"; #write excluded read name and flag to a file
		# print "~",$cigar,"\n";
		# print $fh5 $line,"\n"; #write the sam record to the sam output
		next LINE;
	}
	my $aln_len = 0; 
	foreach my $cur_m_cigar (@m_cigar) {
		$aln_len += $cur_m_cigar; #calculates the total length of the alignment 
	}
	my ($read_ornt,$avoid_bp_left,$avoid_bp_right,$avoid_bp_left_aln,$avoid_bp_right_aln) = ("F",$avoid_bp_5p,$read_len-$avoid_bp_3p,$avoid_bp_5p_aln,$aln_len-$avoid_bp_3p_aln); #sets the read orientation to forward "F" as default , rightmost nucleotide in 5' end (which is in the left) of read sequence to avoid, leftmost nucleotide in 3' end (which is in the right) of read sequence to avoid, rightmost nucleotide in 5' end (which is in the left) of alignment to avoid, leftmost nucleotide in 3' end (which is in the right) of alignment to avoid.
	if ($flag & 16) { #sets the read orientation to reverse "R", rightest nucleotide in 5' end of sequnece  (which is in the left) to avoid, leftmost nucleotide in 3' end (which is in the right) ONLY if read is in reverse orientation. Same is true for the ends of the alignment.
		($read_ornt,$avoid_bp_left,$avoid_bp_right,$avoid_bp_left_aln,$avoid_bp_right_aln) = ("R",$avoid_bp_3p,$read_len-$avoid_bp_5p,$avoid_bp_3p_aln,$aln_len-$avoid_bp_5p_aln);
	}
	my $read_end = 1; # sets the read "end" to be 1 as default.
	if ($flag & 128) { #sets the read "end" to be 2 if the read in the second in pair.
		$read_end = 2;
	}
	my ($pos_seq,$pos_aln,$prev_cigar_n) = (0,0,0); #position while iterating on $seq i.e. position on read (0-based), position within the alignment (i.e. the aligned bps)  indicates whether the previous CIGAR was N (i.e. splicing gap).
	undef(my $md_element); #For parts of the MD string. Set as undefined
	my @cur_read_mms=(); #stores the mismatches types found on this read
	my $cur_read_output=""; #stores the output for the current read - namely, mismatches that passed the mismatch-level filtration.
	my $cigar_for_legnth=$cigar;
	my $aln_length=0; # number of nucleotide involve in the alignment. Matched bases (M) and Inserted bases (I) are counted to assess this value.  
	CIGAR_LENGTH: while ($cigar_for_legnth ne "") { #The CIGAR temporary string is being shortened in each iteration to assess the length of alignment
		my ($cigar_num,$cigar_letter) = ($cigar_for_legnth =~ /(\d+)([MNSHID])/); #CIGAR LETTERS ARE DIFFINED IN THIS LINE AND IN THREE OTHER LINES.
		# print $cigar_num;
		$cigar_for_legnth =~ s/\d+[MNSHID]//; #trimming the first alphabetical of CIGAR. CIGAR LETTERS ARE DIFFINED IN THIS LINE AND IN THREE OTHER LINES.
		if ($cigar_letter =~ m/[MI]/) {
			$aln_length+=$cigar_num;
		}
	}
	CIGAR: while ($cigar ne "") { #The CIGAR string is being shortened in each iteration and the values are being assigned to $cigar_num and $cigar_letter.
		# print "!",$data[5],"\n";
		# print "@",$cigar,"\n";
		my ($cigar_num) = $cigar =~ /(\d+)/; #retrieving the first numeral of CIGAR.
		# print "#",$cigar,"\n";
		$cigar =~ s/\d+//; #trimming the first numeral of CIGAR.
		# print "\$",$cigar,"\n";
		my ($cigar_letter) = $cigar =~ /([MNSHID])/; #retrieving the first alphabetical of CIGAR. CIGAR LETTERS ARE DIFFINED IN THIS LINE AND IN THREE OTHER LINES.
		$cigar =~ s/[MNSHID]//; #trimming the first alphabetical of CIGAR. CIGAR LETTERS ARE DIFFINED IN THIS LINE AND IN THREE OTHER LINES.
		# print "%",$cigar,"\t",$cigar_num,"\t",$cigar_letter,"\n";
		my $next_cigar_n = $cigar =~ /^\d+N/; #checking if the next cigar element (number+letter) is intron gap ("N"). equals 1 if true or undef if false
		$next_cigar_n=0 unless (defined($next_cigar_n)); # if $next_cigar is undef assign 0  
		if ($cigar_letter eq "H") { #read hard-clipping
			; #do nothing. PAY ATTENTION that mismatch position along the read would not be correct in case of hard clipping and would refer to the location in the clipped sequence and not to the original read sequence. 
		} elsif ($cigar_letter eq "S") { #read soft-clipping
			$pos_seq+=$cigar_num; #moving forward on the sequence 
		} elsif ($cigar_letter eq "M") { #match
			# print "^",$cigar_letter,"\n";
			my $original_cigar_m  = $cigar_num; #storing the length of current M CIGAR element 
			# print "&",$md,"\t",$md_element,"\n";
			MD: while (($md ne "") or (defined($md_element))) { #iterating on MD string
				# print "*",$md,"\t",$md_element,"\n";
				if (not defined($md_element)) {
					($md_element) = $md =~ /^(\d+|[ACTGN])/;  #retrieving the first leftest element of the MD string if it's a number
					# print "(",$md,"\t",$md_element,"\n";
					$md =~ s/^\d+|[ACTGN]//; #trimming the first leftest MD element.
					# print ")",$md,"\t",$md_element,"\n";
				}
				if ($md_element =~ m/^\d+/) { #the element is a number 
					# print "-",$md,"\t",$md_element,"\n";
					if ($cigar_num > $md_element) {
						$pos_seq+=$md_element; #moving forward on the sequence.
						$pos_aln+=$md_element; #moving forward on the alignment.
						$coord+=$md_element; #increasing coordinate in agreement with perfectly-matched bases.
						$cigar_num-=$md_element; #moving "forward" on CIGAR string.
						undef($md_element); #unsetting $md_element in order to retrieving the next element in the next iteration MD LOOP. 
					} else {
						$pos_seq+=$cigar_num; #moving forward on the sequence.
						$pos_aln+=$cigar_num; #moving forward on the alignment.
						$coord+=$cigar_num; #increasing coordinate in agreement with perfectly-matched bases.
						$md_element-=$cigar_num; #moving "forward" on MD string.
						# $cigar_num=0; 
						next CIGAR; #moving to the next iteration of CIGAR LOOP in order to retrieving the next CIGAR element.
					}
					# print "_",$cigar_num,"\n";
				} else {
					# print "+",$md,"\t",$md_element,"\n";
					if ($md_element !~ m/^N/) { #the element is a known nucleotide /[ACGT]/
						my $base_qual = ord(substr($qual,$pos_seq,1))-$offset; #retrieve the base quality from the quality string
						my $mm_type=$md_element.substr($seq,$pos_seq,1); #mismatch type is two capital letters without any character in-between: (reference nucleotide)(read nucleotide)
						push (@cur_read_mms,$mm_type) if ($when_apply_max_mms == 0); #save list of ALL mismatches for read-level mismatch-dependent filters which are been applied after the inspection of the whole read
						my ($up_seq,$down_seq)=("","");
						if ($pos_seq - $avoid_homopol >= 0) { #>= since $pos_seq is 0-based
							$up_seq = substr($seq,$pos_seq - $avoid_homopol,$avoid_homopol); #retrieve base-pair upstream to the mismatch for homopolymer check 
						} 
						if ($read_len - $pos_seq > $avoid_homopol) { 
							$down_seq = substr($seq,$pos_seq + 1 ,$avoid_homopol); #retrieve base-pair downstream to the mismatch for homopolymer check 
						}
						# print $up_seq,"_",$down_seq,"\n";
						# print "Voila\n" if ($down_seq =~ m/^([ACGT])(\1){\Q$avoid_homopol_minus_one\E}/);
						#FILTERING (mismatch-level):
						if (($pos_seq >= $avoid_bp_left) and #>= since $pos_seq is 0-based
							($pos_seq < $avoid_bp_right) and 
							($cigar_num > $avoid_bp_ss*$next_cigar_n) and
							($pos_aln >= $avoid_bp_left_aln) and #>= since $pos_seq is 0-based
							($pos_aln < $avoid_bp_right_aln) and 
							($original_cigar_m - $cigar_num >= $avoid_bp_ss*$prev_cigar_n) and #if previous CIGAR element is not N the expression on the left equals 0.
							($base_qual >= $min_qual) and
							((not($up_seq =~ m/^([ACGT])(\1){\Q$avoid_homopol_minus_one\E}/)) and (not($down_seq =~ m/^([ACGT])(\1){\Q$avoid_homopol_minus_one\E}/)))) {
							# my $mm_type=$md_element.substr($seq,$pos_seq,1); #mismatch type is two capital letters without any character in-between: (reference nucleotide)(read nucleotide)
							push (@cur_read_mms,$mm_type) if ($when_apply_max_mms != 0); #save list of only FILTERED-IN mismatches for read-level mismatch-dependent filters which are been applied after the inspection of the whole read
							my $type_check=1; #indicates by default that the current mismatch type should be outputted, i.e. assuming that $mm_types_str is set to "ALL"
							if (@mm_types) { #if TRUE the array is defined, meaning that $mm_types_str is set to value other than "ALL"
								$type_check=0;
								MM_TYPES: foreach my $cur_type (@mm_types) { #checking if the current mismatch is of one of the types that should be outputted 
									if ($cur_type eq $mm_type) {
									$type_check=1;
									last MM_TYPES; #if the type is in the list stops the search
									}
								}
							}
							#FILTERING (mismatch-type):
							if ($type_check == 1) { #current mismatch is one of the types that should be outputted
								# print $mm_type,"\n";
								my $mm_pos = $pos_seq + 1; # sets the mismatch position along the read to $pos_seq which is the correct position if the read was mapped to the positive DNA strand. +1 since $pos_seq is 0-based but $mm_pos is 1-based 
								my $mm_pos_aln = $pos_aln + 1; # sets the mismatch position along the alignment to $pos_aln which is the correct position if the read was mapped to the positive DNA strand. +1 since $pos_aln is 0-based but $mm_pos_aln is 1-based 
								if ($read_ornt eq "R") {
									$mm_pos = $read_len - $pos_seq; # fixes the mismatch position along the read in case the read was mapped to the negative DNA strand. (Read length is "1-based" so adding +1 to $pos_seq is not necessary this time 
									$mm_pos_aln = $aln_len - $pos_aln; # fixes the mismatch position along the alignment in case the read was mapped to the negative DNA strand. (Alignment length is "1-based" so adding +1 to $pos_aln is not necessary this time 
								}
								my $mm_strand = "+"; # sets the strand of the mismatch to positive DNA strand as default. The strand would be remained unchanged in case of non-strand-specific RNA-seq.
								if (($strand_sp==2) or ($strand_sp==1)) {
									if ($strand_sp == $read_end) {
										($mm_strand = $read_ornt) =~ tr/FR/\+\-/; #if the read contains the RNA sequence the strand to which the read was mapped is the strand of the mismatch. F or R annotation in the read level is translated to + or - annotation, respectively, relative to reference sequence.
									} else {
										($mm_strand = $read_ornt) =~ tr/FR/\-\+/; #if the read contains the revere compliment sequence relative to the RNA the strand to which the read was mapped is the opposite strand of the mismatch. F or R annotation in the read level is translated to - or + annotation, respectively, relative to reference sequence.
									}
									$mm_type=&compSeq($mm_type) if $mm_strand=="-"; #change the mismatch type to the complement sequence if it's found on the the negative strand.
								}
								my $bed_name = join(";",$name,$read_ornt,$read_end,$mm_type,$mm_pos,$mm_pos_aln,$aln_length,$flag); # sets the bed name file, which be printed in the 4th field in the output: read name, read "end" (1 or 2), mismatch type, mismatch position along the read, mismatch position along the alignment,alignment length, read flag.
								# print $fh2 join("\t",($chr,$coord-1,$coord,$bed_name,$base_qual,$mm_strand)),"\n"; #prints bed file record (tab delimited): chromosome, start coordinate (0-based),end coordinate(1-based), bed name, bed score: base quality, strand: mismatch strand.
								$cur_read_output .= join("\t",($chr,$coord-1,$coord,$bed_name,$base_qual,$mm_strand))."\n"; #prints bed file record (tab delimited): chromosome, start coordinate (0-based),end coordinate(1-based), bed name, bed score: base quality, strand: mismatch strand.
								# print $cur_read_output;
							}
						}
					}
					$coord++; #increasing coordinate by 1 position (Either if the md element is known nucleotide /[ACGT]/ or unknown nuclotide /N/).
					$pos_seq++; #moving forward on the sequence by 1 position (Either if the md element is known nucleotide /[ACGT]/ or unknown nuclotide /N/).
					$pos_aln++; #moving forward on the sequence by 1 position (Either if the md element is known nucleotide /[ACGT]/ or unknown nuclotide /N/).
					$cigar_num--; #moving "forward" on CIGAR string by 1 position (Either if the md element is known nucleotide /[ACGT]/ or unknown nucleotide /N/).
					undef($md_element); #unsetting $md_element in order to retrieving the next element in the next iteration MD LOOP.
					# print "=",$cigar_num,"\n";
				}
			}
		$prev_cigar_n = 0; #indicates for the next CIGAR LOOP iteration that the previous CIGAR element is not a splicing gap ("N").
		} elsif ($cigar_letter eq "N") { #splicing gap (a.k.a intron)
			$prev_cigar_n = 1; #indicates for the next CIGAR LOOP iteration that the previous CIGAR element is a splicing gap ("N").
			$coord+=$cigar_num; #increasing coordinate in agreement with splicing gap.
		} elsif ($cigar_letter eq "I") { #insetrion (in the read)
			$pos_seq+=$cigar_num; #moving forward on the sequence by the insertion size. $md_element or $md is unchanged since insertions to the reads are not represented in it but only in the CIGAR string. $coord is unchanged since the reference do not include the inserted sequence
		} elsif ($cigar_letter eq "D") { #deletion (in the read)
			$coord+=$cigar_num; #increasing coordinate in agreement with deletion gap. $pos_seq is unchanged since the deleted sequence is not found in the read.
			# print "1",$md,"\n";
			$md =~ s/^\^[ACGTN]{\Q$cigar_num\E}//; #trimming the annotation for deleted sequence from MD string. 
			# print "2",$md,"\n";
			# $md =~ s/^\^[ACGTN]//; #trimming the annotation for deleted sequence from MD string. 
		} else {
			print $fh3 "Read ${name}.${read_end} contains an unsupported CIGAR. Moving to next read.\n";
			print $fh5 $line,"\n"; #write the sam record to the sam output
			next LINE; 
		}
	}
	my $cur_read_mms_num = scalar(@cur_read_mms); #counts the number of mismatches in the read
	my $cur_read_mms_types_num = scalar(&uniq(@cur_read_mms)); #counts the number of mismatches' types in the read
	#FILTERING (read-level mismatch-dependent):
	# print $cur_read_mms_num," ",$cur_read_mms_types_num,"\n";
	unless ((($cur_read_mms_num > $max_mms) and ($max_mms != 0)) or (($cur_read_mms_types_num > $max_mm_types) and ($max_mm_types != 0) and ($cur_read_mms_num >= $min_mms_for_max_mm_types))) {
		print $fh2 $cur_read_output; #print the records for the current read if the filters were passed.
		print $fh5 $line,"\n"; #write the sam record to the sam output
	} else {
		print $fh4 $name,"\t",$flag,"\n"; #write excluded read name and flag to a file
	}
}

close $fh1;
close $fh2;
close $fh3;
close $fh4;
close $fh5;
