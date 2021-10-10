#!/usr/bin/perl -w

#A script by Eli Eisenberg

use Math::CDF qw(:all);
my $prob = pnorm(1.96);

$err=0.001; #sequencing error
$snp_prob=0.0001; #heterozygous snp
$snp_prob2=$snp_prob*$snp_prob; #homozygous snp

$level_cut=0.001;
$score_cut=0.1;

open (my $fh1, "<", $ARGV[0]) or die "Can't open input file: $ARGV[0]\n"; #input file: chr,start,end,type,for each sample: supporting_reads, supporting_reads_plus, total_coverage,  total_coverage_plus
$Nsam=$ARGV[1]; #number of samples in input file 
# $Nsam=47;
#$Nsam=138;

# open LOG,">/tmp/log";

# $maxf=4*$Nsam+4; #number of fields in input
$maxf=$Nsam+4; #number of fields in input

while ($l=<$fh1>){
    chomp($l);
    @a=split("\t",$l);
    next if ($a[0] eq "chrM"); #filters undesired chrs
    next if ($a[0] eq "chrEBV");
    # next if($a[0] =~ /chr.*_alt/);
    # next if($a[0] =~ /chr.*_random/);
    # next if($a[0] =~ /chrUn.*/);
    $score=0;
    $scorep=0;
    $scorem=0;
    $scorep2=0;
    $scorem2=0;
    $n=0;
    $n1=0;
    $n2=0;
    $sup=0;
    $supmm=0;
    $supmm_plus=0;
    $sup_plus=0;

    for($i=4;$i<$maxf;$i++){
		next if ($a[$i] eq "NA"); #skip if no data 
		next if ($a[$i] eq "EX"); #skip if data was excluded (due to presence in WXS) 
		my @cov=split(";",$a[$i]); #split data string for each sample
		next if ($cov[0+2] == 0); #skip if no coverage
		$sup+=$cov[0+2]; #sum total covarage

		next if($cov[0] eq "NA"); #skip if no data
		$n++; #sum #tissues showing editing

#		$pmm = 1-pbinom($cov[0]-1,$cov[0+2],$err);
		$pmm = pbinom($cov[0+2]-$cov[0],$cov[0+2],1-$err);

	#	$psnp1=(1-pbinom($cov[0]-1,$cov[0+2],0.5))*$snp_prob;
		$psnp1=pbinom($cov[0+2]-$cov[0],$cov[0+2],0.5)*$snp_prob;

		$psnp2=(pbinom($cov[0],$cov[0+2],$err))*$snp_prob2;
		
		$p=$pmm; #assign the most probable p to $p
		if($p<$psnp1){$p=$psnp1;}
		if($p<$psnp2){$p=$psnp2;}
		$n1++ if($p==$psnp1);
		$n2++ if($p==$psnp2);

#		if($p<1e-8){print "LOW-p $p $cov[0] $cov[0+2]\n"; }

		$score+= (-1*log($p+1e-300)); 

		if($cov[3]>0){
		    $pmm_plus = pbinom($cov[3]-$cov[1],$cov[3],1-$err);

		    $scorep += -1*log($pmm_plus+1e-300);

		    $psnp1=pbinom($cov[3]-$cov[1],$cov[3],0.5)*$snp_prob;
		    $psnp2=(pbinom($cov[1],$cov[3],$err))*$snp_prob2;
		    if($pmm_plus<$psnp1){$pmm_plus=$psnp1;}
		    if($pmm_plus<$psnp2){$pmm_plus=$psnp2;}
		    
		    $scorep2 += -1*log($pmm_plus+1e-300);
		}
		if($cov[2]>$cov[3]){
		    $match_minus=$cov[2]-$cov[3]-$cov[0]+$cov[1];      
		    $pmm_minus = pbinom($match_minus,$cov[2]-$cov[3],1-$err);

		    $scorem += -1*log($pmm_minus+1e-300);

		    $psnp1=pbinom($match_minus,$cov[2]-$cov[3],0.5)*$snp_prob;
		    $psnp2=(pbinom($cov[0]-$cov[1],$cov[2]-$cov[3],$err))*$snp_prob2;
		    if($pmm_minus<$psnp1){$pmm_minus=$psnp1;}
		    if($pmm_minus<$psnp2){$pmm_minus=$psnp2;}

		    $scorem2 += -1*log($pmm_minus+1e-300);
		}


		$sup_plus+=$cov[0+3];
		$supmm+=$cov[0];
		$supmm_plus+=$cov[0+1];

		# print LOG "$a[0] $a[1] $a[2] $a[3] $i $a[$i] $a[$i+1] $a[$i+2] $pmm $psnp1 $psnp2 $p $score $n1 $n2 $sup $supmm $supmm_plus\n";
    }

    $ed_lev=$supmm/$sup;
    next if($ed_lev<$level_cut);
    $avscore=$score/$Nsam;
    next if($avscore<$score_cut);

    $match=$sup-$supmm;
    $match_plus=$sup_plus-$supmm_plus;

    $p_strand_err= 1;
    $p_strand_err= pbinom($match_plus,$sup_plus,1-$err) if ($sup_plus>0);
    $p2=1;
    $p2 = pbinom($match-$match_plus,$sup-$sup_plus,1-$err) if ($sup>$sup_plus);

    if($p2>$p_strand_err){$p_strand_err=$p2;}

    $min_mm=$supmm-$supmm_plus;
    $sup_min=$sup-$sup_plus; #total coverage on the - strand
    if($min_mm>$supmm_plus){$min_mm=$supmm_plus;$sup_min=$sup_plus}
    $p_strand=2*pbinom($min_mm,$supmm,0.5);
#    $p_strand_err= 1;
#    $p_strand_err= 1-pbinom($min_mm-1,$sup_min,$err) if $min_mm>0;

    $score2=$scorep;
    if($scorem<$score2){$score2=$scorem;}
    $score2=$score2/$Nsam;

    $score3=$scorep2;
    if($scorem2<$score3){$score3=$scorem2;}
    $score3=$score3/$Nsam;

    $score4=-log($p_strand_err+1e-300)/$Nsam;

	#The output format is: location, MM_type, score, #tissues showing editing, # where hetero SNP is the "best" explanation, same for homo SNP, tot cove, plus cov, tot supp, plus supp, old p_val, new p_va
	# print "$a[0] $a[1] $a[2] $a[3] $avscore \t\t\t $n $n1 $n2 $sup $sup_plus $supmm $supmm_plus $p_strand $p_strand_err\n";
	print join("\t",($a[0],$a[1],$a[2],$a[3],$avscore,$n,$n1,$n2,$sup,$sup_plus,$supmm,$supmm_plus,$p_strand,$p_strand_err,$score2,$score3,$score4)),"\n";
}

close $fh1;
