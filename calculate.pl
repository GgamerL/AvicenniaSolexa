#!perl -w
#written by Zhengzhen Wang
#use strict;
#The script is used in calculating Pi, D, S, Dxy and Fst for pairwise populations from ms simulation files.
#all done by 80 genes of 1000bp
#$temp_seqfile = "out";

y $ef_length =80000;
my $n1 = shift;
my $n2 = shift;
my $switch = 0;
my $site = 0;
#$n1 = 2 * $n1;
#$n2 = 2 * $n2;
#for my $n (1..$n1){
#	$pop{$n} = "pop1";
#}
#for my $n (($n1+1) ..($n1+$n2)){
#	$pop{$n} = "pop2";
#}
#read ms sim seq file into array 
my $seekpin = ($n1+$n2+4)*80;
my $l_count = 0;
my $l_seek = $seekpin;


readline(*STDIN);
readline(*STDIN);

while (<>){
	$l_count ++;
#	print "line $l_count\n";
#	if ($l_count > $l_seek){
#		print "Error:find $l_count\>$l_seek in term\n";
#	}
	if (/segsites\:\s(\d+)/){
#		print "seg at $l_count\n";
		$site = $site + $1;
		$site0 = $1;
		if ($site0 == 0){
#			print "Rat found line $l_count, seek $l_seek\t\=>";
			$l_seek = $l_seek - $n1 - $n2;
#			print "  $l_seek\t";
#			if ($l_count > $l_seek){
#				print "Error: find $l_count \> $l_seek\n";
#				
#			}
		}
		$switch = 0;
#		$s_count ++;
	}elsif( /^[0|1]+$/ ){
		s/\s+//;
#	print "zoo $_\n";
		my @line = split //;
		$switch  ++;
		for my $cal(0..($site0 - 1)){
			if ($switch < ($n1+1)){
				if (exists $pop1[$site - $site0 +$cal]){
					$pop1[$site - $site0 + $cal] = $pop1[$site - $site0 + $cal] + $line[$cal];
				}else{
					$pop1[$site - $site0 + $cal] = 0;
				}
			}elsif ($switch <($n1 +$n2)){
				if (exists $pop2[$site - $site0 + $cal]){
					$pop2[$site - $site0 + $cal] = $pop2[$site - $site0 + $cal] + $line[$cal];
				}else{
					$pop2[$site - $site0 + $cal] = 0;
				}
			}
		}
	
	}

###when it come to 80 gene's node###
if ($l_seek eq $l_count){
#	print "==found $l_seek at\t$l_count\t";
#	print "===";
	my $freq1; my $freq2; my $freq3;
	my $p1; my $p2; my $p3;
	my $seg1; my $seg2;
	my $in_pop1; my $in_pop2;
	my $Dxy;
	for my $n (1 ..$site){
	
		$p1 = $pop1[$n-1] / $n1;
		$p2 = $pop2[$n-1] / $n2;
		$p3 = ($p1 * $n1 + $p2 * $n2) / ($n1 + $n2);
		$freq1 += $n1 / ($n1 - 1) * $p1 * (1 - $p1);
		$freq2 += $n2 / ($n2 - 1) * $p2 * (1 - $p2);
		$freq3 += ($n1 + $n2) / ($n1 + $n2 -1) * $p3 *(1 - $p3);
#		$tetest = 1-($freq1+$freq2)/2/$freq3;
#		print "$freq1\t$freq2\t$freq3\t$tetest\n";
		$Dxy += $p1 * (1 - $p2) + $p2 * (1 - $p1);
	
		if ($p1 != 0 && $p1 != 1){
			$in_pop1 += $n1 / ($n1 - 1) * 2 * $p1 * (1 - $p1);
			$seg1 += 1;
		}#else{$in_pop1=1;$seg1=1;}
	
		if ($p2 != 1 && $p2 != 0){
			$in_pop2 += $n2 / ($n2 - 1) * 2 * $p2 * (1 - $p2);
			$seg2 += 1;
		}#else{$in_pop2=1;$seg2=1;}
	
	}
	
	my $a1_1; my $a1_2; my $a2_1; my $a2_2;

	for my $ii (1..($n1 -1)){
		$a1_1 += 1 / $ii;
		$a2_1 += 1 / $ii / $ii;
	}

	for my $ii (1..($n2 -1)){
		$a1_2 += 1 / $ii;
		$a2_2 += 1 / $ii / $ii;
	}
	my $b1_1 = ($n1 + 1) / 3 / ($n1 - 1);
	my $b1_2 = ($n2 + 1) / 3 / ($n2 - 1);
	my $b2_1 = 2 * ($n1 * $n1 + $n1 +3) / 9 / $n1 / ($n1 - 1);
	my $b2_2 = 2 * ($n2 * $n2 + $n2 +3) / 9 / $n2 / ($n2 - 1);
	my $c1_1 = $b1_1 - 1 / $a1_1;
	my $c1_2 = $b1_2 - 1 / $a1_2;
	my $c2_1 = $b2_1 - ($n1 + 2) / $a1_1 / $n1 + $a2_1 / $a1_1 / $a1_1;
	my $c2_2 = $b2_2 - ($n2 + 2) / $a1_2 / $n2 + $a2_2 / $a1_2 / $a1_2;
	my $e1_1 = $c1_1 / $a1_1;
	my $e1_2 = $c1_2 / $a1_2;
	my $e2_1 = $c2_1 / ($a1_1 * $a1_1 + $a2_1);
	my $e2_2 = $c2_2 / ($a1_2 * $a1_2 + $a2_2);
	
	my $fst = 1 - ($freq1 + $freq2) / 2 / $freq3;
	my $dxy = $Dxy / $ef_length * 1000;
	my $pi1 = $in_pop1 / $ef_length * 1000;
	my $pi2 = $in_pop2 / $ef_length * 1000;
	my $theta1 = $seg1 / $a1_1 / $ef_length * 1000;
	my $theta2 = $seg2 / $a1_2 / $ef_length * 1000;
	my $D1 = ($in_pop1 - $seg1 / $a1_1) / sqrt($e1_1 * $seg1 + $e2_1 * $seg1 * ($seg1 - 1));
	my $D2 = ($in_pop2 - $seg2 / $a1_2) / sqrt($e1_2 * $seg2 + $e2_2 * $seg2 * ($seg2 - 1));
	
	@pop1 = ();
	@pop2 = ();
	$site=0;
	$l_count = 0;
	$l_seek = $seekpin;
#	print "seek for $l_seek\n";
	print "$fst\t$dxy\t$pi1\t$pi2\t$theta1\t$theta2\t$seg1\t$seg2\t$D1\t$D2\n";

}

}
