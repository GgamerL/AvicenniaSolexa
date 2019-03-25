#!perl -w
use strict;
use List::Util qw/sum/;

# For calculating the Fst among different populations with the segmental resequencing SNP data
# Input:README (sample size); Depth files; All SNP data files with prefix or full names
# Output: table of Fst
# By Haomin Lyu
# 2016-9-11
# wangzz, 2016-10-12
my $readmefile="/data3/users/haomin/PopSeg/index/README";
my @pops = ('AmNZ', 'AmBS', 'AmCA', 'AmDW', 'AmBB', 'AmBL', 'AmCB', 'AmSB', 'AmWC', 'AmSY', 'AmSS', 'AmBK', 'AmTN', 'AmLS', 'AmPN', 'AmKY');

my %sample;
#my $readmefile = shift @ARGV;
open README,"<","$readmefile"
	or die "Cannot open the $readmefile: $!";
while (<README>) {
	chomp;
	my @readmeline = split /\t/;
	$sample{$readmeline[0]} = $readmeline[3];
}
close README;


my %covsite;
my %effective_len;
foreach my $popgenefile (@pops) {
	$popgenefile="/data3/users/haomin/PopSeg/Depth/"."$popgenefile"."\.depth";
	open Cov,"<","$popgenefile"
		or die "Cannot open the $popgenefile: $!";

	$popgenefile =~ s/.depth//;
	$popgenefile =~ s/.+\///;

	while (<Cov>) {
		chomp;
		my @geneline = split /\t/;
		#my $scapos = join("\t", $geneline[0], $geneline[1]);
		${${$covsite{$popgenefile}}{$geneline[0]}}{$geneline[1]} = 1;
		if($effective_len{$geneline[0]}){
			$effective_len{$geneline[0]}++;
		}else{
			$effective_len{$geneline[0]}=1;
		}
	}

	close Cov;
}


my %snp;
foreach my $popsnpfile (@pops) {
	$popsnpfile="/data2/users/wangzz/SOLEXA/A_marina/haplotype8/"."$popsnpfile"."\.snp\.pro";
	open FILE,"<","$popsnpfile"
		or die "Cannot open the $popsnpfile: $!";

	$popsnpfile =~ s/.snp.pro//;
	$popsnpfile =~ s/.+\///;

	while (<FILE>) {
		chomp;
		my @snpline = split /\t/;
		#my $snplabel = join("\t", $snpline[0], $snpline[1]);
		${${$snp{$popsnpfile}}{$snpline[0]}}{$snpline[1]} = $snpline[7];
	}

	close FILE;
}


#my @pops = ('AmNZ', 'AmBS', 'AmCA', 'AmDW', 'AmBB', 'AmBL', 'AmCB', 'AmSB', 'AmWC', 'AmSY', 'AmSS', 'AmBK', 'AmTN', 'AmLS', 'AmPN', 'AmKY'); 
my $popno = scalar @pops;
for my $popno1 (1..$popno) {
	my $pop1 = $pops[$popno1 - 1];

	for my $popno2 (($popno1+1)..$popno) {
		my $pop2 = $pops[$popno2 - 1];
		open (OUT,">$pop1\_$pop2")||die "$!";
		my $sample1 = $sample{$pop1} * 2;
		#print "$sample1";
		my $sample2 = $sample{$pop2} * 2;
#		print OUT "a";=sort {$a<=>$b} keys %saw2
		my @pp1=keys %{$covsite{$pop1}};
		my @pp2=keys %{$covsite{$pop2}};
		my %hash0;
		my @pp=(@pp1,@pp2);
		 @pp = grep{ ++$hash0{$_} < 2} @pp;
		foreach my $gen_key(@pp){
			next unless exists ${$covsite{$pop1}}{$gen_key} && ${$covsite{$pop2}}{$gen_key};
			#print "yes,$gen_key\n";
			my @covname1 = keys %{${$covsite{$pop1}}{$gen_key}};
			my @covname2 = keys %{${$covsite{$pop2}}{$gen_key}};
			#print "@covname1\n\n";
			my @covname = (@covname1, @covname2);
			my %hash1;
			@covname = grep { ++$hash1{$_} < 2 } @covname; 

			my $length = scalar @covname;

			my @snpname1 = keys %{${$snp{$pop1}}{$gen_key}};
			my @snpname2 = keys %{${$snp{$pop2}}{$gen_key}};
			my @snpname = (@snpname1, @snpname2);
			#print "@snpname\n\n";
			my %hash2;
			@snpname = grep { ++$hash2{$_} < 2 } @snpname;
	
			my ($freq1, $freq2, $freq3);
			foreach my $snpname (@snpname) {
				next unless exists ${${$covsite{$pop1}}{$gen_key}}{$snpname} && ${${$covsite{$pop2}}{$gen_key}}{$snpname};

				my ($p1, $p2, $p3);
				if (${${$snp{$pop1}}{$gen_key}}{$snpname}) { 
					$p1 = ${${$snp{$pop1}}{$gen_key}}{$snpname};
					#print "$gen_key\t$snpname\t$p1\n";
				} else { 
					$p1 = 0;
				}

				if (${${$snp{$pop2}}{$gen_key}}{$snpname}) {
					$p2 = ${${$snp{$pop2}}{$gen_key}}{$snpname};
				} else {
					$p2 = 0;
				}
				next if ($p1==$p2);
				$p3 = ($p1 * $sample1 + $p2 * $sample2) / ($sample1 + $sample2);
	
				$freq1 += $sample1 / ($sample1 - 1) * $p1 * (1 - $p1);
				$freq2 += $sample2 / ($sample2 - 1) * $p2 * (1 - $p2);
				$freq3 += ($sample1 + $sample2) / ($sample1 + $sample2 - 1) * $p3 * (1 - $p3);
				if ($freq3 == 0){
					print "$pop1\t$pop2\t$gen_key\t$snpname\t$sample1\t$sample2\t$p1\t$p2\t$p3\n";
					next;
				}
			}

			#my $fst = 1 - ($freq1 + $freq2)/2/$freq3;
			#print OUT "$gen_key\t$effective_len{$gen_key}\t$fst\n";
		}
	close OUT;
	}
	
}


exit;
