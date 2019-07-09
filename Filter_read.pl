#!/perl -w
use strict;
# This script is written by Haomin Lyu
# This script is used to filter reads of Illumina NGS fastq, according to three standards: 
# 1.average quality of every read; 
# 2.number of sites whose quality below 10; 
# 3.number of 'N' in one read.
# When applied to filter genome resequencing data, three standards set as: >20; <=20; <=5.
# Input: pair-end fastq files (two); Output: two files after filter, and two files containing reads filtered.
# 2013-10-23

open FQ1,"<","$ARGV[0]" || die "Cannot open the file: $!";
open FQ2,"<","$ARGV[1]" || die "Cannot open the file: $!";
open FQ3,">>","$ARGV[2]" || die "Cannot open the file: $!";
open FQ4,">>","$ARGV[3]" || die "Cannot open the file: $!";
open FQ5,">>","$ARGV[4]" || die "Cannot open the file: $!";
open FQ6,">>","$ARGV[5]" || die "Cannot open the file: $!";


# Main Programm
my ($name1,$name2,$qua1,$qua2,$seq1,$seq2,$sym1,$sym2);
while ($name1 = <FQ1>, $seq1 = <FQ1>, $sym1 = <FQ1>, $qua1 = <FQ1>,
		$name2 = <FQ2>, $seq2 = <FQ2>, $sym2 = <FQ2>, $qua2 = <FQ2>) {
	$qua1 =~ s/\s+$//;
	$qua2 =~ s/\s+$//;
	$seq1 =~ s/\s+$//;
	$seq2 =~ s/\s+$//;
	my @quality1 = split //,$qua1;
	my @quality2 = split //,$qua2;
	my @read1 = split //,$seq1;
	my @read2 = split //,$seq2;

# Calculate three judgement standards.
	my ($avequa1, $avequa2, $numbelow1, $numbelow2, $numN1, $numN2) = (0,0,0,0,0,0);
	foreach my $site1 (@read1) {
		if ($site1 eq "N") { $numN1 ++; }
	}
	foreach my $site2 (@read2) {
		if ($site2 eq "N") { $numN2 ++; }
	}
	my ($totalqua1, $totalqua2);
	foreach my $siqua1 (@quality1) {
		my $quality = ord($siqua1) - 33;
		$totalqua1 += $quality;
		if ($quality <= 10) { $numbelow1 ++; }
	}
	foreach my $siqua2 (@quality2) {
		my $quality = ord($siqua2) - 33;
		$totalqua2 += $quality;
		if ($quality <= 10) { $numbelow2 ++; }
	}
	$avequa1 = $totalqua1 / @quality1;
	$avequa2 = $totalqua2 / @quality2;


# Output the filtered reads.
	if (($avequa1 > 20 && $numbelow1 <= 20 && $numN1 <= 5) && ($avequa2 > 20 && $numbelow2 <= 20 && $numN2 <= 5)) {
		print FQ3 "$name1$seq1\n$sym1$qua1\n";
		print FQ4 "$name2$seq2\n$sym2$qua2\n";
	} else {
		print FQ5 "$name1$seq1\n$sym1$qua1\n";
                print FQ6 "$name2$seq2\n$sym2$qua2\n";
	}
}


close FQ1;
close FQ2;
close FQ3;
close FQ4;
close FQ5;
close FQ6;


exit;
