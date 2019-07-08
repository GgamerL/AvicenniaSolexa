#!perl -w
#wangzz 20160926
#this script is used for hap(generated by get_hapoem_seg_rec6.pl) file turned into haplotype file step 1 from single 0&1 marks to haplotype marks for all populations needed.
#input all .hap files for needed population
#  
open (OUT,">result")||die "$!";
my @population_list;
for(@ARGV){
	open(IN,$_) or die $!;
	my $name=$_;
	$name=~s/\.hap//;
	push @population_list,$name;
	while(<IN>){
		s/\s+$//;
		($gene,$vxa,$vxb,$vxc,$vxd,$vxe,$vxf,$vxg,$vxh,$seg,$count,$site,$type,@hap)=split/\t/;
		@s=split/,/,$site;
		$n=@hap;
		$range="$s[0]<$s[-1]";
		for(1..$n){
			my $string;
			($alle,$freq)=split/,/,$hap[$_-1];
			@a=split//,$alle;
			for(0..($count-1)){
				$se=$s[$_]*$a[$_];
				if($se!=0){
					$string.="$se,";
				}
			}
			if(defined ($string) && $string=~/,/){
				$add="$name\t$seg\t$range\t$string\t$freq\n";
				push @{$gene_list{$gene}},$add;
			}else{
				if($freq!=1){
					$add="$name\t$seg\t$range\t,\t$freq\n";
					push @{$gene_list{$gene}},$add;
				}
			}
		}
	}
	close IN;
}
my %type;
foreach $gene(sort keys %gene_list){
	my %site_list;
	my %type_count;
	my @all;
	my %seg;
	my @cut;
	my %saw2;
	my %saw3;
	for(0..scalar@{$gene_list{$gene}}-1){
		($p,$seg,$site)=(split/\t/,@{$gene_list{$gene}}[$_])[0,1,3];
		@s=split/,/,$site;
		$type_count{"$p\t$seg"}+=1;
		push @{$site_list{"$p\t$seg"}},@s;
	}
	for(sort keys %site_list){
		my %saw;
		$popu=(split/\t/)[0];
		$seg{$popu}+=1;
		@saw{@{$site_list{$_}}}=();
		@{$site_list{$_}}=sort{$a<=>$b} keys %saw;
		push @all,@{$site_list{$_}};
		if($seg{$popu}!=1){
			push @cut,@{$site_list{$_}}[0];
		}
	}
	@saw2{@cut}=();
	@cut=sort {$a<=>$b} keys %saw2;
	
	$duan=@cut;
	@saw3{@all}=();
	@all=sort {$a<=>$b} keys %saw3;
	$segment{$gene}=1;
	foreach $i(1..$duan){
		my $string;
		$cut_off=$cut[$i-1];
		for(@all){
			if($_<$cut_off){
			$string.="$_,";
			push @{"$gene$i"},$_;
		}	
	}
	push @{$type{$gene}},$string;
	for(@{"$gene$i"}){
	shift(@all);
	}
	$segment{$gene}+=1;
	}
	my $string2;
	for(@all){
		$string2.="$_,";
	}
	push @{$type{$gene}},$string2;
}
foreach $gen(sort keys %type){
	$n=@{$type{$gen}};
	for(1..$n){
	$ss{"$gen\t$_"}=${$type{$gen}}[$_-1];
	@{$_}=split/,/,${$type{$gen}}[$_-1];
	}
	for(@{$gene_list{$gen}}){
	chomp;
	($po,$seg,$range,$site,$freq)=split/\t/;
	($min,$max)=split/\</,$range;
	@si=split/,/,$site;
	foreach $i(1..$n){
	my $haplotype;
	if(!($max<${$i}[0] || $min> ${$i}[-1])){
		for(@{$i}){
			if($_~~@si){
				$haplotype.="1";
			}else{
				$haplotype.="0";
			}
		}
		$add2="$po\t$haplotype,$freq";
		push @{$hap_list{"$gen\t$i"}},$add2;
	}
	}
}
}
for(sort keys %hap_list){
	print OUT "$_\t$ss{$_}\n";
	my @site_num=split/,/,$ss{$_};
	my $num=@site_num;
	my %hap;
	my %haplotype;
	for(@{$hap_list{$_}}){
			chomp;
			($po,$hap)=split/\t/;
			($type,$freq)=split/,/,$hap;
			$hap{"$po\t$type"}+=$freq;
	}
	for(sort keys %hap){
		chomp;
		($population,$haplotype)=split/\t/;
		$f=$hap{$_};
		$haplotype{$population}.="\t$haplotype,$f";
	}
	for(@population_list){
		if(exists $haplotype{$_}){
			print OUT "$_$haplotype{$_}\n";
		}else{
			my $bu="0"x$num;
			print OUT "$_\t$bu,1\n";
		}
	}
	print OUT "\/\/\n";
}