#!perl -w
#originaly written by Ziwen He,2013
#revised by Zhengzhen Wang, 2015
#input SNP frequency matrix and mapview file
($error,$rec,$overlapreads,$lowest,$seglength,$snpfile,$mapviewfile)=@ARGV[0..6];
$noerror=1-$error;
open F,"<$snpfile";
while(<F>){
	($ch,$po,$a1,$a2)=(split/\s/)[0,1,2,3];
	$ch{"$ch"}.="$po,";
	$a1{"$ch:$po"}=$a1;
	$a2{"$ch:$po"}=$a2
}
open F,"<$mapviewfile";
while(<F>){
	($read,$ch,$po,$flag,$indel,$readlength,$seq)=(split/\t/)[0,1,2,5,6,13,14];
	next unless(exists $ch{$ch} && ($flag<130 || ($flag==130 && $indel==0)));
# next unless(exists $ch{$ch} && $flag==18);
	$seq!~tr/a-z/A-Z/;
	@snp=(split/,/,$ch{$ch});
	$flag{$flag}++;
	$type="";
	for(@snp){
		$snp=$_;
		if($snp>=$po && $snp<$po+$readlength){
			$base=substr($seq,$snp-$po,1);
			if($base eq $a1{"$ch:$snp"}){
				$type.="0"
			}elsif($base eq $a2{"$ch:$snp"}){
				$type.="1"
			}else{
				$type.="."
			}
		}else{
			$type.="."
		}
	}
	$allsnp=@snp;
	next if($type eq "."x$allsnp);
	($read,$pair)=(split/\//,$read)[0,1];
	if($pair eq "1"){
		$allright{$read}=$type
	}else{
		$allleft{$read}=$type
	}
	$allread{$read}=1;
	$chr=$ch;
	last
}
while(<F>){
	$line=$_;
	($read,$ch,$po,$flag,$indel,$readlength,$seq)=(split/\t/)[0,1,2,5,6,13,14];
	next unless(exists $ch{$ch} && ($flag<130 || ($flag==130 && $indel==0)));
# next unless(exists $ch{$ch} && $flag==18);
	$seq!~tr/a-z/A-Z/;
	if($ch eq $chr){
		$flag{$flag}++;
		$type="";
  	for(@snp){
  		$snp=$_;
  		if($snp>=$po && $snp<$po+$readlength){
  			$base=substr($seq,$snp-$po,1);
  			if($base eq $a1{"$ch:$snp"}){
  				$type.="0"
  			}elsif($base eq $a2{"$ch:$snp"}){
  				$type.="1"
  			}else{
  				$type.="."
  			}
  		}else{
  			$type.="."
  		}
  	}
  	next if($type eq "."x$allsnp);
  	($read,$pair)=(split/\//,$read)[0,1];
  	if($pair eq "1"){
  		$allright{$read}=$type
  	}else{
  		$allleft{$read}=$type
  	}
  	$allread{$read}=1
  }else{
  	&HAPO;
  	%flag=();
  	%alltype=();
  	$type="";
		($read,$ch,$po,$flag,$indel,$readlength,$seq)=(split/\t/,$line)[0,1,2,5,6,13,14];
		$chr=$ch;
		$seq!~tr/a-z/A-Z/;
		@snp=(split/,/,$ch{$ch});
		$flag{$flag}++;
  	for(@snp){
  		$snp=$_;
  		if($snp>=$po && $snp<$po+$readlength){
  			$base=substr($seq,$snp-$po,1);
  			if($base eq $a1{"$ch:$snp"}){
  				$type.="0"
  			}elsif($base eq $a2{"$ch:$snp"}){
  				$type.="1"
  			}else{
  				$type.="."
  			}
  		}else{
  			$type.="."
  		}
  	}
  	$allsnp=@snp;
  	next if($type eq "."x$allsnp);
  	($read,$pair)=(split/\//,$read)[0,1];
  	if($pair eq "1"){
  		$allright{$read}=$type
  	}else{
  		$allleft{$read}=$type
  	}
  	$allread{$read}=1
  }
}
&HAPO;

sub HAPO{
	%break=();
	for(1..($allsnp-1)){
		$snp=$_-1;
		$break{$snp}=0
	}
	for(keys %allread){
		$read=$_;
		$type="";
		if(exists $allright{$read}){
			if(exists $allleft{$read}){
				for(1..$allsnp){
					$snp=$_-1;
					$rightsnp=substr($allright{$read},$snp,1);
					$leftsnp=substr($allleft{$read},$snp,1);
					if($rightsnp eq $leftsnp){
						$type.=$rightsnp
					}elsif($rightsnp eq "."){
						$type.=$leftsnp
					}elsif($leftsnp eq "."){
						$type.=$rightsnp
					}else{
						$type.="."
					}
				}
			}else{
				$type=$allright{$read}
			}
		}else{
			$type=$allleft{$read}
		}
		$alltype{$type}++;
		for(keys %break){
			$break=$_;
			unless(substr($type,$break,1) eq "." || substr($type,$break+1,1) eq "."){
				$break{$break}++
			}
		}
	}
	%allright=();
	%allleft=();
	%allread=();
	$piece=1;
	$s=0;
	$le=1;
	for(1..($allsnp-1)){
		$i=$_;
		if($break{$i-1}>=$overlapreads){
			$le++
		}else{
			$nexts=$i;
			&PP;
			$s=$nexts;
			$le=1;
			$piece++
		}
	}
	&PP;
}

sub PP{
	%type=();
	$snp=$le;
	print "$chr\t$flag{1}\t$flag{2}\t$flag{4}\t$flag{8}\t$flag{18}\t$flag{32}\t$flag{64}\t$flag{130}\t$piece\t$snp\t";
	for(1..$le){
		$i=$s+$_-1;
		print "$snp[$i],"
	}
	for(keys %alltype){
		$type=$_;
		$subtype=substr($type,$s,$le);
		if($subtype ne "."x$le){
			$type{$subtype}+=$alltype{$type}
		}
	}
	$max=2**$snp-1;
	$numseg=int(($snp-1)/$seglength)+1;
	$seg=0;
	$segsnp=$seglength;
	$segsnp=$snp if($numseg==1);
	&Singleseg;
	if($numseg==1){
		$numberofhapo=keys %c;
		print "\t$numberofhapo";
		for(keys %c){
			print "\t$_,$c{$_}"
		}
		print "\n";
		return
	}else{
		%seg1=%c
	}
	for(1..($numseg-1)){
		$seg=$_;
		$segsnp=$seglength;
		$segsnp=$snp-($numseg-1)*$seglength if($seg==($numseg-1));
		&Singleseg;
		%seg2=%c;
		$segsnp=($seg+1)*$seglength;
		$segsnp=$snp if($seg==($numseg-1));
		&Dualseg;
		%seg1=%c
	}
	$numberofhapo=keys %seg1;
	print "\t$numberofhapo";
	for(keys %seg1){
		print "\t$_,$seg1{$_}"
	}
	print "\n"	
}

sub Singleseg{
	%p1=();
	%p2=();
	%c=();
	%segtype=();
	$nsegtype=0;
	for(keys %type){
		$type=$_;
		$segtype=substr($type,$seglength*$seg,$segsnp);
		if($segtype ne "."x$segsnp){
			$nsegtype+=2*$type{$type};
			if($segsnp==1){
				$segtype{"$segtype\t0\t0"}+=$type{$type}
			}else{
				for(0..($segsnp-1)){
					$j=$_;
					$typesnp=substr($segtype,$j,1);
					if($typesnp ne "."){
						$reads=$j;
						last
					}
				}
				for(0..($segsnp-1)){
					$j=$segsnp-1-$_;
					$typesnp=substr($segtype,$j,1);
					if($typesnp ne "."){
						$reade=$j;
						last
					}
				}
				if($reads==$reade){
#					$nsegtype-=2*$type{$type};
 					$segtype{"$segtype\t$reads\t$reade"}+=$type{$type}
				}else{
					$allrecp="";
					for($reads..($reade-1)){
  					$j=$_;
  					$i=$seglength*$seg+$j;
  					$recp=(1-(1-$rec)**($snp[$s+$i+1]-$snp[$s+$i]))*(1-$rec)**($snp[$s+$seglength*$seg+$reade]-$snp[$s+$seglength*$seg+$reads]-$snp[$s+$i+1]+$snp[$s+$i]);
  					$allrecp.="$recp,"
  				}
  				$recp=(1-$rec)**($snp[$s+$seglength*$seg+$reade]-$snp[$s+$seglength*$seg+$reads]);
  				$allrecp.="$recp";
  				$segtype{"$segtype\t$reads\t$reade\t$allrecp"}+=$type{$type}
  			}
  		}
  	}
  }
	$segmax=2**$segsnp-1;
	for(0..$segmax){
		$hapo=sprintf("%b",$_);
		$x=$segsnp-length($hapo);
		$hapo="0"x$x.$hapo;
		$c{$hapo}=1/($segmax+1);
		for(keys %segtype){
			$type=$_;			
			($type,$reads,$reade,$allrecp)=(split/\t/,$type)[0..3];
			@sitep=();
			for($reads..$reade){
				$i=$_;
				$typesnp=substr($type,$i,1);
				$haposnp=substr($hapo,$i,1);
				if($typesnp eq "."){
					$sitep[$i]=1
				}elsif($typesnp eq $haposnp){
					$sitep[$i]=$noerror
				}else{
					$sitep[$i]=$error
				}
			}
			if($reads==$reade){
				$p1{"$hapo\t$type"}=$sitep[$reads];
				$p2{"$hapo\t$type"}=$sitep[$reads]
			}else{
				@allrecp=(split/,/,$allrecp);
				$j=-1;
				$recp=$allrecp[$j];
				$allsitep=1;
				for($reads..$reade){
					$i=$_;
					$allsitep*=$sitep[$i]
				}
				$p1{"$hapo\t$type"}+=$recp*$allsitep;
				$p2{"$hapo\t$type"}+=$recp*$allsitep;
				for($reads..($reade-1)){
  				$i=$_;
  				$j++;
  				$recp=$allrecp[$j];
  				$allsitep=1;
  				for($reads..$i){
  					$allsitep*=$sitep[$_]
  				}
  				$p1{"$hapo\t$type"}+=$recp*$allsitep;
  				$allsitep=1;
  				for(($i+1)..$reade){
  					$allsitep*=$sitep[$_]
  				}
  				$p2{"$hapo\t$type"}+=$recp*$allsitep
				}
			}
		}
	}	
	for(1..10000){
		$rep=$_;
		%pc=();
		for(keys %segtype){
			$key=$_;
			$type=(split/\t/,$key)[0];
			$total1=0;
			$total2=0;
			for(keys %c){
				$hapo=$_;
				$total1+=$p1{"$hapo\t$type"}*$c{$hapo};
				$total2+=$p2{"$hapo\t$type"}*$c{$hapo};
			}
			for(keys %c){
				$hapo=$_;
				$pc{$hapo}+=$p1{"$hapo\t$type"}*$c{$hapo}/$total1*$segtype{$key};
				$pc{$hapo}+=$p2{"$hapo\t$type"}*$c{$hapo}/$total2*$segtype{$key};
			}
		}
		if($rep<=100){
  		for(keys %pc){
  			$hapo=$_;
  			$c{$hapo}=$pc{$hapo}/$nsegtype;
  			delete $c{$hapo} if($c{$hapo}<$lowest/($segmax+1))
  		}
  	}elsif($rep%50==0){  
  		$checksame=0;		
  		for(keys %pc){
  			$hapo=$_;
  			unless($c{$hapo}==$pc{$hapo}/$nsegtype){
  				$checksame=1;
  				last
  			}
  		}
  		if($checksame==1){
  			for(keys %pc){
  				$hapo=$_;
  				$c{$hapo}=$pc{$hapo}/$nsegtype;
  				delete $c{$hapo} if($c{$hapo}<$lowest)
  			}
  		}else{
  			last
  		}
  	}else{
		  for(keys %pc){
				$hapo=$_;
				$c{$hapo}=$pc{$hapo}/$nsegtype;
				delete $c{$hapo} if($c{$hapo}<$lowest)
			}
		}
	}
}
		
sub Dualseg{
	%p1=();
	%p2=();
	%c=();
	%segtype=();
	$nsegtype=0;
	for(keys %type){
		$type=$_;
		$segtype=substr($type,0,$segsnp);
		if($segtype ne "."x$segsnp){
			$nsegtype+=2*$type{$type};
			if($segsnp==1){
				$segtype{"$segtype\t0\t0"}+=$type{$type}
			}else{
				for(0..($segsnp-1)){
					$j=$_;
					$typesnp=substr($segtype,$j,1);
					if($typesnp ne "."){
						$reads=$j;
						last
					}
				}
				for(0..($segsnp-1)){
					$j=$segsnp-1-$_;
					$typesnp=substr($segtype,$j,1);
					if($typesnp ne "."){
						$reade=$j;
						last
					}
				}
				if($reads==$reade){
#					$nsegtype-=2*$type{$type};
					$segtype{"$segtype\t$reads\t$reade"}+=$type{$type}
				}else{
					$allrecp="";
					for($reads..($reade-1)){
  					$i=$_;
  					$recp=(1-(1-$rec)**($snp[$s+$i+1]-$snp[$s+$i]))*(1-$rec)**($snp[$s+$reade]-$snp[$s+$reads]-$snp[$s+$i+1]+$snp[$s+$i]);
  					$allrecp.="$recp,"
  				}
  				$recp=(1-$rec)**($snp[$s+$reade]-$snp[$s+$reads]);
  				$allrecp.="$recp";
  				$segtype{"$segtype\t$reads\t$reade\t$allrecp"}+=$type{$type}
  			}
  		}
  	}
  }
	$segmax=2**$segsnp-1;
	for(keys %seg1){
		$hapo1=$_;
		for(keys %seg2){
			$hapo2=$_;
			$hapo="$hapo1$hapo2";
			$c{$hapo}=1/($segmax+1);
  		for(keys %segtype){
  			$type=$_;
  			($type,$reads,$reade,$allrecp)=(split/\t/,$type)[0..3];
  			@sitep=();
  			for($reads..$reade){
  				$i=$_;
  				$typesnp=substr($type,$i,1);
  				$haposnp=substr($hapo,$i,1);
  				if($typesnp eq "."){
  					$sitep[$i]=1
  				}elsif($typesnp eq $haposnp){
  					$sitep[$i]=$noerror
  				}else{
  					$sitep[$i]=$error
  				}
  			}
  			if($reads==$reade){
  				$p1{"$hapo\t$type"}=$sitep[$reads];
  				$p2{"$hapo\t$type"}=$sitep[$reads]
  			}else{
  				@allrecp=(split/,/,$allrecp);
  				$j=-1;
  				$recp=$allrecp[$j];
  				$allsitep=1;
  				for($reads..$reade){
  					$i=$_;
  					$allsitep*=$sitep[$i]
  				}
  				$p1{"$hapo\t$type"}+=$recp*$allsitep;
  				$p2{"$hapo\t$type"}+=$recp*$allsitep;
  				for($reads..($reade-1)){
    				$i=$_;
    				$j++;
    				$recp=$allrecp[$j];
    				$allsitep=1;
    				for($reads..$i){
    					$allsitep*=$sitep[$_]
    				}
    				$p1{"$hapo\t$type"}+=$recp*$allsitep;
    				$allsitep=1;
    				for(($i+1)..$reade){
    					$allsitep*=$sitep[$_]
    				}
    				$p2{"$hapo\t$type"}+=$recp*$allsitep
  				}
  			}
  		}
  	}
	}
	for(1..10000){
		$rep=$_;
		%pc=();
		for(keys %segtype){
			$key=$_;
			$type=(split/\t/,$key)[0];
			$total1=0;
			$total2=0;
			for(keys %c){
				$hapo=$_;
				$total1+=$p1{"$hapo\t$type"}*$c{$hapo};
				$total2+=$p2{"$hapo\t$type"}*$c{$hapo};
			}
			for(keys %c){
				$hapo=$_;
				$pc{$hapo}+=$p1{"$hapo\t$type"}*$c{$hapo}/$total1*$segtype{$key};
				$pc{$hapo}+=$p2{"$hapo\t$type"}*$c{$hapo}/$total2*$segtype{$key};
			}
		}
		if($rep<=100){
  		for(keys %pc){
  			$hapo=$_;
  			$c{$hapo}=$pc{$hapo}/$nsegtype;
  			delete $c{$hapo} if($c{$hapo}<$lowest/($segmax+1))
  		}
  	}elsif($rep%50==0){
  		$checksame=0;
  		for(keys %pc){
  			$hapo=$_;
  			unless($c{$hapo}==$pc{$hapo}/$nsegtype){
  				$checksame=1;
  				last
  			}
  		}
  		if($checksame==1){
  			for(keys %pc){
  				$hapo=$_;
  				$c{$hapo}=$pc{$hapo}/$nsegtype;
  				delete $c{$hapo} if($c{$hapo}<$lowest)
  			}
  		}else{
  			last
  		}
  	}else{
		  for(keys %pc){
				$hapo=$_;
				$c{$hapo}=$pc{$hapo}/$nsegtype;
				delete $c{$hapo} if($c{$hapo}<$lowest)
			}
		}
	}
}
