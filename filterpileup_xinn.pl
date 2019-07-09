#! perl -w 
# This script is written by Xinnian Li
# This script filters the specific-site in a read which have mapped on the reference. The pileup file should contain the position 
# information.

#main program
for (<>){
@line=split(/\t/);
$fline=join "\t",$line[0],$line[1],$line[2];
print $fline,"\t";
&fliter($line[4], $line[5], $line[7], $line[3]);  # input base, quality, position in the read, depth
print $num,"\t","@",@rest,"\t","@",@rqust,"\t",join(',',@position),"\n";
}

sub fliter{
$num=0;
@rest="";
@rqust="";
@position="";
@b=split //,$_[0];
@q=split //,$_[1];
@p=split /,/,$_[2];
for $n(1..$_[3]){
if (ord($q[$n])-64>=20 && $p[$n-1]>=11 && $p[$n-1]<=80
){
push @rest,$b[$n];
push @rqust,$q[$n]; 
push @position,$p[$n-1];
$num++;
}
}
return $num,@rest,@rqust,@position;
}
