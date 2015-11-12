#!/usr/bin/perl
die "usage: perl smith-waterman_aln_two_sequence.pl <sequence 1> <sequence 2>\n" unless @ARGV == 2;

#set penalty
$g=-2;
$match=2;
$mismatch=-1;

alignment($ARGV[0],$ARGV[1]);

sub alignment {

	$S1=$_[0];
	$S2=$_[1];
	$s1_length=length($S1);
	$s2_length=length($S2);
	
	@s1=[];
	@s2=[];
	@A=[];
	
	#store align result
	@r1=[];
	@r2=[];
	
	$counter=0;
	while($s1[$counter]=substr($S1,$counter,1)){$counter=$counter+1};
	$counter=0;
	while($s2[$counter]=substr($S2,$counter,1)){$counter=$counter+1};
	
	
	#initiate align score
	for $i(0..$s1_length){$A[0][$i]=$i*$g};
	for $i(0..$s2_length){$A[$i][0]=$i*$g};
	#calculate align score
	for $j(1..($s2_length)){
		for $i(1..($s1_length)){
			$A[$j][$i]=max(
				$A[$j-1][$i]+$g,
				$A[$j][$i-1]+$g,
				$A[$j-1][$i-1]+singleMatch($s2[$j-1],$s1[$i-1])
			);
		};
	};
	
	#print align score
	#for $j(0..$s2_length){
	#	for $i(0..$s1_length){
	#		printA($A[$j][$i]);
	#		if ($i==$s1_length){print "\n";};
	#	};
	#};
	
	#get align result
	$i=$s1_length;
	$j=$s2_length;
	$k1=0;
	$k2=0;
	
	
	if ($A[$j][$i]>=$A[$j-1][$i] and $A[$j][$i]>=$A[$j][$i-1]) { $r1[$k1]=$s1[$i-1];$r2[$k2]=$s2[$j-1];}
	elsif ($A[$j-1][$i]>=$A[$j][$i] and $A[$j-1][$i]>=$A[$j][$i-1]) { $r1[$k1]=$s1[$i-1];$r2[$k2]=$s2[$j-2]; $j=$j-1;}
	elsif ($A[$j][$i-1]>=$A[$j][$i] and $A[$j][$i]>=$A[$j-1][$i]) { $r1[$k1]=$s1[$i-2];$r2[$k2]=$s2[$j-1]; $i=$i-1;}
	while( $i>1 and $j>1 ){
		$k1=$k1+1;$k2=$k2+1;
		if ( $A[$j-1][$i-1]>=$A[$j-1][$i] and $A[$j-1][$i-1]>=$A[$j][$i-1] ) { $i=$i-1;$j=$j-1; $r1[$k1]=$s1[$i-1]; $r2[$k2]=$s2[$j-1];}
		elsif ( $A[$j-1][$i]>=$A[$j-1][$i-1] and $A[$j-1][$i]>=$A[$j][$i-1] ) { $j=$j-1; $r1[$k1]="-"; $r2[$k2]=$s2[$j-1];}
		elsif ( $A[$j][$i-1]>=$A[$j-1][$i-1] and $A[$j][$i-1]>=$A[$j-1][$i] ) { $i=$i-1; $r2[$k2]="-"; $r1[$k1]=$s1[$i-1];}
	}
	
	#if ($A[$j][$i]>=$A[$j-1][$i] and $A[$j][$i]>=$A[$j][$i-1]) { $r1[$k1]=$s1[$i-1];$r2[$k2]=$s2[$j-1];print "($j,$i)\t($r2[$k2],$r1[$k1])\n"; }
	#elsif ($A[$j-1][$i]>=$A[$j][$i] and $A[$j-1][$i]>=$A[$j][$i-1]) { $r1[$k1]=$s1[$i-1];$r2[$k2]=$s2[$j-2];print "(".($j-1).",$i)\t($r2[$k2],$r1[$k1])\n"; $j=$j-1;}
	#elsif ($A[$j][$i-1]>=$A[$j][$i] and $A[$j][$i]>=$A[$j-1][$i]) { $r1[$k1]=$s1[$i-2];$r2[$k2]=$s2[$j-1];print "($j,".($i-1).")\t($r2[$k2],$r1[$k1])\n"; $i=$i-1;}
	#while( $i>1 and $j>1 ){
	#	$k1=$k1+1;$k2=$k2+1;
	#	if ( $A[$j-1][$i-1]>=$A[$j-1][$i] and $A[$j-1][$i-1]>=$A[$j][$i-1] ) { $i=$i-1;$j=$j-1; $r1[$k1]=$s1[$i-1]; $r2[$k2]=$s2[$j-1]; print "($j,$i)\t($r2[$k2],$r1[$k1])\n"}	
	#	elsif ( $A[$j-1][$i]>=$A[$j-1][$i-1] and $A[$j-1][$i]>=$A[$j][$i-1] ) { $j=$j-1; $r1[$k1]="-"; $r2[$k2]=$s2[$j-1]; print "($j,$i)\t($r2[$k2],$r1[$k1])\n"}
	#	elsif ( $A[$j][$i-1]>=$A[$j-1][$i-1] and $A[$j][$i-1]>=$A[$j-1][$i] ) { $i=$i-1; $r2[$k2]="-"; $r1[$k1]=$s1[$i-1]; print "($j,$i)\t($r2[$k2],$r1[$k1])\n"}
	#}
	
	$len1=@r1;
	$len2=@r2;
	for $i(1..$len1){print $r1[$len1-$i];};
	print "\n";
	for $j(1..$len2){print $r2[$len2-$j];};
	print "\n";

	$identity_base=0;
	for $i(1..$len1){
		if ($r1[$len1-$i] eq $r2[$len2-$i]){
			$identity_base=$identity_base+1;
		} 
	};
	$identity=$identity_base*100/$len1;
	print "Identity score: $identity%\n\n";

	return $identity;
}#alignment sub end


#sub routine to align single letter
sub singleMatch{
	if ($_[0] eq $_[1]){return $match;}
	else {return $mismatch;}
}

sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}



