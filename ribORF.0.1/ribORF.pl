if ($#ARGV < 2) {
	print "##ribORF predication"."\n";
	print "usage: perl ribORF.pl readCorrectedFile candidateORFFile outputFile [orfLengthCutoff] [orfReadCutoff]"."\n";
	print "readCorrectedFile: input read mapping file after offset correction;"."\n";
	print "candidateORFFile: candidate ORFs, genePred format;"."\n";
	print "outputDir: output directory, with files reporting testing parameters and predicted translating probability;"."\n";
	print "orfLengthCutoff [optional]: cutoff of ORF length (nt), default: 12;"."\n";
	print "orfReadCutoff [optional]: cutoff of supported read numbe, default: 11."."\n";
	exit;
}

my $readfile=$ARGV[0];
my $orfFile=$ARGV[1];
my $outputDir=$ARGV[2];
my $orfLengthCutoff;
my $orfReadCutoff;
if (exists ($ARGV[3])) {
	$orfLengthCutoff=$ARGV[3];
} else {
	$orfLengthCutoff=11;
}
if (exists ($ARGV[4])) {
	$orfReadCutoff=$ARGV[3];
} else {
	$orfReadCutoff=10;
}

my %read;
open (IN, "$readfile");
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\D+/, $s1[5];
	my @s3=split /\d+/, $s1[5];
	if ($s1[1] == 0) {
		my $k=$s1[2].":"."+".":".($s1[3]-1);
		$read{$k}++;
	} elsif ($s1[1] == 16) {
		my $loc=$s1[3];
		for (my $i=0; $i<$#s3; $i++) {
			$loc+=$s2[$i];
		}
		my $k=$s1[2].":"."-".":".($loc-2);
		$read{$k}++;
	}
}
close IN;

open (AN, "$orfFile");
open (OUT, ">$outputDir/input.parameters.txt");
print OUT "geneID"."\t"."transcriptID"."\t"."chrom"."\t"."strand"."\t"."codon5"."\t"."codon3"."\t"."length"."\t"."read.num"."\t"."f1"."\t"."f2"."\t"."f3"."\t"."entropy"."\t"."MAXentropy"."\t"."pme"."\n";
my $dist=3;
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[9];
	my @s3=split /,/, $s1[10];
	my $len1=0;
	my $len2=0;
	my $tot=0;
	my @val;
	my @per;
	my @post;
	for (my $i=0; $i<=2; $i++) {
		$per[$i]=0;
	}
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			push @post, $j;
		}
	}
	if ($s1[3] eq '-') {
		@post=reverse(@post);
	}
	if ($#post >= 5) {
	$len1=$#post+1;
	splice @post, 0, 3;
	#splice @post, -3;
	for (my $m=0; $m<=$#post; $m++) {
		my $k=$s1[2].":".$s1[3].":".$post[$m];
		$len2++;
		if (exists ($read{$k})) {
			$tot+=$read{$k};
			$val[int(($len2-1)/$dist)]+=$read{$k};
			$per[($len2-1)%3]+=$read{$k};
		}
	}
	if ($tot >= $orfReadCutoff && $len2 >= $orfLengthCutoff) {
		my $ent=0;
		my $ten=0;
		my $a=int(($len2+2)/$dist);
		my $b=$tot;
		my $t1=int(($a+$b-1)/$b);
		my @val2;
		for ($i=0; $i<=$#val; $i++) {
			if ($val[$i] > 0) {
				$val2[int($i/$t1)]+=$val[$i];
			}
		}
		for ($i=0; $i<=$#val2; $i++) {
			if ($val2[$i] > 0) {
				my $p=$val2[$i]/($tot);
				$ent+=$p*log(1/$p);
			}
		}
		my $t2=int($a/$t1);
		my $d1=int($b/$t2);
		my $d2=$b%$t2;
		my @va;
		for (my $i=0; $i<$t2; $i++) {
			$va[$i]=$d1;
		}
		for (my $j=0; $j<$d2; $j++) {
			$va[$j]++;
		}
		for (my $i=0; $i<=$#va; $i++) {
			if ($va[$i] > 0) {
				$p=$va[$i]/($tot);
				$ten+=$p*log(1/$p);
			}
		}
		my $per;
		if ($ten == 0) {
			$per=0;
		} else {
			$per=$ent/$ten;
		}
		if ($tot > 0) {
			my $out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$s1[3]."\t".$s1[4]."\t".$s1[5]."\t".$len1."\t".$tot;
			$out.="\t".sprintf("%.3f", $per[0]/$tot);
			$out.="\t".sprintf("%.3f", $per[1]/$tot);
			$out.="\t".sprintf("%.3f", $per[2]/$tot);
			$out.="\t".sprintf("%.3f", $ent);
			$out.="\t".sprintf("%.3f", $ten);
			$out.="\t".sprintf("%.3f", $per);
			print OUT $out."\n";
		}
	}
	}
}
close AN;
close OUT;

open (RI, ">$outputDir/ribORF.learning.R");
print RI "library (\"e1071\")"."\n";
print RI "load (\"ribORF.rda\")"."\n";
print RI "A1 <- read.table (\"$outputDir/input.parameters.txt\", sep=\"\\t\", header=T)"."\n";
print RI "f1 <- A1[,9]"."\n";
print RI "f2 <- A1[,10]"."\n";
print RI "pme <- A1[,14]"."\n";
print RI "A2 <- cbind(f1,f2,pme)"."\n";
print RI "L1 <- sprintf(\"%.3f\", predict(pred1, A2, decision.values = TRUE))"."\n";
print RI "out <- cbind(A1, pvalue=L1)"."\n";
print RI "write.table (out, \"$outputDir/pred.pvalue.parameters.txt\", sep=\"\\t\", quote=F, row.names=F, col.names=T)"."\n";
close RI;

my $com1="Rscript $outputDir/ribORF.learning.R";
system ($com1);
my $com2="rm $outputDir/input.parameters.txt";
#system ($com2);

