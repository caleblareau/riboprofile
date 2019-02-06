if ($#ARGV < 3) {
	print "## Calcualte read distribution around CDS start and stop codons"."\n";
	print "usage: perl readDist.pl readFile geneFile outputDir readLength [leftNum] [rightNum]"."\n";
	print "readFile: read mapping file, SAM format;"."\n";
	print "geneFile: canonical protein-coding ORF annotation, genepred format;"."\n";
	print "outputDir: output directory;"."\n";
	print "readLength: specified RPF length;"."\n";
	print "leftNum [optional]: N nucleotides upstream start codon and downstream stop codon, default: 30;"."\n";
	print "rightNum [optional]: N nucleotides downstream start codon and upstream stop codon, default: 50."."\n";
	exit;
}

my $readfile=$ARGV[0];
my $geneFile=$ARGV[1];
my $outputDir=$ARGV[2];
my $readLength=$ARGV[3];
my $start;
my $end;
if ($ARGV[4]) {
	$start=$ARGV[4];
} else {
	$start=30;
}
if ($ARGV[5]) {
	$end=$ARGV[5];
} else {
	$end=50;
}

open (IN, "$readfile");
my %read;
my $tot=0;
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\D+/, $s1[5];
	my @s3=split /\d+/, $s1[5];
	if (length($s1[9]) == $readLength) {
	if ($s1[1] == 0) {
		my $k=$s1[2].":"."+".":".($s1[3]-1);
		$read{$k}++;
		$tot++;
	} elsif ($s1[1] == 16) {
		my $loc=$s1[3];
		for (my $i=0; $i<$#s3; $i++) {
			$loc+=$s2[$i];
		}
		my $k=$s1[2].":"."-".":".($loc-2);
		$read{$k}++;
		$tot++;
	}
	}
}
close IN;
	
my %count5;
my %count3;
open (AN, "$geneFile");
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[$#s1-1];
	my @s3=split /,/, $s1[$#s1];
	my @val;
	my $k=-1;
	my $codon5;
	my $codon3;
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			$k++;
			push @val, $j;
			if ($s1[3] eq '+') {
				if ($j==$s1[6]) {
					$codon5=$k;
				} elsif ($j==($s1[7]-3)) {
					$codon3=$k;
				}
			} else {
				if ($j==($s1[6]+2)) {
					$codon5=$k;
				} elsif ($j==($s1[7]-1)) {
					$codon3=$k;
				}
			} 
		}
	}
	if ($codon5 > $start && ($codon3-$codon5) > $end && ($#val-$codon3) > $start) {
		$count5{"0"}++;
		$count3{"0"}++;
		for (my $i=-$start; $i<=$end; $i++) {
			my $k;
			if ($s1[3] eq '+') {
				$k=$s1[2].":".$s1[3].":".$val[$codon5+$i];
			} else {
				$k=$s1[2].":".$s1[3].":".$val[$codon3-$i];
			}
			if (exists ($read{$k})) {
				$count5{$i+$start+1}+=$read{$k};
			}
		} 
		for (my $i=-$end; $i<=$start; $i++) {
			my $k;
			if ($s1[3] eq '+') {
				$k=$s1[2].":".$s1[3].":".$val[$codon3+$i];
			} else {
				$k=$s1[2].":".$s1[3].":".$val[$codon5-$i];
			}
			if (exists ($read{$k})) {
				$count3{$i+$end+1}+=$read{$k};
			}
		} 
	}
}
close AN;

if (! exists ($count5{"0"})) {
	print "Error: no hits!"."\n";
	exit;
}

open (OUT, ">$outputDir/read.dist.sample.$readLength.txt");
my $out="start.codon";
for (my $i=1; $i<=($end+$start+1); $i++) {
	if (! exists ($count5{$i})) {
		$count5{$i}=0;
	}
	$out.="\t".sprintf("%.3e", $count5{$i}/$count5{"0"}/$tot*1000000);
}
print OUT $out."\n";

my $out="stop.codon";
for (my $i=1; $i<=($end+$start+1); $i++) {
	if (! exists ($count3{$i})) {
		$count3{$i}=0;
	}
	$out.="\t".sprintf("%.3e", $count3{$i}/$count5{"0"}/$tot*1000000);
}
print OUT $out."\n";
close OUT;

open (RI, ">$outputDir/readDist.plot.$readLength.R");
print RI "A <- read.table (\"$outputDir/read.dist.sample.$readLength.txt\", sep=\"\\t\")"."\n";
print RI "loc1 <- c(-$start:$end)"."\n";
print RI "loc2 <- c(-$end:$start)"."\n";
print RI "B1 <- apply(A[1,2:ncol(A)], 2, sum)"."\n";
print RI "B2 <- apply(A[2,2:ncol(A)], 2, sum)"."\n";
print RI "pdf (file=\"$outputDir/plot.readDist.$readLength.pdf\")"."\n";
print RI "plot(loc1, B1, xlab=\"Around start codon\", ylab=\"RPM\", type=\"h\", lwd=2, ylim=c(0,max(c(B1,B2))))"."\n";
print RI "plot(loc2, B2, xlab=\"Around stop codon\", ylab=\"RPM\", type=\"h\", lwd=2, ylim=c(0,max(c(B1,B2))))"."\n";
print RI "dev.off()"."\n";
close RI;

my $comm="Rscript $outputDir/readDist.plot.$readLength.R";
system ($comm);

