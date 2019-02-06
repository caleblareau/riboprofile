if ($#ARGV < 2) {
	print "## remove adapter sequence"."\n";
	print "usage: remove.adaptor.pl fastqFile adapterSequence outputFile [readLengthCutoff]"."\n";
	print "fastqFile: raw fastq read sequences;"."\n";
	print "adapterSequence: 5' end sequence of adapter, 10nt is recommended;"."\n";
	print "outputFile: output file;"."\n";
	print "readLengthCutoff [optional]: read length cutoff after trimming adapters, default is 15 nt."."\n";
	exit;
}

my $file=$ARGV[0];
my $adapter=$ARGV[1];
my $outputFile=$ARGV[2];
my $cutoff;
if ($ARGV[3]) {
	$cutoff=$ARGV[3];
} else {
	$cutoff=15;
}

open (IN, "$file");
open (OUT, ">$outputFile");

my $num=0;
my $out;
my $m;
while (<IN>) {
	chomp;
	$num++;
	if ($num%4==1) {
		$out=$_."\n";
	} elsif ($num%4==2) {
		if ($_ =~ /$adapter/) {
			my @s=split /$adapter/, $_;
			$m=$s[0];
		} else {
			$m=$_;
		}
		$out.=$m."\n";
	} elsif ($num%4==3) {
		$out.=$_."\n";
	} elsif ($num%4==0) {
		if (length($m) >= $cutoff) {
			print OUT $out.substr($_, 0, length($m))."\n";
		}
	} 
}
close IN;
close OUT;
