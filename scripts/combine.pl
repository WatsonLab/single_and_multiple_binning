#!/usr/bin/perl

use File::Basename;

my $o;
my @names;
my $tadcounter = 0;

my @files = @ARGV;

foreach $file (@files) {

	open(IN, $file);
	while(<IN>) {
		last;
	}

	$tadcounter++;

	my $base = basename $file, ".txt";
	$base = $base . ".bam";
	push(@names, $base);
	while(<IN>) {
		chomp();

		my($cn,$cl,$tad,$bam,$var) = split(/\t/);

		$o->{$cn}->{contigInfo}->{cl} = $cl;
		$o->{$cn}->{contigInfo}->{tad} += $tad;

		$o->{$cn}->{bamInfo}->{$base}->{bam} = $bam;
                $o->{$cn}->{bamInfo}->{$base}->{var} = $var;

	}
	close IN;
}

# print titles
print "contigName\tcontigLen\ttotalAvgDepth";
foreach $name (@names) {
	print "\t$name\t$name" . "-var";
}
print "\n";

my @sortedcn = sort mysortfunction keys %{$o};

foreach $cn (@sortedcn) {
	
	my $cr = $o->{$cn}->{contigInfo};
	my $br = $o->{$cn}->{bamInfo};


	print $cn, "\t", $cr->{cl}, "\t", $cr->{tad} / $tadcounter;

	foreach $name (@names) {
		print "\t", $br->{$name}->{bam}, "\t", $br->{$name}->{var};
	}

	print "\n";
}


sub mysortfunction {

	my($a1,$a2) = split(/_/, $a);
	my($b1,$b2) = split(/_/, $b);

	$a1 =~ s/^k//;
	$b1 =~ s/^k//;

	return $a1 <=> $b1 || $a2 <=> $b2;

}
