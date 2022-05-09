#!/usr/bin/perl

use File::Basename;

my $checkm_tsv = shift;
my $location = shift;
my $destination = shift;

system("mkdir -p $destination");

open(IN, $checkm_tsv);
while(<IN>) {
	next if m/^Bin Id/;

	my @d = split(/\t/);

	my $comp = $d[11];
	my $cont = $d[12];

	my $bin = $d[0];

	print "$bin, $comp, $cont\n";

	if ($comp>=80 && $cont<=10) {
		print "$bin, $comp, $cont\n";
		print "Moving $location/$bin.fa to $destination/$base\n";
		system("cp $location/$bin.fa $destination/$base");
	}	
}
close IN;
