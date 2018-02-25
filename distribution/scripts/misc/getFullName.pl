#!/usr/bin/env perl
use strict;
use warnings;

my $taxonInfo = "/Users/trvinh/work/R_projects/phyloprofile_app/latest/data/taxonomyMatrix.txt";
my $map = "ncbi_mapping.txt";

open(TAX,$taxonInfo);
my %tax;
foreach my $line(<TAX>){
	chomp($line);
	my @tmp = split(/\t/,$line); # abbrName        ncbiID  fullName   rank
	$tax{$tmp[1]} = $tmp[2];
}

open(MAP,$map);
open(OUT,">$map.new");
foreach my $line(<MAP>){
	chomp($line);
	my @tmp = split(/\t/,$line);
	unless($tax{$tmp[1]}){ print $tmp[1];<>;}
	print OUT $line,"\t",$tax{$tmp[1]},"\n";
}
close (MAP);
close (OUT);
exit;
