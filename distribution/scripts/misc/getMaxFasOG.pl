#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
calculate max/mean FAS for lca OG group
=cut

open(IN,"/Users/trvinh/work/OLD/distribution/lca.list.distribution.FAS");

my %stat;
my %maxSpec;
my $c = 1;
foreach my $line(<IN>){
	chomp($line);
	my @tmp = split(/\t/,$line);
	unless($stat{$tmp[0]}){
		$stat{$tmp[0]} = $tmp[@tmp-1];
		$maxSpec{$tmp[0]} = $tmp[1];
	} else {
		if($stat{$tmp[0]} < $tmp[@tmp-1]){
			$stat{$tmp[0]} = $tmp[@tmp-1];
			$maxSpec{$tmp[0]} = $tmp[1];
		}
	}
	print "$c\n";$c++;
}
close (IN);

open(OUT,">maxFAS.out");
foreach(sort keys %stat){
	print OUT "$_\t$stat{$_}\t$maxSpec{$_}\n";
}
close (OUT);

print "FINISHED! Check output maxFAS.out\n";

exit;
