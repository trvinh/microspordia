#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

### check if there is any abbr gene names in the domain files
my $fol = "/Users/trvinh/work/OLD/distribution/mDomain_files";

my @files = glob("$fol/*.domains");
my %gene;
my $c = 1;
foreach my $file(@files){
	open(IN,$file);
	foreach my $line(<IN>){
		chomp($line);
		my @tmp = split(/\t/,$line);
		unless($tmp[1] =~ /@/){
			my @protTMP = split(/\:/,$tmp[1]);
			$gene{$protTMP[0]} = $file;
		}
	}
	close (IN);
	print $c,"/",scalar(@files),"\t",$file,"\n";
	$c++;
}
open(OUT,">errorID.list");
foreach(keys %gene){
	print OUT $_,"\t",$gene{$_},"\n";
}
close (OUT);
exit;
