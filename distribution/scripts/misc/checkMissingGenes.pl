#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

my $matrix = "/Users/trvinh/work/OLD/distribution/lca.list.distribution.FASmatrix";
my $fol = "/Users/trvinh/work/OLD/distribution/mDomain_files";
my $pair = "/Users/trvinh/work/OLD/distribution/calcFASpairs.list";

open(MATRIX,$matrix);
open(OUT,">missingGene.list");
open(PAIR,">calcFASagain.list");
my %spec;
foreach my $line(<MATRIX>){
	chomp($line);
	if($line =~ /^OG/){
		my @tmp = split(/\t/,$line);
		my $group = shift(@tmp);

		if(-e "$fol/$group.domains"){
			open(IN,"$fol/$group.domains");
			my @domains = <IN>;
			close (IN);
			my $domains = join("",@domains);

			foreach my $item(@tmp){
				my @itemTMP = split(/#/,$item); 	# D.yakuba@7245@1113#0.49426229299
				unless($itemTMP[0] eq "NA"){
					unless($domains =~ /\t\Q$itemTMP[0]\E\t/){
						if($itemTMP[1] > 0){
							print OUT "$itemTMP[0]\t$group\n";
							my $pairs = `less $pair | grep $group | grep $itemTMP[0]`;
							print PAIR $pairs;
						}
					}
				}
			}
		} else {
			# print OUT "no domain file for $group\n";
		}

		print $group,"\n";
	}
}
close (MATRIX);
close (OUT);
close (PAIR);
exit;
