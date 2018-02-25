#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;
use POSIX;

# get a list of LCA groups together with their KOs (lca_ko.tab)
# used for importing into database
# v2.0 (16.06.2015)
# include KO annotation level (1 for microsporidia, 2-5 for fungi, 6 for metazoa & co,...)
# and color code for KEGG Mapping
# v3.0 (30.05.2016)
# include FAS scores
# v4.0 (16.08.2016)
# replace FAS scores by confident values
# confValue = FAS * exp(-alpha*dist)

sub usage {
    my $msg = shift;
    print "example: perl getKoTab.pl -i orthologs.list.KO.SUMMARY\n";
    print "-i\tinput list of orthologs groups and their KOs\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");

open(IN,$input) || die "Cannot open $input!!\n";
my @input = <IN>;
close (IN);

# OUTPUT
open(OUT,">$input.ko_tab");

# MAIN
my %color_code = (
	0 => '#AE250A',
	1 => '#A8390B',
	2 => '#A24E0C',
	3 => '#9C630D',
	4 => '#97780E',
	5 => '#918D0F',
	6 => '#8BA110',
	7 => '#86B611',
	8 => '#80CB12',
	9 => '#7AE013',
	10 => '#75F514',
);
my $i = 10001;
foreach my $line(@input) {
	if(length($line)>2){
		chomp ($line);	# groupID	koID#FAS#distance

		my @tmp = split(/\t/,$line);	# OG_1263	K11964#0.32909988351#0.619982244352527	K13513#0.99206150709#0.427564784152742
		my $groupID = shift (@tmp);
		foreach my $ko(@tmp){
			my @ko_tmp = split(/#/,$ko);
			my $koID = $ko_tmp[0];
			my $dist = $ko_tmp[2];
			my $fas = $ko_tmp[1];

#			my $confLv = $level;
#			my $confLv = ceil($level*0.5 + (11-ceil($fas*10))*0.5);
#			my $confVl = $fas*exp(-0.25*$dist);
#			my $colorLv = 0;
#			if($confVl ne 0){
#				$colorLv = int($confVl*10 + $confVl*10/abs($confVl*10*2));
#			}

      my $colorLv = "";
      if($dist eq "NA"){
        $colorLv = 0;
      } else {
        $colorLv = 10-int($dist*10);
      }
			my $fgcolor = "#000000";
			if($colorLv < 7) {$fgcolor = "#FFFFFA";}

#			print OUT $i,"\t",$groupID,"\t",$koID,"\t",$confVl,"\t",$color_code{$colorLv},",",$fgcolor,"\t",$fas,"\n";
      # unless($color_code{$colorLv}){print "HERE $groupID - $colorLv";<>;}
			print OUT $i,"\t",$groupID,"\t",$koID,"\t",$dist,"\t",$color_code{$colorLv},",",$fgcolor,"\t",$fas,"\n";
			$i++;
		}
	}
}

print "Check output file $input.ko_tab!\n";
close (OUT);
exit;
