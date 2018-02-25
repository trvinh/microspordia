#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

# convert a wide format phyloprofile file into long format
# 11.07.2017

sub usage {
    my $msg = shift;
    print "example: perl wide2long.pl -i demo/test.main -f FAS -s traceability\n";
    print "-i\tWide-format phyloprofile input\n";
		print "-f\tscore 1\n";
		print "-s\tscore 2\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f,$opt_s);
getopts('i:f:s:');

my $fileIn = ($opt_i) ? $opt_i : usage("ERROR: No input phyloprofile file given\n");
my $stScore = ($opt_f) ? $opt_f : usage("ERROR: No name for first score given\n");
my $ndScore = ($opt_s) ? $opt_s : usage("ERROR: No name for second score given\n");

### MAIN
open(IN,$fileIn) || die "Cannot open $fileIn!\n";
my @in = <IN>;
close (IN);

my $header = shift(@in);
chomp(my @header = split(/\t/,$header));

open(OUT,">$fileIn.long") || die "Cannot create $fileIn.long!\n";
print OUT "geneID	ncbiID	orthoID	$stScore	$ndScore\n";

foreach my $line(@in){
	chomp($line);
	my @tmp = split(/\t/,$line);
	my $geneID = $tmp[0];

	for(my $i = 1; $i < scalar(@header); $i++){
		my @hitTMP = split(/#/,$tmp[$i]);
#		print "$geneID\t$header[$i]\t$hitTMP[0]\t$hitTMP[1]\t$hitTMP[2]\n";<>;
		print OUT "$geneID\t$header[$i]\t$hitTMP[0]\t$hitTMP[1]\t$hitTMP[2]\n";
	}
}
close (OUT);

print "Finished!! Check output $fileIn.long\n";
exit;
