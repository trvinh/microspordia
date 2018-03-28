#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
parse GO annotation from blast2go for list of ortholog groups
28.03.2018
=cut

sub usage {
    my $msg = shift;
    print "example: perl parseGO.pl -i blast2go.output -l microsOnly.list\n";
    print "-i\tOutput from blast2go\n";
    print "-l\tList of ortholog groups\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_l);
getopts('i:l:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No blast2go file given\n");
my $orthoList = ($opt_l) ? $opt_l : usage("ERROR: No list of orthogroups given\n");

### get GO annotations
open(GO,$input) || die "Cannot open $input!\n";
my %go; 	# $go{protID} = "goID\tDesc"
foreach my $line(<GO>){
	unless($line =~ /^(Tags)/){
		chomp($line);
		my @tmp = split(/\t/,$line);
		unless($tmp[@tmp-1] eq "no GO terms" or $tmp[@tmp-1] eq "no IPS match" ){
			# print $tmp[2],"\t",$tmp[@tmp-2],"\t",$tmp[@tmp-1];<>;
			$go{$tmp[2]} = $tmp[@tmp-2]."\t".$tmp[@tmp-1];
		}
	}
}
close (GO);

### add GO info to ortholog groups
open(LIST,$orthoList) || die "Cannot open $orthoList!\n";
open(OUT, ">$orthoList.GO");
foreach my $line(<LIST>){
	chomp($line);
	my @tmp = split(/\t/,$line);
	my $groupID = shift(@tmp);
	my $out = "";
	foreach my $prot(@tmp){
		$prot =~ s/:/_/g;
		if($go{$prot}){
			$out .= $go{$prot}."\t";
		}
	}
	$out =~ s/\t$//;
	if(length $out > 0){
		print OUT $groupID,"\t",$out,"\n";
	}
}
close (LIST);
close (OUT);
print "FINISHED!! Check output in $orthoList.GO\n";
exit;
