#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

sub usage {
    my $msg = shift;
    print "example: perl long2fasta.pl -i input.long -f fasta-folder\n";
    print "-i\tinput.main.long file\n";
		print "-f\tFasta folder\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f);
getopts('i:f:');

my $input = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $FasFol = ($opt_f) ? $opt_f : usage("ERROR: No fasta folder given\n");

open(IN,$input) || die "Cannot open $input!\n";
foreach my $line(<IN>){
		chomp($line);
		my @tmp = split(/\t/,$line);
		if($tmp[2] ne "NA" && $tmp[2] =~ /@/){
			# print"@tmp\n";
			my @specTMP = split(/@/,$tmp[2]);
			open(FAS,"$FasFol/$specTMP[0]\@$specTMP[1]\.fa") || die "Cannot open $FasFol/$specTMP[0]\@$specTMP[1]\.fa\n";
			my @fas = <FAS>;
			close (FAS);
			my $fas = join("",@fas);
			my @allSeq = split(/>/,$fas);
			foreach my $seq(@allSeq){
				chomp($seq);
				if($seq =~ /$tmp[2]\n/){
					# print $seq;<>;
					my @seqTMP = split(/\n/,$seq);
					my $header = join('|',@tmp);
					print ">$header\n$seqTMP[1]\n";
				}
			}
		}
}
close (IN);

exit;
