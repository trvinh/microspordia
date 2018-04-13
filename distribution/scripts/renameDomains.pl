#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

sub usage {
    my $msg = shift;
    print "example: perl renameDomains.pl -i domain_files -m mappingID.txt -s renameTaxa.pl\n";
    print "-i\tInput domain folder\n";
		print "-m\tInput mapping file\n";
    print "-s\tScript renameTaxa.pl\n";
    die $msg."\n";
}

# global variables
our($opt_i, $opt_s, $opt_m);
getopts('i:s:m:');

# sanity checks
my $folder = ($opt_i) ? $opt_i : usage("ERROR: No input folder given\n");
my $mapIn = ($opt_m) ? $opt_m : usage("ERROR: No mapping file given\n");
my $script = ($opt_s) ? $opt_s : usage("ERROR: No script file given\n");

# my $folder = "/Users/trvinh/work/OLD/distribution/mDomain_files";
my @files = glob("$folder/*.domains");

foreach my $file(@files){

	my $cmd = "perl $script -m $mapIn -i $file -o $file.renamed";
	system($cmd);
	system("mv $file.renamed $file");
	print $file,"\n";
}
exit;
