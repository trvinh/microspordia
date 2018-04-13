#!/usr/bin/env perl
use Cwd;
use Cwd 'abs_path';
use Getopt::Std;

=description
run calcPfam.pl script for an input folder containing multiple lists
Version 1.0 (01.03.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl runCalcPfam.pl -i input_folder -t o\n";
    print "-i\tInput folder\n";
    print "-t\tType of proteins: (o)rphan or (h)omologous\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_t);
getopts('i:t:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input folder given\n");
my $type = ($opt_t) ? $opt_t : usage("ERROR: No protein type (homologous or orphan) given\n");

my $ext = "";
my $script = "/home/vinh/Desktop/data/project/iniAnalysis/scripts";
if($type eq "h" or $type eq "H"){
	$ext = "homologous";
	$script .= "/calcPfamHomolog.pl";
} elsif($type eq "o" or $type eq "O"){
	$ext = "orphan";
	$script .= "/calcPfamOrphan.pl";
} else{
	usage("ERROR: -t has to be o or h!\n");
}

### get all input files
my @file = glob(abs_path($input)."/*.$ext");

foreach my $file(@file){
	### get file name
	my @fileTMP = split(/\//,$file); # filename = $fileTMP[@fileTMP-1]
	print $fileTMP[@fileTMP-1],"\n";

	### run script
	my $cmd = "perl $script -i $file";
	print $cmd;<>;
#	system($cmd);
#	system("echo $cmd | qsub -V -S /bin/bash -cwd -j y -N $fileTMP[@fileTMP-1] -r y -q all.q");
}

exit;

