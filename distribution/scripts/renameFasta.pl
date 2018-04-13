#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

sub usage {
    my $msg = shift;
    print "example: perl renameTaxa.pl -i inputFolder -m names.tab\n";
    print "-i\tInput folder\n";
    print "-m\tName mapping file\n";
    die $msg."\n";
}

# global variables
our($opt_i, $opt_m);
getopts('i:m:');

# sanity checks
my $folder = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $nameIn = ($opt_m) ? $opt_m : usage("ERROR: No name mapping file given\n");



# parse names
open(NAME,$nameIn) || die "Cannot open $nameIn!!\n";
my %name;
foreach my $line(<NAME>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# astph   Astpho
	$name{$tmp[0]} = $tmp[1];
}
close (NAME);

### MAIN
my @files = glob("$folder/*.fa");
mkdir("$folder/new");
foreach my $file(@files){
	open (INPUT, "$file") or die ("Cannot open $file file!!\n");
	my @in = <INPUT>;
	close (INPUT);
	my $in = join("",@in);
	my $new = $in;

	my $spec = "";
	foreach my $abbr (keys %name){
		if($new =~ /$abbr\:/){
			$new =~ s/$abbr\:/$name{$abbr}\@/g;
			$spec = $name{$abbr};
		}
	}

	# print output
	open(OUT,">$file.renamed");
	print OUT $new;
	close (OUT);

	# print "Finished!! Check output $output\n";

	system("mv $file.renamed $folder/new/$spec.fa");
	print $file,"\n";
}
exit;
