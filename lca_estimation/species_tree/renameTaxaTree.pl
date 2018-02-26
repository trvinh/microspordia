#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

sub usage {
    my $msg = shift;
    print "example: perl renameTaxaTree.pl -i input.tree -m mapping_file.txt -c column\n";
    print "-i\tInput domain file\n";
    print "-m\tName mapping file\n";
    print "-c\tColumn contains full names\n";
    die $msg."\n";
}

# global variables
our($opt_i, $opt_c, $opt_m);
getopts('i:c:m:');

# sanity checks
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $nameIn = ($opt_m) ? $opt_m : usage("ERROR: No name mapping file given\n");
my $col = ($opt_c) ? $opt_c : usage("ERROR: No column given\n");

# parse names
open(NAME,$nameIn) || die "Cannot open $nameIn!!\n";
my %name;
foreach my $line(<NAME>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# aspni	A.nidulans	Aspergillus nidulans

	$name{$tmp[0]} = $tmp[$col-1];
}
close (NAME);

### MAIN
open (INPUT, "$input") or die ("Cannot open $input file!!\n");
my @in = <INPUT>;
close (INPUT);
my $in = join("",@in);
my $new = $in;

foreach my $abbr (keys %name){
	if($new =~ /$abbr\:/){
		$new =~ s/$abbr\:/$name{$abbr}\:/g;
	}
}

# print output
open(OUT,">$input.renamed");
print OUT $new;
close (OUT);

print "Finished!! Check output $input.renamed\n";
exit;
