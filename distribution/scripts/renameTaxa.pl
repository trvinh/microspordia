#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

sub usage {
    my $msg = shift;
    print "example: perl renameTaxa.pl -i input.domains -m ncbi_mapping -o output\n";
    print "-i\tInput domain file\n";
    print "-m\tName mapping file\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
our($opt_i, $opt_o, $opt_m);
getopts('i:o:m:');

# sanity checks
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $nameIn = ($opt_m) ? $opt_m : usage("ERROR: No name mapping file given\n");
my $output = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");

# parse names
open(NAME,$nameIn) || die "Cannot open $nameIn!!\n";
my %name;
foreach my $line(<NAME>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# thetr_4772_1	529818	Thecamonas trahens

	my $name = $tmp[2];
	my @nameTMP = split(/\s/,$tmp[2]);
	if(scalar @nameTMP > 1){
		if(length($nameTMP[1]) > 3){
			$name = substr($nameTMP[0],0,1).".".$nameTMP[1];
		}
	}

	$name{$tmp[0]} = $name."@".$tmp[1];
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
		$new =~ s/$abbr\:/$name{$abbr}\@/g;
	}
}

# print output
open(OUT,">$output");
print OUT $new;
close (OUT);

print "Finished!! Check output $output\n";
exit;
