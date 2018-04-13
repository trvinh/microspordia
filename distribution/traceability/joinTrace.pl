#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

sub usage {
    my $msg = shift;
    print "example: perl joinTrace.pl -i lca.list.distribution.phyloprofile -t vinhOgs_traceprofile.mod.txt -o output\n";
    print "-i\tInput phyloprofile file\n";
    print "-t\Input traceability matrix file\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
our($opt_i, $opt_t, $opt_o);
getopts('i:o:t:');

# sanity checks
my $file = ($opt_i) ? $opt_i : usage("ERROR: No phyloprofile file given\n");
my $arpit = ($opt_t) ? $opt_t : usage("ERROR: No traceability file given\n");
my $output = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");


# my $arpit = "vinhOgs_traceprofile.mod.txt";
open(ARPIT,$arpit);
my @arpit = <ARPIT>;
close (ARPIT);
my @arpitHeader = split(/\t/,$arpit[0]);

my %trace;
for(my $i = 1; $i < scalar(@arpit); $i++){
	chomp(my $line = $arpit[$i]);
	my @tmp = split(/\t/,$line);

	for(my $n=1; $n < scalar @arpitHeader; $n++){
		$trace{"$tmp[0]#$arpitHeader[$n]"} = $tmp[$n];
#		print "$tmp[0]#$arpitHeader[$n]#$tmp[$n]";<>;
	}
}
#print(rand());<>;

my $micros = "ncbi31281;ncbi40302;ncbi993615;ncbi103449;ncbi27973;ncbi6035;ncbi876142;ncbi1288291;ncbi70536;ncbi278021;ncbi586133;";
my $selTaxa = $arpit[0].$micros;
$selTaxa =~ s/\s+/;/g;

my $file = "lca.list.distribution.FASmatrix_abbr";
open(IN,$file);
my @in = <IN>;
close (IN);
my $header = $in[0];
my @header = split(/\t/,$header);

open(OUT,$output);
for(my $i=0; $i<scalar @in; $i++){
	chomp(my $line = $in[$i]);
	my @tmp = split(/\t/,$line);

	my $newLine = "";
	for(my $n=0; $n < scalar(@header); $n++){
		if($selTaxa =~ /$header[$n];/){
			if($i == 0){
				#print $header[$n];
				$newLine .= $header[$n]."\t";
			} else {
				unless($tmp[$n] =~ /^OG_/){
					my $info = $tmp[$n];
					if($info eq "NA"){
						$info = "NA#NA";
					}
					if($trace{"$tmp[0]#$header[$n]"}){
#						print "$tmp[$n]#",$trace{"$tmp[0]#$header[$n]"};<>;
						$newLine .= $info."#".$trace{"$tmp[0]#$header[$n]"}."\t";
					} else {
						#print "$tmp[$n]#",rand();<>;
						$newLine .= $info."#".rand()."\t";
					}
				} else {
					#print $tmp[$n];
					$newLine .= $tmp[$n]."\t";
				}
			}
		}
	}
	$newLine =~ s/\t$//;
	print OUT $newLine,"\n";
}
close (OUT);

exit;
