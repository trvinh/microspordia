#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
create a matrix for phyloprofile tool
from a list of ortholog groups
12.06.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl phyloProfile.pl -i ortholog.list.NEW -t taxa.list\n";
    print "-i\tOutholog groups list\n";
    print "-t\tFile contain ncbi taxonomy info of all taxa under the study\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_t);
getopts('i:t:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input ortholog list given\n");
my $taxaFile = ($opt_t) ? $opt_t : usage("ERROR: No file contain list of all taxa given\n");

### output
open(OUT,">$input.phyloprofile");
print OUT "geneID";

#### get super taxa name
open(TAXA,$taxaFile) || die "Cannot open $taxaFile!!\n";
my @taxa = <TAXA>;
close (TAXA);
my %ncbiID;	# $ncbiID{abbrName} = "ncbi"ncbiID (e.g. ncbi1355, ncbi8943,...)

foreach my $line(@taxa){
	chomp($line);
	if($line =~ /\d/ && length($line)>1){
		my @lineTMP = split(/\t/,$line);	# abbrName        ncbiID
		$ncbiID{$lineTMP[0]} = "ncbi".$lineTMP[1];
	}
}

## print list of all taxa into output file (header)
foreach my $abbrName(sort keys %ncbiID){
	print OUT "\t$ncbiID{$abbrName}";
}
print OUT "\n";

### MAIN ###
### read ortholog input file, get their max FAS for each species in @allSpecies list and write to a matrix output
open(IN,$input) || die "Cannot open $input!!\n";
my @in = <IN>;
close (IN);

my $c = 1;
foreach my $group(sort @in){
	chomp($group);
	my @prots = split(/\t/,$group);
	my $groupID = shift(@prots);
	print OUT "$groupID";

	### check for the present of every species in @allSpeices list
	foreach my $specID(sort keys %ncbiID){
#		print $specID,"\n";
		## if this group has ortholog in this species, than get FAS score
		if($group =~ /$specID(_)*(.)*?\:(.)+?\t/){
			my $orthoID = $&;	$orthoID=~ s/\t$//;
#			print $id," - ",$maxID{$id}."\$".$maxFas{$id},"\n";<>;
			print OUT "\t$orthoID#1#NA";
		}
		## else, return NA
		else {
#			print "NA\n";
			print OUT "\tNA";
		}
	}

	print OUT "\n";
	print $c,"/",scalar(@in),"\n"; $c++;
}

close (OUT);

print "FINISHED!! Check output in $input.phyloprofile\n";
exit;

