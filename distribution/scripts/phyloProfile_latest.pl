#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
create matrix of FAS scores for all ortholog groups.
23.05.2016
=cut

sub usage {
    my $msg = shift;
    print "example: perl phyloProfile.pl -i lca.list.distribution -m lca.list.distribution.mFAS -t ncbiID_mapping.txt\n";
    print "-i\tOutholog groups list\n";
    print "-m\tFile contain max/mean FASs\n";
    print "-t\tFile contain ncbi taxonomy ID of all taxa under the study\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_m,$opt_t);
getopts('i:m:t:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input ortholog list given\n");
my $taxaFile = ($opt_t) ? $opt_t : usage("ERROR: No file contain list of all taxa given\n");
my $fasFile = ($opt_m) ? $opt_m : usage("ERROR: No FAS file given\n");

### get FACT scores
open(FAS,$fasFile) || die "Cannot open $fasFile!!\n";
my @fas = <FAS>;
close (FAS);
my %fas;	# $fas{groupID#protID} = meax_FAS
foreach my $fasLine (@fas){
	chomp ($fasLine);
	my @fasLineTMP = split(/\t/,$fasLine);	# OG_1001 ciosa_2525_21:27575     0.71355995648
	my $id = $fasLineTMP[0]."#".$fasLineTMP[1];
	$fas{$id} = $fasLineTMP[2];
#	print $id," - ",$fasLineTMP[2];<>;
}

### output
open(OUT,">$input.phyloprofile");
print OUT "geneID";

#### get super taxa name
open(TAXA,$taxaFile) || die "Cannot open $taxaFile!!\n";
my @taxa = <TAXA>;
close (TAXA);
my %ncbiID;	# $ncbiID{abbrName} = "ncbi"ncbiID (e.g. ncbi1355, ncbi8943,...)
# my %ncbiName; #$ncbiID{abbrName} = fullName

foreach my $line(@taxa){
	chomp($line);
	if(length($line)>1){
		my @lineTMP = split(/\t/,$line);	# taxaName    ncbiID	fullName
		$ncbiID{$lineTMP[0]} = "ncbi".$lineTMP[1];
		# my @nameTMP = split(/\s/,$lineTMP[2]);
		# my $name = $lineTMP[2];
		# if(scalar @nameTMP > 1){
		# 	if(length($nameTMP[1]) > 3){
		# 		$name = substr($nameTMP[0],0,1).".".$nameTMP[1];
		# 	}
		# }
		# $ncbiName{$lineTMP[0]} = $name."@".$lineTMP[1];
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

	### get max FAS score for each species in this group
	my %maxFas;	# $maxFas{SPECIES} = maxFAS
	my %maxID;	# ID of protein that has maxFAS
	foreach my $ortho(@prots){
		my @orthoTMP = split(/:/,$ortho);
		my $id = $groupID."#".$orthoTMP[0];

		if($fas{"$groupID#$ortho"} >= 0){
			unless($maxFas{$id}){
#				print "$id - ",$fas{"$groupID#$ortho"};<>;
				$maxFas{$id} = $fas{"$groupID#$ortho"};
				$maxID{$id} = $ortho;
			} else {
				if($maxFas{$id} < $fas{"$groupID#$ortho"}){
					$maxFas{$id} = $fas{"$groupID#$ortho"};
					$maxID{$id} = $ortho;
				}
			}
		} else {
			print "NO mFAS for $id!!\n";<>;
		}
	}

	### check for the present of every species in @allSpeices list
	foreach my $specID(sort keys %ncbiID){
#		print $specID,"\n";
		## if this group has ortholog in this species, than get FAS score
		if($group =~ /$specID(_)*(.)*?\:/){
			my $specID = $&; $specID =~ s/\:$//;
			my $id = $groupID."#".$specID;	#$id=~ s/\:$//;
			# print $id," - ",$maxID{$id}."\t".$maxFas{$id},"\n";<>;
			my $maxGene = $maxID{$id}; #$maxGene =~ s/$specID:/$ncbiName{$specID}@/;
			print OUT "\t$maxGene\#$maxFas{$id}";
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
