#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
input is list of otholog groups
using raxml to calculate gene tree for each group
then calculate patristic distances for each proteins based on that tree
v1.0 (05.07.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl patristicDistCalc.pl -i orthologs.list -f fasta_folder -o out\n";
    print "-i\tlist of ortholog groups\n";
    print "-f\tfolder contains all fasta files\n";
    print "-o\toutput folder\n";
    die $msg."\n";
}



# global variables
our($opt_i, $opt_f, $opt_o);
getopts('i:f:o:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input orthologs file given\n");
open (INPUT, "$input") or die ("Cannot open $input file!!\n");
my @in = <INPUT>;
close (INPUT);

my $fastaFol = ($opt_f) ? $opt_f : usage("ERROR: No fasta folder given\n");
if($fastaFol =~ /\/$/){$fastaFol =~ s/\/$//;}

### prepare directories
my $dir = cwd();
if(!(-d ("$dir/tmp/"))){  # if folder doesnt exist create it
	mkdir("$dir/tmp");
}
my $output = ($opt_o) ? $opt_o : usage("ERROR: No output folder given\n");
if(!(-d ("$dir/$output/"))){  # if folder doesnt exist create it
	mkdir("$dir/$output");
}

##### MAIN
my $n = 1;
chdir("$dir/tmp");
foreach my $line(@in){
	chomp($line);
	my @tmp = split(/\t/,$line);
	# groupID
	my $groupID = shift(@tmp);

	# get fasta for proteins in this group and save them into $groupID.fasta
	open(FA,">$dir/tmp/$groupID.fasta") || die "Cannot create $dir/tmp/$groupID.fasta!!\n";
	my %tmpIDs;	# used to store original protein IDs
	my $index = 0;
	foreach my $prot (@tmp){
		my $fasta = getFasta($prot);
		my @fastaTMP = split(/\n/,$fasta);	# $fastaTMP[0] = ID, $fastaTMP[1] = Seq

		$tmpIDs{"seq$index"} = $fastaTMP[0];

		### replace U by X to prevent problem with raxml
		my $seq = $fastaTMP[1];
		$seq =~ s/U/X/g;
		print FA ">seq$index\n",$seq,"\n";
#		print FA ">seq$index\n",$fastaTMP[1],"\n";
		$index ++;
	}
	close (FA);
#	print "FINISH!\n";<>;

	### do tree reconstruction
	system("mafft-linsi --anysymbol $dir/tmp/$groupID.fasta > $dir/tmp/$groupID.fasta.aln");	# do alignment
	system("perl /home/vinh/Desktop/scripts/fasta/joinLines.pl -i $dir/tmp/$groupID.fasta.aln");	# join multiple lines into one line per sequence
	system("perl /home/vinh/Desktop/data/project/tree_scripts/degapper.pl -in=$dir/tmp/$groupID.fasta.aln.edited -limit=0.3");	# do de-gapping
	system("perl /home/vinh/Desktop/data/project/tree_scripts/MFAtoPHY.pl -i $dir/tmp/$groupID.fasta.aln.edited.proc");	# convert fa into phylip format
	system("raxmlHPC-PTHREADS -s $dir/tmp/$groupID.fasta.aln.edited.proc.phy -m PROTGAMMAILGF -n $groupID -p 10 -T 4");	# run raxml

	### calculate patristic distance
	system("python $dir/calculatePatristicDistances.py $dir/tmp/RAxML_result.$groupID $dir/tmp/$groupID.dist");

	open(OUT,">$dir/$output/$groupID.dist");
	open(my $fh,"$dir/tmp/$groupID.dist");
	while(my $row = <$fh>){
		chomp($row);
		my @rowTMP = split(/\t/,$row);
		print OUT $tmpIDs{$rowTMP[0]},"\t",$tmpIDs{$rowTMP[1]},"\t",$rowTMP[2],"\t",$rowTMP[3],"\n";
	}
	close (OUT);
	### remove all temporary files
	system("rm $dir/tmp/$groupID.*");
  system("rm $dir/tmp/RAxML*.$groupID");

	print "$groupID finished...\n";
#	<>;
}

exit;


### get sequence from multifasta file for a protein (species_name:proteinID)
sub getFasta{
	my $desc = $_[0];	# $desc = SPECIES:PROTEIN_ID
	#print $desc;<>;
	# get species name
	my @tmp = split(/\:/,$desc);
	my $species = $tmp[0];  ## $tmp[0] = species name, $tmp[1] = protein_id
	# open multifasta file for this species
	my $file = "$fastaFol/$species.fasta";
	unless(-e $file){
		$file = "$fastaFol/$species.fa";
	}
	open(FILE, $file) or die "Cannot open the multiple fasta file of this $species ($file)!!!\n";
	my @file = <FILE>;
	# split into many fasta
	$file = "\n";
	$file .= join("",@file);
	my @genes = split("\n>",$file);
	# search for proteinID
	my $result = "";
	$desc =~ s/\./\\./g; # to find exact match in case the desc have "." in between, e.g. "cal:CaO19.5501"
	foreach my $gene (@genes){
		if(length($gene)>0 and ($gene =~ m/$desc\|/ or $gene =~ m/$desc\n/)){
			$result = $gene;
			last;
		}
	}
#	print "FASTA=",$result;<>;
	close (FILE);
	return $result;
}
