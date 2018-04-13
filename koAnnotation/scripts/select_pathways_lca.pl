#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

use DBI;

## get pathways from a list for LCA proteins
## v1.0 (2015.05.29)

sub usage {
    my $msg = shift;
    print "example: perl select_pathways_lca.pl -i pathways list -o output\n";
    print "-i pathways list. e.g. pathways.list\n";
    print "-o output file. e.g. lca_pathways.out\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_o);
getopts('i:o:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No pathways list given\n");

open(IN,"$input") || die "Cannot open $input file !!\n";
my @in = <IN>;
close (IN);

## output files
my $output = ($opt_o) ? $opt_o : usage("ERROR: No output file name given\n");#"lca_pathways.out";
open(OUT,">$output") || die "Cannot create $output!!\n";
open(SUM,">$output.summary") || die "Cannot create $output.summary!!\n";
print SUM "#PathwayID\t#PathwayName\t#Group\t#Category\t#Number of LCA proteins\t#Number of LCA KOs\t#Total KOs\t#Number of LCA reactions\t#Total reactions\n";
open(REAC,">$output.rn");

# MySQL database configurations
my $database = "pathways_kegg";
my $host = "172.17.100.150";
my $username = "vinh";
#print "Enter password: ";
#system('/bin/stty', '-echo');
#my $password = <STDIN>;
my $password = "bemun2018";
#chomp($password);
print "\n";
system('/bin/stty', 'echo');

# connect to MySQL database
my $dbh = DBI->connect("DBI:mysql:database=$database;host=$host",$username, $password, {'RaiseError' => 1});

### list using for stat
my %lcaKOType;my %totalKOType;my %lcaRNType;my %totalRNType;
my @allType = ("Amino acid metabolism", "Metabolism of other amino acids",
		"Cell growth and death",
		"Carbohydrate metabolism",
		"Cell motility",
		"Energy metabolism",
		"Folding, sorting and degradation",
		"Lipid metabolism",
		"Metabolism of cofactors and vitamins",
		"Membrane transport",
		"Metabolism of terpenoids and polyketides",
		"Nucleotide metabolism",
		"Replication and repair",
		"Signal transduction",
		"Transport and catabolism",
		"Transcription",
		"Translation");

# read list of pathways
foreach my $line (@in){
	chomp($line);

	### get info of this pathway (e.g. ID = 00051	Name = Fructose and mannose metabolism	Group = Carbohydrate metabolism		Category = Metabolism)
	my @tmp = split(/\t/,$line);
	my $pathID = $tmp[0];
	my $pathName = $tmp[1];
	my $pathGroup = $tmp[2];
	my $pathCategory = $tmp[3];
	print $pathID,"\t",$pathGroup,"\n";

	### get all LCA proteins for this pathway
	my %ko; my %lcaProt;
	# query data
	my $sql = "select ko_map.map_id,ko_map.ko_id,lca_ko.lca_id,lca_ko.dist,lca_ko.fas from ko_map,lca_ko where ko_map.map_id=$pathID and ko_map.ko_id=lca_ko.ko_id;";
#	print "$sql\n";
	my $sth = $dbh->prepare($sql);

	# excute the query
	$sth->execute();
	while (my @row = $sth->fetchrow_array()){
#		print "@row #####\n";<>;
		$ko{$row[1]} = 1;
		unless($lcaProt{$row[2]}){
			$lcaProt{$row[2]} = $row[1]."#".$row[4]."#".$row[3];
		} else {
			$lcaProt{$row[2]} .= "\t".$row[1]."#".$row[4]."#".$row[3];
		}
	}
	$sth->finish();

	### print LCA KOs_Pathways result
#	print "all LCA of this pathway $pathID:\n";
	print OUT "### $pathID\t$pathName\n";
	foreach my $lcaID(keys %lcaProt){
#		print $lcaID,"\t",$lcaProt{$lcaID},"\n";
		print OUT $lcaID,"\t",$lcaProt{$lcaID},"\n";
	}
	print OUT "\n";

	### get all LCA together with their reactions for this pathway
	my %rn;
	my $rnSQL = "select distinct(lca_ko.lca_id),rn_map.rn_id from rn_map,rn_ko,lca_ko where rn_map.map_id=$pathID and rn_map.rn_id=rn_ko.rn_id and lca_ko.ko_id=rn_ko.ko_id;";
	$sth = $dbh->prepare($rnSQL);
	$sth->execute();
	while (my @row = $sth->fetchrow_array()){
#		print "@row #####\n";
		unless($rn{$row[0]}){
			$rn{$row[0]} = $row[1];
		} else {
			$rn{$row[0]} .= "\t".$row[1];
		}
	}
	$sth->finish();

	# print LCA Reaction_Pathway result
	print REAC "### $pathID\t$pathName\n";
	foreach my $rnID (sort keys %rn){
		print REAC $rnID,"\t",$rn{$rnID},"\n";
	}
	print REAC "\n";

	# get total number of KOs in this pathway
	my $countAllKoSQL = "select count(distinct ko_map.ko_id) from ko_map where ko_map.map_id=$pathID";
	$sth = $dbh->prepare($countAllKoSQL);
	$sth->execute();
	my $totalKO = "";
	while (my @row = $sth->fetchrow_array()){
#		print "@row #####\n";
		$totalKO = $row[0];
	}
	$sth->finish();

	# get number of LCA reactions take part in this pathway
	my $countlcaRnSQL = "select count(distinct rn_map.rn_id) from rn_map,rn_ko,lca_ko where rn_map.map_id=$pathID and rn_map.rn_id=rn_ko.rn_id and lca_ko.ko_id=rn_ko.ko_id";
	$sth = $dbh->prepare($countlcaRnSQL);
	$sth->execute();
	my $countlcaRN = "";
	while (my @row = $sth->fetchrow_array()){
#		print "@row #####\n";
		$countlcaRN = $row[0];
	}
	$sth->finish();

	# get total number of Reactions in this pathway
	my $countAllRnSQL = "select count(rn_map.rn_id) from rn_map where rn_map.map_id=$pathID;";
	$sth = $dbh->prepare($countAllRnSQL);
	$sth->execute();
	my $totalRN = "";
	while (my @row = $sth->fetchrow_array()){
#		print "@row #####\n";
		$totalRN = $row[0];
	}
	$sth->finish();


#	print "number of KOs:\t";
#	print scalar(keys %ko),"\t";
#	print "number of LCAs:\t";
#	print scalar(keys %lcaProt),"\t";
#	print "Total KOs = $totalKO\n";
	print SUM $line,"\t",scalar(keys %lcaProt),"\t",scalar(keys %ko),"\t",$totalKO,"\t",$countlcaRN,"\t",$totalRN,"\n";

	### stat
	$lcaKOType{$pathGroup} += scalar(keys %ko);
	$lcaRNType{$pathGroup} += $countlcaRN;
	$totalKOType{$pathGroup} += $totalKO;
	$totalRNType{$pathGroup} += $totalRN;

#	<>;
}
close (OUT);
close (REAC);
close (SUM);

open(STAT,">$output.stat");
print STAT "type\tlcaKO\ttotalKO\tlcaRN\ttotalRN\n";
foreach (@allType){
#	print $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";<>;
	print STAT $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";
}
close (STAT);
# Disconnect from the database.
$dbh->disconnect();

exit;
