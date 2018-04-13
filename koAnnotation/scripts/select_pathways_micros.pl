#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

use DBI;

## get pathways from a list for a specific species (enccu, enche, encin or nosce)
## v1.0 (2015.06.11)

sub usage {
    my $msg = shift;
    print "example: perl select_pathways_lca.pl -i pathways list -s species name\n";
    print "-i pathways list\n";
    print "-s species name (enccu, enche, encin, nosce)\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_s);
getopts('i:s:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No pathways list given\n");

open(IN,"$input") || die "Cannot open $input file !!\n";
my @in = <IN>;
close (IN);

my $spec = ($opt_s) ? $opt_s : usage("ERROR: No species name given\n");
#unless ($spec eq "enccu" or $spec eq "encin" or $spec eq "enche" or $spec eq "nosce"){
#	usage("ERROR: species name is not available!!");
#}

## output files
my $output = "$spec\_pathways.out";
open(OUT,">$output") || die "Cannot create $output!!\n";
open(SUM,">$output.summary") || die "Cannot create $output.summary!!\n";
print SUM "#PathwayID\t#PathwayName\t#Group\t#Category\t#Number of $spec proteins\t#Number of $spec KOs\t#Total KOs\t#Number of $spec reactions\t#Total reactions\n";
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

	### get all proteins for this pathway
	my %ko; my %prot;
	# query data
	my $sql = "select ko_map.map_id,ko_map.ko_id,$spec\_ko.id from ko_map,$spec\_ko where ko_map.map_id=$pathID and ko_map.ko_id=$spec\_ko.ko_id;";
	# print "$sql\n";
	my $sth = $dbh->prepare($sql);

	# excute the query
	$sth->execute();
	while (my @row = $sth->fetchrow_array()){
		# print "@row #####\n";<>;
		$ko{$row[1]} = 1;
		unless($prot{$row[2]}){
			$prot{$row[2]} = $row[1];
		} else {
			$prot{$row[2]} .= "\t".$row[1];
		}
	}
	$sth->finish();

	### print KOs_Pathways result
#	print "all proteins of this pathway $pathID:\n";
	print OUT "### $pathID\t$pathName\n";
	foreach my $protID(keys %prot){
#		print $protID,"\t",$prot{$protID},"\n";
		print OUT $protID,"\t",$prot{$protID},"\n";
	}
	print OUT "\n";

	### get all proteins together with their reactions for this pathway
	my %rn;
	my $rnSQL = "select distinct($spec\_ko.id),rn_map.rn_id from rn_map,rn_ko,$spec\_ko where rn_map.map_id=$pathID and rn_map.rn_id=rn_ko.rn_id and $spec\_ko.ko_id=rn_ko.ko_id;";
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

	# get number of this species' reactions take part in this pathway
	my $countlcaRnSQL = "select count(distinct rn_map.rn_id) from rn_map,rn_ko,$spec\_ko where rn_map.map_id=$pathID and rn_map.rn_id=rn_ko.rn_id and $spec\_ko.ko_id=rn_ko.ko_id";
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
#	print scalar(keys %prot),"\t";
#	print "Total KOs = $totalKO\n";
	print SUM $line,"\t",scalar(keys %prot),"\t",scalar(keys %ko),"\t",$totalKO,"\t",$countlcaRN,"\t",$totalRN,"\n";

	### stat
	$lcaKOType{$pathGroup} += scalar(keys %ko);
	$lcaRNType{$pathGroup} += $countlcaRN;
	$totalKOType{$pathGroup} += $totalKO;
	$totalRNType{$pathGroup} += $totalRN;
#	<>;
}

print SUM "\n#########################################\n";
foreach (@allType){
#	print $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";<>;
	print SUM $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";
}

close (OUT);
close (REAC);
close (SUM);

open(STAT,">$output.stat");
print STAT "type\tlcaProt\ttotalKO\tlcaRN\ttotalRN\n";
foreach (@allType){
#	print $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";<>;
	print STAT $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";
}
close (STAT);
# Disconnect from the database.
$dbh->disconnect();

exit;
