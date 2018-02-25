#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

## compare pathway statistics between LCA and extant microsporidia
## v1.0 (2015.12.18)

sub usage {
    my $msg = shift;
    print "example: perl compare_pathways.pl -i lca_pathways.out.summary -m micros_pathways -o pathways_compare.stat\n";
    print "-i statistic summary of LCA\n";
    print "-m folder contains extant micros pathways\n";
    print "-o output file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_o,$opt_m);
getopts('i:o:m:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No lca_pathways.out.SUMMARY given\n");
my $microsFol = ($opt_m) ? $opt_m : usage("ERROR: No folder micros_pathways given\n");
my $output = ($opt_o) ? $opt_o : usage("ERROR: No output file name given\n");#"lca_pathways.out";

### list of pathway types
my %lcaKOType; my %totalKOType; my %microsKOType;
my %lcaRNType; my %totalRNType; my %microsRNType;
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
my %allType = map { $_ => 1 } @allType;

### get micros KOs for each pathway
my %microsName = ("enche" => "E.hellem","encin" => "E.intestinalis","enccu" => "E.cuniculi","nosce" => "N.ceranae");
my @microsList = sort keys %microsName;
my %microsKO;
my %microsRN;
foreach my $micros(sort @microsList){
  open(MICROS,"$microsFol/$micros\_pathways.out.summary") || die "Cannot open $microsFol/$micros\_pathways.out.summary\n";
  foreach my $line(<MICROS>){
    if($line =~ /\d+/){
      chomp($line);
      my @tmp = split(/\t/,$line);  # #PathwayID	#PathwayName	#Group	#Category	#Number of enche proteins	#Number of enche KOs	#Total KOs	#Number of enche reactions	#Total reactions

      $microsKO{"$tmp[0]#$micros"} = $tmp[5];
      $microsRN{"$tmp[0]#$micros"} = $tmp[7];
    }
  }
  close (MICROS);
}

open(OUTKO,">$output.KO") || die "Cannot creat $output.KO!\n";
print OUTKO "#PathwayID\t#PathwayName\t#Group\t#Category\t#Total KOs\tLCA\tEnche\tEncin\tEnccu\tNosce\tminDiff\n";
open(OUTRN,">$output.RN") || die "Cannot creat $output.RN!\n";
print OUTRN "#PathwayID\t#PathwayName\t#Group\t#Category\t#Total RNs\tLCA\tEnche\tEncin\tEnccu\tNosce\tminDiff\n";
my %cat;  # $cat{group} = category
open(IN,"$input") || die "Cannot open $input file !!\n";
foreach my $line(<IN>){
  if($line =~ /\d+/){
    chomp($line);
    my @tmp = split(/\t/,$line);	# [0]PathwayID	[1]PathwayName	[2]Group	[3]Category	[4]LCA proteins	[5]LCA KOs	[6]Total KOs	[7]LCA reactions	[8]Total reactions
    print OUTKO "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[6]\t$tmp[5]";
    print OUTRN "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[8]\t$tmp[7]";
    my $maxKO=0;
    my $maxRN=0;
		$cat{$tmp[2]} = $tmp[3];

		### stat
		$lcaKOType{$tmp[2]} += $tmp[5];
		$lcaRNType{$tmp[2]} += $tmp[7];
		$totalKOType{$tmp[2]} += $tmp[6];
		$totalRNType{$tmp[2]} += $tmp[8];

    foreach(sort @microsList){
      print OUTKO "\t",$microsKO{"$tmp[0]#$_"};
      if($maxKO < $microsKO{"$tmp[0]#$_"}){
        $maxKO = $microsKO{"$tmp[0]#$_"};
      }
      print OUTRN "\t",$microsRN{"$tmp[0]#$_"};
      if($maxRN < $microsRN{"$tmp[0]#$_"}){
        $maxRN = $microsRN{"$tmp[0]#$_"};
      }
			### stat
			my $tmpID = $_."#".$tmp[2];
			$microsKOType{$tmpID} += $microsKO{"$tmp[0]#$_"};
			$microsRNType{$tmpID} += $microsRN{"$tmp[0]#$_"};
    }
    if($maxKO > 0){
      print OUTKO "\t",$tmp[5]/$maxKO,"\n";
    } else {
      print OUTKO "\t",0,"\n";
    }

    if($maxRN > 0){
      print OUTRN "\t",$tmp[7]/$maxRN,"\n";
    } else {
      print OUTRN "\t",0,"\n";
    }
  }
}
close (IN);

## print stat
open(STAT,">$output.stat");

### long format
print STAT "category\tpathway\tcount\tsource\ttype\n";
foreach my $type(@allType){
	print STAT $cat{$type},"\t",$type,"\t",$lcaKOType{$type},"\t","01_LCA_Microsporidia","\t","KO","\n";
	for(my $i=0; $i < scalar @microsList; $i++){
		print STAT $cat{$type},"\t",$type,"\t",$microsKOType{"$microsList[$i]#$type"},"\t0",$i+2,"_",$microsName{$microsList[$i]},"\t","KO","\n";
	}
	print STAT $cat{$type},"\t",$type,"\t",$lcaRNType{$type},"\t","01_LCA_Microsporidia","\t","RN","\n";
	for(my $i=0; $i < scalar @microsList; $i++){
		print STAT $cat{$type},"\t",$type,"\t",$microsRNType{"$microsList[$i]#$type"},"\t0",$i+2,"_",$microsName{$microsList[$i]},"\t","RN","\n";
	}
}
close (STAT);

### wide format
# print STAT "type\ttotalKO\tlcaKO";
# foreach(sort @microsList){
# 	print STAT "\t",$_,"KO";
# }
# print STAT "\ttotalRN\tlcaRN";
# foreach(sort @microsList){
# 	print STAT "\t",$_,"RN";
# }
# print STAT "\n";
#
# foreach (@allType){
# #	print $_,"\t",$lcaKOType{$_},"\t",$totalKOType{$_},"\t",$lcaRNType{$_},"\t",$totalRNType{$_},"\n";<>;
# 	print STAT $_,,"","\t",$totalKOType{$_},"\t",$lcaKOType{$_};
# 	foreach my $key(sort keys %microsKOType){
# 		my @idTMP = split(/#/,$key);
# 		if($idTMP[1] eq $_){
# 			print STAT "\t",$microsKOType{$key};
# 		}
# 	}
# 	print STAT "\t",$totalRNType{$_},"\t",$lcaRNType{$_};
# 	foreach my $key(sort keys %microsRNType){
# 		my @idTMP = split(/#/,$key);
# 		if($idTMP[1] eq $_){
# 			print STAT "\t",$microsRNType{$key};
# 		}
# 	}
# 	print STAT "\n";
# }
# close (STAT);

print "FINISHED!! Check output\n$output.KO\n$output.RN\n$output.stat\n";
exit;
