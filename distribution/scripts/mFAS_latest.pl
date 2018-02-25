#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;
use List::Util qw(sum max min);
use Scalar::Util qw(looks_like_number);

=desc
calculate max (or mean [optional]) FAS score for each ortholog in its group
v1.0 (03.03.2016)
and get domain architectures of corresponding protein pairs that have MAX score (mean score will not work because it is impossible to get "mean" architecture!!)
(28.10.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl phyloScore_FAS.pl -i orthologGroups.list -f orthologGroups.list.FAS -d domain_folder -m max\n";
    print "-i\tOutholog groups list\n";
    print "-f\tFAS file\n";
    print "-d\tDomains folder\n";
    print "-m\tmax / mean\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f,$opt_d,$opt_m);
getopts('i:f:d:m:');

# sanity checks;
my $orthoFile = ($opt_i) ? $opt_i : usage("ERROR: No input ortholog list given\n");
my $fasFile = ($opt_f) ? $opt_f : usage("ERROR: No FAS file given\n");
my $domainFol = ($opt_d) ? $opt_d : usage("ERROR: No domain folder given\n");
my $calc = ($opt_m) ? $opt_m : usage("ERROR: You have to set -m as max or mean\n");

#unless($calc eq "max" or $calc eq "mean"){
#	usage("-m has to be max or mean!!\n");
#}
$calc = "max";

### open and parse FAS output file
open(FAS,$fasFile) || die "Cannot open $fasFile!!\n";
my @fasIN = <FAS>;
close (FAS);

my %fas;	# $fas{groupID#protID} = FAS1;FAS2;FAS3;...
foreach my $fasLine(@fasIN){
	chomp($fasLine);
#	print $fasLine;<>;
	my @tmp = split(/\t/,$fasLine);	# OG_1501 acama_4692_1:4376       encin_5517_1:Eint_081890        0.93730814314
	my $id = $tmp[0]."#".$tmp[1];	# $id = groupID#protID
	my $score = $tmp[3];

	unless($fas{$id}){
		$fas{$id} = $score."#".$tmp[2];
	} else {
		$fas{$id} .= ";".$score."#".$tmp[2];
	}
}

# ### open domain file
# open(DOMAIN,$domainFile) || die "Cannot open $domainFile!!\n";
# my @domainIN = <DOMAIN>;
# close (DOMAIN);
# my %domain;
# open(NEW,">$domainFile.mod");
# foreach my $domainLine(@domainIN){
# 	chomp($domainLine);
# 	my @tmp = split(/\t/,$domainLine);	# OG_1001#annca_339:H312_00043#annca_339:H312_01361	annca_339:H312_01361	pfam_DDE_Tnp_IS1595	39	124	1.0
#
# 	unless($tmp[0]){
# 	#print "HERE:$domainLine-$tmp[0]###";<>;
# 	}
# 	elsif(scalar @tmp ne 6){
# 	# print $domainLine;<>;
# 	}
# 	else {print NEW $domainLine,"\n";
#
# 		if(length($tmp[0]) > 2){
# 			unless($domain{$tmp[0]}){
# 				$domain{$tmp[0]} = $domainLine;
# 			} else {
# 				$domain{$tmp[0]} .= "\n".$domainLine;
# 			}
# 		}
# 	}
# }
# print "PARSING DONE!\n";<>;
unless(-d "$domainFol/mDomain_files"){
	mkdir("$domainFol/mDomain_files");
}

### calculate max/mean FAS for orthologs
my %mScore; # $mScore{GroupID_ProtID} = max/mean FAS

open(MEAX,">$orthoFile.mFAS");	# used to store #groupID #protID #mScore
# open(MDOMAIN,">$orthoFile.mDomains");
# print MDOMAIN "seedID\torthoID\tfeature\tstart\tend\tweight\n";

my $c = 0;
foreach my $id(sort keys %fas){
	my @id = split(/#/,$id);
	my @allScore = split(/;/,$fas{$id});	# FAS;FAS;FAS

	# calculate max/mean
	my $meax = 0;
	my $mID = "";
=mean
	if($calc =~ /mean/i){
		$meax = sum(@allScore)/scalar(@allScore);
	}
=cut
	if($calc =~ /max/i){
		foreach (@allScore){
#			print "here:",$_;<>;
			my @scoreTMP = split(/\#/,$_);
#print $meax," - ",$scoreTMP[0];<>;
unless(looks_like_number($scoreTMP[0])){print "NOT NUMBER $_ - $scoreTMP[0]\n";<>;}
			if($meax < $scoreTMP[0]){
				$meax = $scoreTMP[0];
				$mID = $id."#".$scoreTMP[1];
			}
		}
#		print $meax," - ",$dist,"\n";<>;
	}

	# save to mScore output
	# print $meax," - ",$mID;<>;
	my $seedID = "NA";
	if($meax == 0){$meax = '0.0';}
	else {
		my @mIDtmp = split(/#/,$mID);
		$seedID = $mIDtmp[@mIDtmp-1];
	}
	print MEAX $id[0],"\t",$id[1],"\t",$meax,,"\t",$seedID,"\n";

	# save max_architecture into MDOMAIN
	if(length($mID) >2){
		# if($domain{$mID}){
		# 	print MDOMAIN $domain{$mID},"\n";
		# } else {
		# 	print "NO architectures for $mID\n";<>;
		# }
		my @idTMP = split(/#/,$mID);
		open(DOMAIN,"$domainFol/$idTMP[0].domains") || die "Cannot open $domainFol/$idTMP[0].domains\n";
		open(mDOMAIN,">>$domainFol/mDomain_files/$idTMP[0].domains");
		foreach my $line(<DOMAIN>){
			if($line =~ /$mID\t/){
				# print $line,"\n";
				my @lineTMP = split(/#/,$line);
				my @info = split(/\t/,$lineTMP[2]); shift(@info); my $info = join("\t",@info);
				my $newLine = $lineTMP[0]."#".$lineTMP[1]."\t".$info;
				print mDOMAIN $newLine;
			}
		}
	}
	$c++;
	print $c,"/",scalar(keys %fas),"\n";
}
close MEAX;
close (mDOMAIN);
close (DOMAIN);

print "DONE! Check results in:\n\t$orthoFile.mFAS\n\t*.domains files in $domainFol/mDomain_files\n";
exit;
