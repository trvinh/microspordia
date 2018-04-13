#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

## do KO annotation for list of orthologs
## v2.0 (2015.06.15)
## weight KO annotation for LCA proteins based on taxonomy distance (micros > fungi > metazoa > ...)
# v3.0 (21.10.2015)
# add FAS score for each orthologs
# v4.0 (24.05.2016)
# add FAS score for each orthologs
# anno stat for each taxon
# v5.0 (11.07.2016)
# replace taxonomy distance by patristic distance
# v6.0 (21.09.2016)
# taking FAS cutoff into account

sub usage {
    my $msg = shift;
    print "example: perl koAnnotation.pl -i orthologs.list -f orthologs.list.meanFAS -d PatristicDistance folder -c koThreshold.txt\n";
    print "-i\tinput list of orthologs\n";
    print "-f\tinput file contains (max/mean) FAS scores of all orthologs\n";
    print "-d\tfolder contains distance results\n";
    print "-c\tKO-ref threshold file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_f,$opt_c,$opt_d);
getopts('i:f:c:d:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $cutoffIn = ($opt_c) ? $opt_c : usage("ERROR: No KO cutoff file given\n");
my $paDistFol = ($opt_d) ? $opt_d : usage("ERROR: No distance folder given\n");

my $koPath = "/home/vinh/Desktop/data/project/KEGG_annotation/koMap";
my $refseqPath = "/home/vinh/Desktop/data/project/KEGG_annotation/refseq";
my $uniprotPath = "/home/vinh/Desktop/data/project/KEGG_annotation/uniprot";
#my $paDistFol = "/home/vinh/Desktop/data/project/patristicDistCalc/out";

### open and parse FAS output file
#my $qualCheckFile = "/home/vinh/Desktop/data/orthologs/hamstr/results_20150701_FINAL/distribution/rozal/qualCheck/qualCheck.list";
my $qualCheckFile = ($opt_f) ? $opt_f : usage("ERROR: No FAS scores file given\n");
open(QUAL,$qualCheckFile) || die "Cannot open $qualCheckFile!!\n";
my @qual = <QUAL>;
close (QUAL);

my %fas;	# $fas{groupID#protID} = FAS;
foreach my $qualLine(@qual){
	chomp($qualLine);
	my @tmp = split(/\t/,$qualLine);	# OG_1000	mja:MJ1072	0.170314283823513
	my $id = $tmp[0]."#".$tmp[1];
	$fas{$id} = $tmp[2];
}

### open and parse KO threshold file
open(CUTOFF,$cutoffIn) || die "Cannot open $cutoffIn!\n";
my @cutoff = <CUTOFF>;
close (CUTOFF);

my %cutoff;
foreach my $line(@cutoff){
	chomp($line);
	my @tmp = split(/\t/,$line);	# K00011  0.999109342281437
	$cutoff{$tmp[0]} = $tmp[1];
}

### open list of ortholog
open(IN,$input) || die "Cannot open $input!!\n";
my @in=<IN>;
close (IN);

my %lackUniProt;	# used to store protein which cannot be found in uniprotKB file (it should not be happened)
## output files
open(OUT,">$input.KO") || die "Cannot create $input.KO!!\n";
open(SUMMARY,">$input.KO.SUMMARY") || die "Cannot create $input.KO.SUMMARY!!\n";
open(LOG,">$input.KO.LOG") || die "Cannot create $input.LOG!!\n";
open(STAT,">$input.taxonKO") || die "Cannot create $input.taxonKO!!\n";		### KO list for each taxon
my %taxonKO;	# $taxonKO{taxonID} = KOlist

=level
# taxonomy level
my %level = (
	"enche_5516_1" => 1,"encin_5517_1" => 1,"enccu_2934_1" => 1,"nosce_4242_1" => 1,"entbi_4241_1" => 1,"antlo_2712" => 1,"edhae_4124_1" => 1,"vavcu_5255_1" => 1,"nempa_5256_1" => 1,"annca_339" => 1,"vitco_50505" => 1,
	"lacbi_2940_1" => 2,"pucgr_1921_1" => 2,
	"schpo_2340_1" => 3,
	"neucr_1906" => 4,"aspni_2095_1" => 4,
	"sacce_2336_1" => 5,"ago" => 5,"canal_2931_1" => 5,

	"hsa" => 6,"mmu" => 6,"rno" => 9,"dre" => 6,"dme" => 6,"cel" => 6,"nemve_2309_1" => 6,"ampqu_4652_1" => 6,
	"monbr_1569_1" => 6,

	"ehi" => 7,

	"trybr_3330_1" => 8,"ngr" => 8,
	"ath" => 8,"cre" => 8,
	"psoj" => 8,"pfa" => 8,"cho" => 8,

	"mja" => 9,"ape" => 9,

	"eco" => 10,"nme" => 10,"hpy" => 10,"bsu" => 10,"lla" => 10,"mge" => 10,"mtu" => 10,"syn" => 10,"aae" => 10
);
=cut
my $allmicrosID = "enche_5516_1;encin_5517_1;enccu_2934_1;nosce_4242_1;entbi_4241_1;vitco_50505;annca_339;antlo_2712;edhae_4124_1;vavcu_5255_1;nempa_5256_1";

my $c=0;
foreach my $line(@in){
	### process each ortholog group (each line)
	if(length($line)>2){
		chomp($line);
#print $line,"\n";
		my @items = split(/\t/,$line);

		### get groupID
		my $groupID = shift(@items);
	#	print $groupID,"\n";
		my %koMember;	# ko for each protein
#		my %koGroup;	# ko(s) for the whole group (LCA protein)
#		my %koFas;	# max FAS score for each groupKO
#    my %koDist;	# min distance for each groupKO
    my %groupKO;  # $groupKO{koID} = minDist#FAS

		### open patristic distance file, get max/min distance (for scaling) and all pairwise distance betwenn non-micros and micros species
		## (for each non-micros taxon, get the smallest dist between that taxon vs all micros taxa)
		my %paDist;
		my $flag = 1;	# $flag = 0 if groupID.dist file not exists
		if(-e "$paDistFol/$groupID.dist"){
			open(DIST, "$paDistFol/$groupID.dist");
			my @dist = <DIST>;
			close (DIST);

			if(scalar @dist > 1){
				foreach my $line(@dist){
					chomp($line);
					my @lineTMP = split(/\t/,$line);	# annca_339:H312_02162	edhae_4124_1:EDEG_04217	1.329581852

					### get species name for 1.item and 2.item
					my @first = split(/:/,$lineTMP[0]); my $firstID = shift(@first);
					my @second = split(/:/,$lineTMP[1]); my $secondID = shift(@second);

					### save first item to %paDist, only if its partner is microsporidia prot
					if($allmicrosID =~ /$secondID/){
						unless($paDist{$lineTMP[0]}){
							$paDist{$lineTMP[0]} = $lineTMP[2];
						} else {
							## update distance of this protein if the existing one is larger
							if($paDist{$lineTMP[0]} > $lineTMP[2]){
								$paDist{$lineTMP[0]} = $lineTMP[2];
							}
						}
					}

					### save second item to %paDist, only if its partner is microsporidia prot
					if($allmicrosID =~ /$firstID/){
						unless($paDist{$lineTMP[1]}){
							$paDist{$lineTMP[1]} = $lineTMP[2];
						} else {
							## update distance of this protein if the existing one is larger
							if($paDist{$lineTMP[1]} > $lineTMP[2]){
								$paDist{$lineTMP[1]} = $lineTMP[2];
							}
						}
					}
				}
			} else {
				$flag = 0;
			}
		} else {
			$flag = 0;
		}

    ## get min and max distances of this group
    my @allValues = sort {$a<=>$b} values %paDist;
    my $minDist = shift @allValues;
    my $maxDist = pop @allValues;

		### read each ortholog and get KO
		foreach my $prot(@items){
			### get species name
			my @tmp = split(/:/,$prot);
			my $specID = $tmp[0];	# species name
			my $protID = $tmp[1];	# protein ID (ortholog ID)
	#		print $specID,"\n";
			my $koID = ""; my $ecID = "";

			if(-e "$koPath/$specID.ko"){
				### get KO
				open(KO,"$koPath/$specID.ko") || die "Cannot open $koPath/$specID.ko !!\n";
				my @koList = <KO>;
				close (KO);
				my $koList = "\n"; $koList = join("",@koList); $koList .= "\n";

				if($koList =~ /\:$protID\tK\d{5}/){
					my $hit = $&;
					my @hit = split(/\t/,$hit);	# ape:APE_0166.1  K00798  2.5.1.17

					$koID = $hit[1];
				}
			} else {
				### DO SOMETHING HERE....
#				print "NO KO FILE!!\n";
			}

			## save KO and EC for each ortholog in this group into %koMember
			unless($koMember{$prot}){
				if(length($koID)>0){
					$koMember{$prot} = $koID;
				} else {
					$koMember{$prot} = "";
				}
			} else {
				unless($koMember{$prot} =~ /$koID/){
					print "THIS $prot already saved but does not contains $koID!!\n";<>;
				}
			}

			## save KO for this taxon into %taxonKO
			unless($taxonKO{$specID}){
				if(length($koID)>0){
					$taxonKO{$specID} = $koID.";";
				}
			} else {
				if(length($koID)>0 and $taxonKO{$specID} !~ /$koID\;/){
					$taxonKO{$specID} .= $koID.";";
				}
			}

			## save KO for group
			my $cutoff = 0.0;
			if($cutoff{$koID}){ $cutoff = $cutoff{$koID};}

      unless($fas{"$groupID#$prot"}){
      	print "NO FAS FOR $groupID # $prot\n";<>;
      } else {
        if(length($koID)>0 and $fas{"$groupID#$prot"} >= $cutoff){
  				my $scaledDist = "NA";
          if(defined $paDist{$prot}){
            $scaledDist = ($paDist{$prot}-$minDist)/($maxDist - $minDist);;
          }
  				if($scaledDist eq 0){$scaledDist = "0.0"}
#         print $prot," - ",$paDist{$prot}," - ",$scaledDist;<>;
#          if($scaledDist lt 0){print "HERE: $prot\tpa=$paDist{$prot}\tmin=$minDist\tmax=$maxDist\t=>\t$scaledDist";<>;}

          unless($groupKO{$koID}){
            $groupKO{$koID} = $fas{"$groupID#$prot"}."#".$scaledDist;
          } else {
            my @groupKOtmp = split(/#/,$groupKO{$koID});
            if($groupKOtmp[1] gt $scaledDist){
              $groupKO{$koID} = $fas{"$groupID#$prot"}."#".$scaledDist;
            } elsif($groupKOtmp[1] eq $scaledDist){
              if($groupKOtmp[0] < $fas{"$groupID#$prot"}){
                $groupKO{$koID} = $fas{"$groupID#$prot"}."#".$scaledDist;
              }
            }
          }
  			}
      }
		}

		### print OUTPUT
		print OUT "### ",$groupID,"\t";
		print SUMMARY $groupID,"\t";
		my @groupKOs;
		foreach (sort keys(%groupKO)){
			my $kotmp = $_."#".$groupKO{$_};
			push (@groupKOs,$kotmp);
		}

		my $groupKOs = join("\t",@groupKOs); $groupKOs =~ s/^\t//; $groupKOs =~ s/\t$//; $groupKOs =~ s/\t\t+/\t/g;
		print OUT $groupKOs,"\n";
		print SUMMARY $groupKOs,"\n";
		foreach my $ortho (sort(keys(%koMember))){
			my $fas = $fas{"$groupID#$ortho"};
			my $scaledDist = "NA";
			if($paDist{$ortho}){
				$scaledDist = ($paDist{$ortho}-$minDist)/($maxDist - $minDist);
			}

			if($fas){
#				print $ortho,"\tKO=",$koMember{$ortho},"#####\tmFAS=",$fas,"\tCutoff=",$cutoff{$koMember{$ortho}},"\tDist=",$scaledDist,"\n";
				print OUT $ortho,"\t",$koMember{$ortho},"\t",$fas,"\t",$scaledDist,"\n";
			} else {
#				print $ortho,"\tKO=",$koMember{$ortho},"\tFAS=0","\tCutoff=",$cutoff{$koMember{$ortho}},"\tDist=",$scaledDist,"\n";
				print OUT $ortho,"\t",$koMember{$ortho},"\t","\t",$scaledDist,"\n";
			}
#			<>;
		}
		print OUT "\n";

#		print "next group";<>;
		$c++;
		print "$c\n";
	}
}
close (OUT);
close (SUMMARY);

### print KO list for each taxon
foreach my $taxon(sort keys %taxonKO){
	print STAT $taxon,"\t",$taxonKO{$taxon},"\n";
}
close (STAT);

### print LOG file
foreach my $lackUniSpec (keys(%lackUniProt)){
	print LOG $lackUniSpec,"\t",$lackUniProt{$lackUniSpec},"\n";
}
close (LOG);
exit;

sub getKO {
	my ($keggID,$koList) = @_;
	my $koID = ""; my $ecID = "";
	if($koList =~ /\n$keggID(.)*?\n/){
		my @koHit = split(/\t/,$&);
		$koID = $koHit[1];
		$ecID = $koHit[2]; $ecID =~ s/\n//g;
	#	print $koID," - ",$ecID,"\n";
	}
	return ($koID,$ecID);
};
