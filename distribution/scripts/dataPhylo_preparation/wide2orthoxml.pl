#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

# convert a wide format phyloprofile file into orthoXML format
# 27.06.2017

sub usage {
    my $msg = shift;
    print "example: perl wide2orthoxml.pl -i demo/test.main -n taxonNamesReduced.txt -a newTaxa.txt -f FAS -s traceability\n";
    print "-i\tWide-format phyloprofile input\n";
    print "-n\tFile contains all NCBI IDs, their names and ranks (reduced file)\n";
    print "-a\tFile contains newly added taxa IDs, their names and ranks\n";
		print "-f\tscore 1\n";
		print "-s\tscore 2\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_n,$opt_a,$opt_f,$opt_s);
getopts('i:n:a:f:s:');

my $fileIn = ($opt_i) ? $opt_i : usage("ERROR: No input phyloprofile file given\n");
my $nameIN =  ($opt_n) ? $opt_n : usage("ERROR: No taxonNamesFull file given\n");
my $nameNewIN = ($opt_a) ? $opt_a : usage("ERROR: No newly added taxa file given\n");
my $stScore = ($opt_f) ? $opt_f : usage("ERROR: No name for first score given\n");
my $ndScore = ($opt_s) ? $opt_s : usage("ERROR: No name for second score given\n");

#### get list of all species name and their NCBI taxon id
my %id;		# $id{NAME} = ID
my %name;	# $name{$id} = NAME
my %rank;	# $rank{$id} = RANK
my %parent;	# $parent{$id} = PARENT_ID

open(NAME,"$nameIN") || die "Cannot open $nameIN!!\n";
foreach my $line(<NAME>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# id  name  rank  parentID
	$id{$tmp[1]} = $tmp[0];
	$name{$tmp[0]} = $tmp[1];
	$rank{$tmp[0]} = $tmp[2];
	$parent{$tmp[0]} = $tmp[3];
}
close (NAME);

#### get the info for newly added taxa (if necessary)
open(NEW,"$nameNewIN") || die "Cannot open $nameNewIN!!\n";
my @nameNEW = <NEW>;
close (NEW);

if(scalar @nameNEW > 1){
	foreach my $line(@nameNEW){
		chomp($line);
    if(length $line > 2){
      my @tmp = split(/\t/,$line);	# id  name  rank  parentID
      $id{$tmp[1]} = $tmp[0];
      $name{$tmp[0]} = $tmp[1];
      $rank{$tmp[0]} = $tmp[2];
      $parent{$tmp[0]} = $tmp[3];
    }
	}
}
print "Parsing taxonNamesFull.txt done!!\n";

#### get species IDs and their genes
my %specGenes; # $specGenes{$specID} = gene1;gene2;gene3;...

open(IN,$fileIn) || die "cannot open $fileIn!!\n";
my @file = <IN>;
close (IN);

chomp(my $header = shift(@file));
my @allSpec = split(/\t/,$header);

foreach my $line(@file){
	chomp($line);
	my @ortho = split(/\t/,$line);
	for(my $i=1; $i < scalar(@ortho); $i++){
		unless($ortho[$i] =~ /^NA/){
			my @tmp = split(/#/,$ortho[$i]);	# aedae_350_1:16034#0.99980944#0.935407215
			unless($specGenes{$allSpec[$i]}){
				$specGenes{$allSpec[$i]} = ";".$tmp[0].";";
			} else {
				unless($specGenes{$allSpec[$i]} =~ /;$tmp[0];/){
					$specGenes{$allSpec[$i]} .= $tmp[0].";";
				}
			}
		}
	}
}

###############################
####### WRITE XML FILE ########
###############################
open(OUT,">$fileIn.xml") || die "cannot create $fileIn.xml\n";

##### HEADER
print OUT "<?xml version=\"1.0\" encoding=\"utf-8\"?>
<orthoXML xmlns=\"http://orthoXML.org/2011/\" version=\"0.3\" origin=\"phyloprofile\"
  originVersion=\"1.0\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"
  xsi:schemaLocation=\"http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd\">

  <notes>
    OrthoXML file converted from PhyloProfile App.
  </notes>\n\n";

##### SPECIES
# <species name="Caenorhabditis elegans" NCBITaxId="6239">
# 	<database name="WormBase" version="Caenorhabditis-elegans_WormBase_WS199_protein-all.fa"
# 		<genes>
# 			<gene id="1" geneId="WBGene00000962" protId="CE23997" />
# 			<gene id="5" geneId="WBGene00006801" protId="CE43332" />
# 		</genes>
# 	</database>
# </species>

foreach my $specID(sort keys %specGenes){
	my $specIDmod = $specID; $specIDmod =~ s/ncbi//;
unless($name{$specIDmod}){
	print "NO NAME FOR $specIDmod\n";<>;
}
	print OUT "  <species name=\"$name{$specIDmod}\" NCBITaxId=\"$specIDmod\">
    <database name=\"undef\" version=\"undef\">
      <genes>\n";
	my @genes = split(/;/,$specGenes{$specID});
	for(@genes){
		if(length($_) > 0){
			print OUT "        <gene id=\"$_\" protId=\"$_\" />\n"
		}
	}
  print OUT "      </genes>
    </database>
  </species>\n\n";
}

##### SCORES
# <scores>
# 	<scoreDef id="bit" desc="BLAST score in bits of seed orthologs" />
# 	<scoreDef id="inparalog" desc="Distance between edge seed ortholog" />
# 	<scoreDef id="bootstrap" desc="Reliability of seed orthologs" />
# </scores>
print OUT "  <scores>
    <scoreDef id=\"$stScore\" desc=\"$stScore\" />
    <scoreDef id=\"$ndScore\" desc=\"$ndScore\" />
  </scores>\n\n";

##### ortholog groups
# <groups>
# 	<orthologGroup id="3">
# 	  <geneRef id="5">
# 	    <score id="inparalog" value="1" />
# 	    <score id="bootstrap" value="1.00" />
# 	  </geneRef>
# 	  <geneRef id="6">
# 	    <score id="bootstrap" value="1.00" />
# 			<score id="inparalog" value="1" />
# 	  </geneRef>
# 	  <geneRef id="7">
# 	    <score id="inparalog" value="0.4781" />
# 	  </geneRef>
# 	</orthologGroup>
# 	...
# </groups>

print OUT "  <groups>\n";
foreach my $line(@file){
	chomp($line);
	if(length($line)>1){
		my @tmp = split(/\t/,$line);	# OG_1001  gene1#fas#tracae gene3#fas#trace ...

		my $groupID = shift(@tmp);
		print OUT "    <orthologGroup id=\"$groupID\">\n";

		foreach my $prot(@tmp){
			my @protTMP = split(/#/,$prot);
			unless($protTMP[0] eq "NA"){
				print OUT "      <geneRef id=\"$protTMP[0]\">
				<score id=\"$stScore\" value=\"$protTMP[1]\" />
				<score id=\"$ndScore\" value=\"$protTMP[2]\" />
			</geneRef>\n";
			}
		}
		print OUT "    </orthologGroup>\n";
	}
}

print OUT "  </groups>\n";

##### FOOTER
print OUT "\n</orthoXML>";

exit;
