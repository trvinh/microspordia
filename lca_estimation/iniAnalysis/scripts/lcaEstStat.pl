#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use Getopt::Std;
use IO::Handle;
use File::Path;
use List::Util qw(sum);

=desc
calculate the distribution of each microsporidia taxa in ortholog groups (initial and LCA)
count groups that have 1,2,5,10 or more paralogs of each species (check for whole genome duplication)
can be used to get list of all ORTHOLOGS and ORPHANS proteins for each Micros
2016.02.23
=cut

sub usage {
    my $msg = shift;
    print "example: perl lcaEstStat.pl -i orthomcl_micro_orig.list.NEW -l lca.list\n";
    print "-i\tList of initial ortholog groups\n";
    print "-l\tList of LCA ortholog groups\n";
    die $msg."\n";
}


# global variables
our($opt_i,$opt_l);
getopts('i:l:');

# sanity checks
my $iniFile = ($opt_i) ? $opt_i : usage("ERROR: No initial orthogroups given\n");
my $lcaFile = ($opt_l) ? $opt_l : usage("ERROR: No LCA orthogroups given\n");

### folder contains all micros fasta files
my $faDir = "/share/project/vinh/dataset/micros";
#my $microsName = "/share/project/vinh/dataset/microsTaxa.list";
my $dir = abs_path($iniFile);
my @dir = split(/\//,$dir);
pop(@dir);
$dir = join("/",@dir);

### MAIN
my %orthos;	# contains all proteins that have orthologs of each taxa. e.g. $orthos{enche} = enche:1,enche:2,enche:5,enche:9,...
my %orphan;	# contains orphan proteins of each taxa. e.g. $orthos{enche} = enche:1,enche:2,enche:5,enche:9,...

my %iniProt; # $iniProt{protID} = groupID
my %lcaProt; # $lcaProt{protID} = groupID

### micros list
my @allMicros = ("enche_5516_1","encin_5517_1","enccu_2934_1","nosce_4242_1","entbi_4241_1","vitco_50505","annca_339","antlo_2712","edhae_4124_1","vavcu_5255_1","nempa_5256_1");
my %micros = (
	"enche_5516_1"=>"E.hellem",
	"encin_5517_1"=>"E.intestinallis",
	"enccu_2934_1"=>"E.cuniculi",
	"nosce_4242_1"=>"N.ceranae",
	"entbi_4241_1"=>"E.bieneusi",
	"vitco_50505"=>"V.corneae",
	"annca_339"=>"A.algerae",
	"antlo_2712"=>"A.locustae",
	"edhae_4124_1"=>"E.aedis",
	"vavcu_5255_1"=>"V.culicis",
	"nempa_5256_1"=>"N.parisii"
);

my %paralogs = (
	"enche_5516_1"=>0,
	"encin_5517_1"=>0,
	"enccu_2934_1"=>0,
	"nosce_4242_1"=>0,
	"entbi_4241_1"=>0,
	"vitco_50505"=>0,
	"annca_339"=>0,
	"antlo_2712"=>0,
	"edhae_4124_1"=>0,
	"vavcu_5255_1"=>0,
	"nempa_5256_1"=>0
);
my $max = 0;	# maximal paralogs in one group

### initial orthologs list
open(INI,abs_path($iniFile)) || die "Cannot open $iniFile!!\n";
my @iniList = <INI>;
close (INI);
print "parsing initial orthogroups...";
foreach my $iniLine (@iniList){
	chomp($iniLine);
	my @spec = getSpecList2($iniLine);
	my @iniLineTMP = split(/\t/,$iniLine);	# $iniLineTMP[0] is groupID
#	print $iniLine,"\n","@spec\n";<>;

	# save into %orthos and $iniProt
	foreach my $spec(@spec){
		my @specTMP = split(/;/,$spec);	# $specTMP[0]=taxaName; $specTMP[1]= prot1,prot2,prot3,...
		my $specName = shift(@specTMP);
		my $allParalogs = join(";",@specTMP);

		### save orthologous proteins for this species into %ortho
		unless($orthos{$specName}){
			$orthos{$specName} = $allParalogs;
		} else {
			$orthos{$specName} .= ";".$allParalogs;
		}

		### save initial proteins into %iniProt
		foreach my $prot(@specTMP){
			unless($iniProt{$prot}){
				$iniProt{$prot} = $iniLineTMP[0];
			} else {
				print "$prot in $iniLineTMP[0] already exists in $iniProt{$prot}\n";<>;
			}
		}

		### count number of paralogs for this species in this group and save into %paralogs
		unless($paralogs{$specName}){
			$paralogs{$specName} = scalar(@specTMP)."#".$iniLineTMP[0];
		} else {
			$paralogs{$specName} .= ";".scalar(@specTMP)."#".$iniLineTMP[0];
		}
		if($max < scalar(@specTMP)){
			$max = scalar(@specTMP);
		}
	}
}
print "DONE\n";

### LCA orthologs list
open(LCA,abs_path($lcaFile)) || die "Cannot open $lcaFile!!\n";
my @lcaList = <LCA>;
close (LCA);
print "parsing LCA orthogroups...";
foreach my $lcaLine (@lcaList){
	chomp($lcaLine);
	my @spec = getSpecList2($lcaLine);
	my @lcaLineTMP = split(/\t/,$lcaLine);
#	print $lcaLine,"\n","@spec\n";<>;

	# save into %orthos and $iniProt
	foreach my $spec(@spec){
		my @specTMP = split(/;/,$spec);	# $specTMP[0]=taxaName; $specTMP[1]= prot1,prot2,prot3,...
		my $specName = shift(@specTMP);
		my $allParalogs = join(";",@specTMP);

		### save initial proteins into %lcaProt
		foreach my $prot(@specTMP){
			unless($lcaProt{$prot}){
				$lcaProt{$prot} = $lcaLineTMP[0];
			} else {
				print "$prot in $lcaLineTMP[0] already exists in $lcaProt{$prot}\n";<>;
			}
		}
	}
}
print "DONE\n";
print "Max paralogs count = $max\n";

### COUNTING
my $outFile =abs_path($iniFile).".stat"; 
open(OUT,">$outFile");
print OUT "Taxon\tLCA_homologous_prots\tNon-LCA_homologous_prots\tOrphan_prots\tLCA OG\tNon-LCA OG\n";

my $paraOut = abs_path($iniFile).".para";
open(PARA,">$paraOut");
print PARA "Taxon\tParalogs\tGroups\n";
open(PARAID,">$paraOut"."ID");

foreach my $micros(@allMicros){
	print $micros,"\n";
	my $iniGroup = 0;	# number of initial orthogroups
	my $lcaGroup = 0;	# number of LCA orthogroups
	my $allGene = 0;	# total number of proteins of each taxon
	my $orthoLCA = 0;	# number of LCA orthologous proteins
	my $orthoIni = 0;	# number of initial orthologous proteins

	### get number of LCA orthologous proteins
	my @allOrtho = split(/;/,$orthos{$micros});	# number of homologous proteins
	foreach my $ortho(@allOrtho){
		if($lcaProt{$ortho}){
			$orthoLCA ++;
		} else {
			$orthoIni ++;
		}
	}

	### get all prot IDs from fasta file
	my @allFA = getIDfasta("$faDir/$micros\.fa");

	### calc number of ini groups and lca groups for each species
	my %iniGroups; my %lcaGroups;
	foreach my $seq (@allFA){
		if(length($seq) > 0){
			my @seqTMP = split(/\n/,$seq);
			my $seqID = shift(@seqTMP);
#			print $seqID;<>;

			if($lcaProt{$seqID}){
				$lcaGroups{$lcaProt{$seqID}} = 1;
			}

			if($iniProt{$seqID}){
				$iniGroups{$iniProt{$seqID}} = 1;
			}
			### get orphan proteins for this species
			else {
				unless($orphan{$micros}){
					$orphan{$micros} = $seqID;
				} else {
					$orphan{$micros} .= ";".$seqID;
				}
			}
			$allGene ++;
#			print "\n";<>;
		}
	}

#	print "orthogene = ",scalar(@allOrtho),"\nLCA ortho gene = ",$orthoLCA,"\nIni ortho gene = ",$orthoIni,"\nini groups = ",scalar(keys %iniGroups),"\nlca groups = ",scalar(keys %lcaGroups),"\n";<>;
	print OUT "$micros{$micros}\t",$orthoLCA,"\t",$orthoIni,"\t",$allGene-scalar(@allOrtho),"\t",scalar(keys %lcaGroups),"\t",scalar(keys %iniGroups)-scalar(keys %lcaGroups),"\n";

	### calc number of groups, in which this species has 1,2,3,...11...$max paralogs
#	print $paralogs{$micros},"\n";
	my @paralogs = split(/;/,$paralogs{$micros});
	my %stat;	# $stat{1} = OG_1001;OG_1003;...
	foreach my $para(@paralogs){
#		print $para;<>;
		my @paraTMP = split(/#/,$para);
		unless($stat{$paraTMP[0]}){
			$stat{$paraTMP[0]} = $paraTMP[1];
		} else {
			$stat{$paraTMP[0]} .= ";".$paraTMP[1];
		}
	}

	print PARAID "### $micros{$micros}\n";
	my $one = 0; my $two = 0; my $five = 0; my $ten = 0; my $tweenty = 0; my $more = 0;
	foreach(sort {$a<=>$b} keys %stat){
		my @tmp = split(/;/,$stat{$_});
#		print "$_ paralogs (in ",scalar(@tmp)," groups): $stat{$_}\n";
		print PARAID "$_ paralogs (in ",scalar(@tmp)," groups): $stat{$_}\n";
		if($_ eq 1){$one = scalar(@tmp);}
		elsif($_ eq 2){$two = scalar(@tmp);}
		elsif($_ > 2 && $_ < 6){$five += scalar(@tmp);}
		elsif($_ > 5 && $_ < 11){$ten += scalar(@tmp);}
		elsif($_ > 10 && $_ < 21){$tweenty += scalar(@tmp);}
		elsif($_ > 20){$more += scalar(@tmp);}
	}

	print PARA $micros{$micros},"\t1\t",$one,"\n";
	print PARA $micros{$micros},"\t2\t",$two,"\n";
	print PARA $micros{$micros},"\t3-5\t",$five,"\n";
	print PARA $micros{$micros},"\t6-10\t",$ten,"\n";
	print PARA $micros{$micros},"\t11-20\t",$tweenty,"\n";
	print PARA $micros{$micros},"\t>20\t",$more,"\n";
#	<>;
}

### Plot 1 (column, number of proteins and orthogroups)
open(RSCRIPT,">$dir/stat.R");
print RSCRIPT "setwd(\"$dir\")
data <- read.table(\"$outFile\", sep=\'\\t\',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data\$Taxon <- ordered(data\$Taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))
melted <- melt(data, \"Taxon\")

png(\"$outFile.png\", width=1600, height=800)
melted\$cat <- ''
melted[melted\$variable == 'LCA_homologous_prots',]\$cat <- \"#Genes\"
melted[melted\$variable == 'Non.LCA_homologous_prots',]\$cat <- \"#Genes\"
melted[melted\$variable == 'Orphan_prots',]\$cat <- \"#Genes\"
melted[melted\$variable == 'LCA.OG',]\$cat <- \"#Groups\"
melted[melted\$variable == 'Non.LCA.OG',]\$cat <- \"#Groups\"

p = ggplot(melted, aes(x = cat, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    facet_grid(~ Taxon) +
    theme_set(theme_gray(base_size = 15))
p = p+labs(x=\"\", y=\"\")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPT);
system("Rscript $dir/stat.R");

### plot 2 (line, paralogs)
open(RSCRIPT2,">$dir/para.R");
print RSCRIPT2 "setwd(\"$dir\")
data <- read.table(\"$paraOut\", sep=\'\\t\',header=T)

library(ggplot2)
library(RColorBrewer)

data\$Paralogs <- ordered(data\$Paralogs, levels = c(\"1\",\"2\",\"3-5\",\"6-10\",\"11-20\",\">20\"))

png(\"$paraOut.png\", width=1600, height=800)
p = ggplot(data, aes(x=Paralogs, y=Groups, group=Taxon, colour=Taxon)) +
  geom_line() +
  geom_point()
p = p+labs(x=\"Number of paralogs\", y=\"Number of homologous groups\")
p = p+theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),legend.text=element_text(size=15))
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPT2);
system("Rscript $dir/para.R");

### print list of orphan and homologous proteins for each species
unless(-d "$dir/homolog_orphan_dir"){
	mkdir("$dir/homolog_orphan_dir");
}
foreach my $spec (keys %orthos){
#	print $micros{$spec};<>;
	open(ORTHO,">$dir/homolog_orphan_dir/$spec.homologous") || die "Cannot create $dir/homolog_orphan_dir/$spec.homologous!!\n";
	my $allOrtho = $orthos{$spec}; $allOrtho =~ s/;/\n/g;
	print ORTHO $allOrtho;
	close (ORTHO);

	open(ORPHAN,">$dir/homolog_orphan_dir/$spec.orphan") || die "Cannot create $dir/homolog_orphan_dir/$spec.orphan!!\n";
	my $allOrphan = $orphan{$spec}; $allOrphan =~ s/;/\n/g;
	print ORPHAN $allOrphan;
	close (ORPHAN);
}

### FINISH
print "FINISHED!! Check results in\n\t$iniFile.stat\n\t$iniFile.stat.png\n\t$iniFile.para\n\t$iniFile.para.png\n";
close (OUT);
close (PARA);
close (PARAID);
exit;


### get all species of one ortholog group
### v2: together with their all paralogs
sub getSpecList2 {
	my $line = $_[0];
	my %spec_list;
	my @tmp = split(/\t/,$line);
	shift (@tmp);
	foreach my $t(@tmp){
		my @tmp2 = split(/:/,$t);
		unless($spec_list{$tmp2[0]}){
			$spec_list{$tmp2[0]} = $t;
		} else {
			$spec_list{$tmp2[0]} .= ";".$t;
		}
	}
	my @spec = keys(%spec_list);
	my @results; # (specName;#paralogs)
	foreach (@spec){
		my $out = $_.";".$spec_list{$_};
		$out =~ s/,$//;
		push(@results,$out);
	}
	return @results;	# $results[0] = SPECIES;prot1;prot2;prot3
}


### get all protein IDs of a fasta file
sub getIDfasta{
	my $file = $_[0];

	open(FILE,$file) || die "Cannot open $file !!\n";
	my @file = <FILE>;
	close (FILE);

	my %id_list;
	foreach my $line(@file){
		if($line =~ /^>/){
			chomp($line);
			$line =~ s/>//;
			$id_list{$line} = 1;
		}
	}

	return keys(%id_list);
}
