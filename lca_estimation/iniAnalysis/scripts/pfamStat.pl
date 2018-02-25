#!/usr/bin/env perl
use Cwd;
use Cwd 'abs_path';
use Getopt::Std;
use List::Util qw(sum);

=description
do comparative analysis of pfam annotations between homologous and orphan proteins
Version 1.0 (01.03.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl runCalcPfam.pl -i input_folder\n";
    print "-i\tInput folder\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input folder given\n");

### list of all micros taxa
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

### MAIN
my %homoPfam; my %homoNoPfam;	# $homoPfam{$micros} = len1;len2;len3;... (lenghts of all proteins that have pfam)
my %orphanNoPfam; my %orphanNewPfam; my %orphanHomoPfam;	# s.o.

my $dir = cwd();
#open(LENNOPFAM,">$dir/lenNoPfam_density.txt") || die "Cannot create $dir/lenNoPfam_Density.txt!!\n";

foreach my $micros(@allMicros){
	print $micros,"\n";

	###### HOMOLOGOUS ######
	### get homologous proteins that have and don't have pfam annotations and their sequence lengths
	### and save them to separated output files
	my $homoLen = abs_path($input)."/".$micros.".homologous.LEN";
	my %homoLen = getLen($homoLen);

	my $homoPfam = abs_path($input)."/".$micros.".homologous.PFAM";
	open(HOMOPFAM,$homoPfam) || die "Cannot open $homoPfam!!\n";
	my @homoPfam = <HOMOPFAM>;
	close (HOMOPFAM);

	open(HOMOPFAM,">$homoLen\_Pfam");
	open(HOMONOPFAM,">$homoLen\_noPfam");
	foreach my $homoPfamLine(@homoPfam){
		chomp $homoPfamLine;
		my @homoPfamLineTMP = split(/\t/,$homoPfamLine);	# annca_339:H312_01906	noPFAM

		### save to %homoPfam and %homoNoPfam
		if($homoPfamLineTMP[1] eq "noPFAM"){
			unless($homoNoPfam{$micros}){
				$homoNoPfam{$micros} = $homoLen{$homoPfamLineTMP[0]};
			} else {
				$homoNoPfam{$micros} .= ";".$homoLen{$homoPfamLineTMP[0]};
			}
			print HOMONOPFAM $homoPfamLineTMP[0],"\t",$homoLen{$homoPfamLineTMP[0]},"\n";
		} else {
			unless($homoPfam{$micros}){
				$homoPfam{$micros} = $homoLen{$homoPfamLineTMP[0]};
			} else {
				$homoPfam{$micros} .= ";".$homoLen{$homoPfamLineTMP[0]};
			}
			print HOMOPFAM $homoPfamLineTMP[0],"\t",$homoLen{$homoPfamLineTMP[0]},"\n";
		}
	}

#	my @numberHomoNoPfam = split(/;/,$homoNoPfam{$micros});
#	my @numberHomoPfam = split(/;/,$homoPfam{$micros});
#	print scalar(@numberHomoPfam),"\t",scalar(@numberHomoNoPfam),"\n";<>;
#	print $homoNoPfam{$micros};<>;

	###### ORPHAN ######
	### get number of orphan proteins that have and don't have pfam annotations and their sequence lengths
	### and save them to separated output files
	my $orphanLen = abs_path($input)."/".$micros.".orphan.LEN";
	my %orphanLen = getLen($orphanLen);

	my $orphanPfam = abs_path($input)."/".$micros.".orphan.PFAM";
	open(ORPHANPFAM,$orphanPfam) || die "Cannot open $orphanPfam!!\n";
	my @orphanPfam = <ORPHANPFAM>;
	close (ORPHANPFAM);

	open(ORPHANNOPFAM,">$orphanLen\_noPfam");
	open(ORPHANNEWPFAM,">$orphanLen\_newPfam");
	open(ORPHANHOMOPFAM,">$orphanLen\_homologPfam");
	foreach my $orphanPfamLine(@orphanPfam){
		chomp $orphanPfamLine;
		my @orphanPfamLineTMP = split(/\t/,$orphanPfamLine);	# enche_5516_1:EHEL_050090	newPFAM	PF01227.18

		### save to %orphanNoPfam, %orphanNewPfam and %orphanHomoPfam
		if($orphanPfamLineTMP[1] eq "noPFAM"){
			unless($orphanNoPfam{$micros}){
				$orphanNoPfam{$micros} = $orphanLen{$orphanPfamLineTMP[0]};
			} else {
				$orphanNoPfam{$micros} .= ";".$orphanLen{$orphanPfamLineTMP[0]};
			}
			print ORPHANNOPFAM $orphanPfamLineTMP[0],"\t",$orphanLen{$orphanPfamLineTMP[0]},"\n";
		} elsif($orphanPfamLineTMP[1] eq "newPFAM"){
			unless($orphanNewPfam{$micros}){
				$orphanNewPfam{$micros} = $orphanLen{$orphanPfamLineTMP[0]};
			} else {
				$orphanNewPfam{$micros} .= ";".$orphanLen{$orphanPfamLineTMP[0]};
			}
			print ORPHANNEWPFAM $orphanPfamLineTMP[0],"\t",$orphanLen{$orphanPfamLineTMP[0]},"\n";
		} elsif($orphanPfamLineTMP[1] eq "homologPFAM"){
			unless($orphanHomoPfam{$micros}){
				$orphanHomoPfam{$micros} = $orphanLen{$orphanPfamLineTMP[0]};
			} else {
				$orphanHomoPfam{$micros} .= ";".$orphanLen{$orphanPfamLineTMP[0]};
			}
			print ORPHANHOMOPFAM $orphanPfamLineTMP[0],"\t",$orphanLen{$orphanPfamLineTMP[0]},"\n";
		}
	}

#	my @numberOrphanNoPfam = split(/;/,$orphanNoPfam{$micros});
#	my @numberOrphanNewPfam = split(/;/,$orphanNewPfam{$micros});
#	my @numberOrphanHomoPfam = split(/;/,$orphanHomoPfam{$micros});
#	print scalar(@numberOrphanNoPfam),"\t",scalar(@numberOrphanNewPfam),"\t",scalar(@numberOrphanHomoPfam),"\n";<>;
}
close (HOMOPFAM);
close (HOMONOPFAM);
close (ORPHANNOPFAM);
close (ORPHANNEWPFAM);
close (ORPHANHOMOPFAM);

### print OUTPUT and do analysis
open(PFAMOUT,">$dir/orthomcl_micros_orig.list.pfam") || die "Cannot create $dir/orthomcl_micros_orig.list.pfam!!\n";
print PFAMOUT "Taxon\thomologous_Pfam\thomologous_noPfam\torphan_homologPfam\torphan_newPfam\torphan_noPfam\n";

open(LENMED,">$dir/length_median.stat") || die "Cannot create $dir/length_median.stat!!\n";
open(LENMEAN,">$dir/length_mean.stat") || die "Cannot create $dir/length_mean.stat!!\n";
print LENMED "Type\tTaxon\tmedianLength\n";
print LENMEAN "Type\tTaxon\tmeanLength\n";

foreach my $micros(@allMicros){
	### print number of pfam / noPfam proteins
	my @homoNoPfamLen = split(/;/,$homoNoPfam{$micros});
	my @homoPfamLen = split(/;/,$homoPfam{$micros});

	my @orphanNoPfamLen = split(/;/,$orphanNoPfam{$micros});
	my @orphanNewPfamLen = split(/;/,$orphanNewPfam{$micros});
	my @orphanHomoPfamLen = split(/;/,$orphanHomoPfam{$micros});

	print PFAMOUT $micros{$micros},"\t",scalar(@homoPfamLen),"\t",scalar(@homoNoPfamLen),"\t",scalar(@orphanHomoPfamLen),"\t",scalar(@orphanNewPfamLen),"\t",scalar(@orphanNoPfamLen),"\n";

	### print length output
	print LENMED "homologous_Pfam\t",$micros{$micros},"\t",median(@homoPfamLen),"\n";
	print LENMED "homologous_noPfam\t",$micros{$micros},"\t",median(@homoNoPfamLen),"\n";
	print LENMED "orphan_homologPfam\t",$micros{$micros},"\t",median(@orphanHomoPfamLen),"\n";
	print LENMED "orphan_newPfam\t",$micros{$micros},"\t",median(@orphanNewPfamLen),"\n";
	print LENMED "orphan_noPfam\t",$micros{$micros},"\t",median(@orphanNoPfamLen),"\n";

	print LENMEAN "homologous_Pfam\t",$micros{$micros},"\t",mean(@homoPfamLen),"\n";
	print LENMEAN "homologous_noPfam\t",$micros{$micros},"\t",mean(@homoNoPfamLen),"\n";
	print LENMEAN "orphan_homologPfam\t",$micros{$micros},"\t",mean(@orphanHomoPfamLen),"\n";
	print LENMEAN "orphan_newPfam\t",$micros{$micros},"\t",mean(@orphanNewPfamLen),"\n";
	print LENMEAN "orphan_noPfam\t",$micros{$micros},"\t",mean(@orphanNoPfamLen),"\n";
}
close (PFAMOUT);
close (LENMED);
close (LENMEAN);

### Plot PFAM (column)
open(RSCRIPT,">$dir/pfam.R");
print RSCRIPT "setwd(\"$dir\")
data <- read.table(\"$dir/orthomcl_micros_orig.list.pfam\", sep=\'\\t\',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data\$Taxon <- ordered(data\$Taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))
melted <- melt(data, \"Taxon\")

png(\"$dir/orthomcl_micros_orig.list.pfam.png\", width=1600, height=800)
melted\$cat <- ''
melted[melted\$variable == 'homologous_Pfam',]\$cat <- \"homologous_prots\"
melted[melted\$variable == 'homologous_noPfam',]\$cat <- \"homologous_prots\"
melted[melted\$variable == 'orphan_homologPfam',]\$cat <- \"orphans\"
melted[melted\$variable == 'orphan_newPfam',]\$cat <- \"orphans\"
melted[melted\$variable == 'orphan_noPfam',]\$cat <- \"orphans\"

p = ggplot(melted, aes(x = cat, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    facet_grid(~ Taxon) +
    theme_set(theme_gray(base_size = 18)
)

p = p+labs(x=\"\", y=\"\")

p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())

p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPT);
system("Rscript $dir/pfam.R");

### plot length (line)
open(RSCRIPTMED,">$dir/length_median.R");
print RSCRIPTMED "setwd(\"$dir\")
data <- read.table(\"$dir/length_median.stat\", sep=\'\\t\',header=T)

library(ggplot2)
library(RColorBrewer)

data\$Taxon <- ordered(data\$Taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))

png(\"$dir/length_median.stat.png\", width=1600, height=800)
p = ggplot(data, aes(x=Taxon, y=medianLength, group=Type, colour=Type)) +
  geom_line(colour = \"gray90\") +
  geom_point(size=3)
p = p+labs(x=\"\", y=\"Sequence length (median)\")
p = p+theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPTMED);
system("Rscript $dir/length_median.R");

open(RSCRIPTMEAN,">$dir/length_mean.R");
print RSCRIPTMEAN "setwd(\"$dir\")
data <- read.table(\"$dir/length_mean.stat\", sep=\'\\t\',header=T)

library(ggplot2)
library(RColorBrewer)

data\$Taxon <- ordered(data\$Taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))

png(\"$dir/length_mean.stat.png\", width=1600, height=800)
p = ggplot(data, aes(x=Taxon, y=meanLength, group=Type, colour=Type)) +
  geom_line(colour = \"gray90\") +
  geom_point(size=3)
p = p+labs(x=\"\", y=\"Sequence length (mean)\")
p = p+theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPTMEAN);
system("Rscript $dir/length_mean.R");

exit;

sub getLen{
	my $file = $_[0];
	open(LEN,$file) || die "Cannot open $file!!\n";
	my @allLen = <LEN>;
	close (LEN);
	my %len;
	foreach my $line(@allLen){
		chomp($line);
		my @tmp = split(/\t/,$line);
		$len{$tmp[0]} = $tmp[1];
	}
	return %len;
}

sub median {
	my @a = sort {$a <=> $b} @_;
	my $length = scalar @a;
	return undef unless $length;
	($length % 2)
		? $a[$length/2]
		: ($a[$length/2] + $a[$length/2-1]) / 2.0;
}


sub mean{
	return sum(@_)/scalar(@_);
}
