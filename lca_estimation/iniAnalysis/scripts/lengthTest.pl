#!/usr/bin/env perl
use Cwd;
use Cwd 'abs_path';
use Getopt::Std;
use List::Util qw(sum);

=description
do Mann Whitney U-test (Mann-Whitney-Wilcoxon Test) for comparing lengths between homologous and orphan proteins
*** T-test is used only for normal distributed data
Version 1.0 (02.03.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl lengthTest.pl -i /home/vinh/Desktop/data/project/iniAnalysis/results/homolog_orphan_dir\n";
    print "-i\tFolder contains length files\n";
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
my $dir = cwd();
unless(-d "$dir/lengthTest"){
	mkdir("$dir/lengthTest");
}
open(GENERAL,">$dir/lengthTest/homologous_vs_orphan.LEN.pvalue");	# stat for all proteins (with and without pfam)
print GENERAL "Taxon\thomologous_Mean\torphan_Mean\tP_value\n";
open(DIVIDENOPFAM,">$dir/lengthTest/noPfam_homologous_vs_orphan.LEN.pvalue");	# divide into pfam and noPfam proteins
print DIVIDENOPFAM "Taxon\thomologous_Mean\torphan_Mean\tP_value\n";
open(DIVIDEPFAM,">$dir/lengthTest/Pfam_homologous_vs_orphan.LEN.pvalue");	# divide into pfam and noPfam proteins
print DIVIDEPFAM "Taxon\thomologous_Mean\torphan_Mean\tP_value\n";

open(GENERALR,">$dir/lengthTest/homologous_vs_orphan.LEN.list");	# used for plotting
print GENERALR "type\ttaxon\tlength\tpvalue\n";
open(DIVIDENOPFAMR,">$dir/lengthTest/noPfam_homologous_vs_orphan.LEN.list");	# used for plotting
print DIVIDENOPFAMR "type\ttaxon\tlength\tpvalue\n";
open(DIVIDEPFAMR,">$dir/lengthTest/Pfam_homologous_vs_orphan.LEN.list");	# used for plotting
print DIVIDEPFAMR "type\ttaxon\tlength\tpvalue\n";

foreach my $micros(@allMicros){
	print $micros,"\n";

	#################### FOR ALL CASES (with and without PFAM annotations) ############################
	### get length of homologous and orphan proteins for this species
	my $homoFile = abs_path($input)."/".$micros.".homologous.LEN";
	chomp(my $homoLen = `cut -f 2 $homoFile`);
	$homoLen =~ s/\n/,/g;
#	print $homoLen;<>;
	my $orphanFile = abs_path($input)."/".$micros.".orphan.LEN";
	chomp(my $orphanLen = `cut -f 2 $orphanFile`);
	$orphanLen =~ s/\n/,/g;
#	print $orphanLen;<>;

	### perform T-test
	my $ttest_out = parse_utest($homoLen,$orphanLen,$dir);
	print GENERAL $micros{$micros},"\t",$ttest_out,"\n";

	### create R file for plotting
	my @ttest_out = split(/\t/,$ttest_out);
	my $pvalue = $ttest_out[@ttest_outTMP-1];

	my @homoLen = split(/,/,$homoLen);
	my @orphanLen = split(/,/,$orphanLen);
	foreach my $homoGeneral(@homoLen){ print GENERALR "homologous\t$micros{$micros}\t$homoGeneral\t$pvalue\n";}
	foreach my $orphanGeneral(@orphanLen){ print GENERALR "orphan\t$micros{$micros}\t$orphanGeneral\t$pvalue\n";}
	###################################################################################################

	############################ FOR PROTEINS that DO NOT HAVE PFAM annotations #######################
	### get length of homologous and orphan proteins for this species
	my $homoFileNoPfam = abs_path($input)."/".$micros.".homologous.LEN_noPfam";
	chomp(my $homoLenNoPfam = `cut -f 2 $homoFileNoPfam`);
	$homoLenNoPfam =~ s/\n/,/g;
#	print $homoLenNoPfam;<>;
	my $orphanFileNoPfam = abs_path($input)."/".$micros.".orphan.LEN_noPfam";
	chomp(my $orphanLenNoPfam = `cut -f 2 $orphanFileNoPfam`);
	$orphanLenNoPfam =~ s/\n/,/g;
#	print $orphanLenNoPfam;<>;

	### perform T-test
	$ttest_out = parse_utest($homoLenNoPfam,$orphanLenNoPfam,$dir);
	print DIVIDENOPFAM $micros{$micros},"\t",$ttest_out,"\n";

	### create R file for plotting
	@ttest_out = split(/\t/,$ttest_out);
	$pvalue = $ttest_out[@ttest_outTMP-1];

	my @homoLenNoPfam = split(/,/,$homoLenNoPfam);
	my @orphanLenNoPfam = split(/,/,$orphanLenNoPfam);
	foreach my $homoNoPfam(@homoLenNoPfam){ print DIVIDENOPFAMR "homologous\t$micros{$micros}\t$homoNoPfam\t$pvalue\n";}
	foreach my $orphanNoPfam(@orphanLenNoPfam){ print DIVIDENOPFAMR "orphan\t$micros{$micros}\t$orphanNoPfam\t$pvalue\n";}

	###################################################################################################

	############################ FOR PROTEINS that HAVE PFAM annotations #######################
	### get length of homologous and orphan proteins for this species
	my $homoFilePfam = abs_path($input)."/".$micros.".homologous.LEN_Pfam";
	chomp(my $homoLenPfam = `cut -f 2 $homoFilePfam`);
	$homoLenPfam =~ s/\n/,/g;
#	print $homoLenPfam;<>;

	my $orphanFileNewPfam = abs_path($input)."/".$micros.".orphan.LEN_newPfam";
	chomp(my $orphanLenNewPfam = `cut -f 2 $orphanFileNewPfam`);
	$orphanLenNewPfam =~ s/\n/,/g;
#	print $orphanLenNewPfam;<>;
	my $orphanFileHomoPfam = abs_path($input)."/".$micros.".orphan.LEN_homologPfam";
	chomp(my $orphanLenHomoPfam = `cut -f 2 $orphanFileHomoPfam`);
	$orphanLenHomoPfam =~ s/\n/,/g;
#	print $orphanLenHomoPfam;<>;
	my $orphanLenPfam = $orphanLenNewPfam.",".$orphanLenHomoPfam;

	### perform T-test
	### perform T-test
	$ttest_out = parse_utest($homoLenPfam,$orphanLenPfam,$dir);
	print DIVIDEPFAM $micros{$micros},"\t",$ttest_out,"\n";
	
	### create R file for plotting
	@ttest_out = split(/\t/,$ttest_out);
	$pvalue = $ttest_out[@ttest_outTMP-1];

	my @homoLenPfam = split(/,/,$homoLenPfam);
	my @orphanLenPfam = split(/,/,$orphanLenPfam);
	foreach my $homoPfam(@homoLenPfam){ print DIVIDEPFAMR "homologous\t$micros{$micros}\t$homoPfam\t$pvalue\n";}
	foreach my $orphanPfam(@orphanLenPfam){ print DIVIDEPFAMR "orphan\t$micros{$micros}\t$orphanPfam\t$pvalue\n";}
###################################################################################################
}
close (GENERAL);
close (GENERALR);
close (DIVIDENOPFAM);
close (DIVIDERNOPFAM);
close (DIVIDEPFAM);
close (DIVIDERPFAM);

#################################################### PLOTTING ####################################################

####### for all proteins (with and without Pfam annotations) #########
open(RSCRIPT,">$dir/lengthTest/homologous_vs_orphan.LEN.R");
print RSCRIPT "setwd(\"$dir/lengthTest\")
data <- read.table(\"$dir/lengthTest/homologous_vs_orphan.LEN.list\", sep=\'\\t\',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data\$taxon <- ordered(data\$taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))

png(\"$dir/lengthTest/homologous_vs_orphan.LEN.png\", width=800, height=800)

p = ggplot(data, aes(x=type, y=length,fill=type))+
    geom_boxplot()+
    facet_wrap(~taxon,ncol=6,nrow=2)+
    theme_set(theme_gray(base_size = 18)
)

p = p+labs(x=\"\", y=\"Length\")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15))
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPT);
system("Rscript $dir/lengthTest/homologous_vs_orphan.LEN.R");

####### for proteins that have Pfam annotations #########
open(RSCRIPT,">$dir/lengthTest/Pfam_homologous_vs_orphan.LEN.R");
print RSCRIPT "setwd(\"$dir/lengthTest\")
data <- read.table(\"$dir/lengthTest/Pfam_homologous_vs_orphan.LEN.list\", sep=\'\\t\',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data\$taxon <- ordered(data\$taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))

png(\"$dir/lengthTest/Pfam_homologous_vs_orphan.LEN.png\", width=800, height=800)

p = ggplot(data, aes(x=type, y=length,fill=type))+
    geom_boxplot()+
    facet_wrap(~taxon,ncol=6,nrow=2)+
    theme_set(theme_gray(base_size = 18)
)

p = p+labs(x=\"\", y=\"Length\")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15))
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPT);
system("Rscript $dir/lengthTest/Pfam_homologous_vs_orphan.LEN.R");

####### for proteins that DON'T Pfam annotations #########
open(RSCRIPT,">$dir/lengthTest/noPfam_homologous_vs_orphan.LEN.R");
print RSCRIPT "setwd(\"$dir/lengthTest\")
data <- read.table(\"$dir/lengthTest/noPfam_homologous_vs_orphan.LEN.list\", sep=\'\\t\',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data\$taxon <- ordered(data\$taxon, levels = c(\"E.hellem\",\"E.intestinallis\",\"E.cuniculi\",\"N.ceranae\",\"E.bieneusi\",\"V.corneae\",\"A.algerae\",\"A.locustae\",\"E.aedis\",\"V.culicis\",\"N.parisii\"))

png(\"$dir/lengthTest/noPfam_homologous_vs_orphan.LEN.list.png\", width=800, height=800)

p = ggplot(data, aes(x=type, y=length,fill=type))+
    geom_boxplot()+
    facet_wrap(~taxon,ncol=6,nrow=2)+
    theme_set(theme_gray(base_size = 18)
)

p = p+labs(x=\"\", y=\"Length\")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15))
p = p + scale_fill_brewer(palette=\"Set2\")
p
dev.off()
";
close (RSCRIPT);
system("Rscript $dir/lengthTest/noPfam_homologous_vs_orphan.LEN.R");

###################################################################################################################

print "FINISHED!! Check results in: $dir/lengthTest/\n";
exit;

sub parse_utest{
	my ($homoLen,$orphanLen,$dir) = @_;

	open(R,">$dir/utest.R") || die "Cannot create $dir/utest.R!!\n";
	print R "
		homologous = c($homoLen)
		orphan = c($orphanLen)
		wilcox.test(homologous,orphan)
		print(c(mean(homologous),mean(orphan)))
	";
	close (R);
	my $u_test = `Rscript $dir/utest.R`;

	$u_test =~ s/\s*$//;
#	print "\n\n#########################\n$u_test######\n";#<>;

	### get P-Value and mean lengths
	my $pvalue = 0;
	if($u_test =~ /p-value(.)+?\n/){
		chomp($pvalue = $&);
		$pvalue =~ s/p-value = //;
#		print "###$pvalue###";<>;
	} else {
		print "ERROR U-Test for $micros!!\n";<>;
	}

	my @u_test = split(/\n/,$u_test);
	my $mean = $u_test[@u_test-1]; $mean =~ s/\[1\]//; $mean =~ s/^\s*//; $mean =~ s/\s+/\t/;
#	print "#####",$mean,"####";<>;

	my $result = $mean."\t".$pvalue;
#	print $result;<>;
	return $result;
}



