#!/usr/bin/perl
use forks;
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Std;
use Array::Utils qw(:all);
use Parallel::ForkManager;

=description
do pfam annotation for list of proteins
Version 1.0 (29.02.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl getPfam.pl -i ortholog.list\n";
    print "-i\tOrtholog list\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No ortholog list given\n");
open(IN,$input) || die "Cannot find $input!!\n";
my @in = <IN>;
close (IN);

### paths
#my $microsAnno = "/home/vinh/Desktop/data/project/FACT_annotation/microsAnnotation";	# folder contains annotation files
my $factDir = "/home/vinh/Desktop/programs/oneseq_distri/bin/fas";

### get all PFAMs of homologous proteins
my $homologPFAM = "/home/vinh/Desktop/data/project/iniAnalysis/homolog_orphan_dir/allHomologous.PFAM";
open(HOMOPFAM,$homologPFAM) || die "Cannot open $homologPFAM!!\n";
my @homologPFAM = <HOMOPFAM>;
close (HOMOPFAM);
my %homologPFAM;
foreach my $homoPfam (@homologPFAM){
	chomp($homoPfam);
	$homologPFAM{$homoPfam} = 1;
}
#print scalar(keys %homologPFAM);<>;

### working dir
my $dir = abs_path($input);
my @dir = split(/\//,$dir); pop(@dir);
$dir = join("/",@dir);
#print $dir;<>;
unless(-d "$dir/tmp"){
	mkdir("$dir/tmp");
}

# output file
open(PFAM,">$input.PFAM") || die "Cannot create $input.PFAM!!\n";
open(LEN,">$input.LEN") || die "Cannot create $input.LEN!!\n";

my $c = 0;
my %pfam; # $pfam{$protID} = pfam_Accession
foreach my $protID (@in){
	if(length($protID)>2){
		chomp($protID);	#  annca_339:H312_02262
		my @protIDTMP = split(/:/,$protID);
		my $specID = $protIDTMP[0];

		### get fasta seq for this protein
		open(FASTA,">$dir/tmp/$protID.fa") || die "Cannot create $dir/tmp/$protID.fa!!\n";
		my $fasta = getFasta($protID);
		print FASTA ">",$fasta;
		close (FASTA);

		### get sequence length
		my $len = getLength($fasta);
#		print $len;<>;
print "HIER\n";<>;
		### do Hmmscan agains pfam_A.hmm
		system("hmmscan -E 0.001 --noali --tblout $dir/tmp/$protID.fa.out $factDir/Pfam/Pfam-hmms/Pfam-A.hmm $dir/tmp/$protID.fa");
		open(HMM,"$dir/tmp/$protID.fa.out") || die "Cannot open $dir/tmp/$protID.fa.out!!\n";
		my @hmm = <HMM>;
		close (HMM);
print "check $dir/tmp/$protID.fa.out\n";<>;
		### parse hmm output to get pfam accession numbers
		my $hit = 0; # $hit=0 if this protein has no pfam annotation
		my $tag = "newPFAM";	# $tag = noPFAM, newPFAM or homologPFAM
		foreach my $hmmLine(@hmm){
			unless($hmmLine =~ /^#/){
#				print $hmmLine,"\n";
				my $hmmLineTMP = $hmmLine;
				$hmmLineTMP =~ s/\s+/\t/g;
				my @hmmLineTMP = split(/\t/,$hmmLineTMP);
#				print "HRERE: ",$hmmLineTMP[1],"\n";
				unless($pfam{$hmmLineTMP[2]}){
					$pfam{$protID} = $hmmLineTMP[1];
				} else {
					$pfam{$protID} .= ";".$hmmLineTMP[1];
				}

				if($homologPFAM{$hmmLineTMP[1]}){
					$tag = "homologPFAM";
				}
				$hit = 1;
			}
		}

		## remove fasta and hmmscan output
		system("rm $dir/tmp/$protID.fa");
		system("rm $dir/tmp/$protID.fa.out");

		### print output
#		print "$protID\t$tag\t$pfam{$protID}\n";<>;
		if($hit==0){
			print PFAM "$protID\tnoPFAM\tnoPFAM\n";
		} else {
			print PFAM "$protID\t$tag\t$pfam{$protID}\n";
		}

		print LEN "$protID\t$len\n";
#		<>;
	}
}
close (PFAM);
close (LEN);

print "FINISHED. Check outputs in\n\t$input.PFAM\n\t$input.LEN\n";
exit;


### get sequence from multifasta file for a protein (species_name:proteinID)
sub getFasta{
	my $desc = $_[0];	# $desc = SPECIES:PROTEIN_ID
	#print $desc;<>;
	# get species name
	my @tmp = split(/\:/,$desc);
	my $species = $tmp[0];  ## $tmp[0] = species name, $tmp[1] = protein_id
	# open multifasta file for this species
	my $faFolder = "/share/project/vinh/dataset/allSpeciesFasta";
	my $file = "";
	if(-e "$faFolder/$species.fasta"){
		$file = "$faFolder/$species.fasta";
	} elsif(-e "$faFolder/$species.fa") {
		$file = "$faFolder/$species.fa";
	} else {
		die "No fasta file found for $species!!\n";
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
	#print "FASTA=",$result;<>;
	close (FILE);
	return $result;
}

sub getLength{
	my $fasta = $_[0];	# input is a fasta sequence
	my @tmp = split(/\n/,$fasta);
	my $len = length($tmp[1]);
	return $len;
}
