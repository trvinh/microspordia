#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

sub usage {
    my $msg = shift;
    print "example: perl filterDomainFile.pl -i <folder contains domain files> -m phyloprofileInput.long\n";
    print "-i\tDomain folder\n";
		print "-m\tPhyloProfile input file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_m);
getopts('i:m:');

my $folder = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $matrixIn = ($opt_m) ? $opt_m : usage("ERROR: No phyloprofile input file given\n");

open(MATRIX,$matrixIn) || die "Cannot open $matrixIn!!\n";
my %geneID;
foreach my $line(<MATRIX>){
	chomp ($line);
	# print $line;<>;
	my @tmp = split(/\t/,$line);
	if($tmp[2] ne "NA"){
		# print $line;<>;
		$geneID{$tmp[2]} = 1;
	}
}

my @files = glob("$folder/*.*");
mkdir("$folder/new");

foreach my $fileIn(@files){
	print "$fileIn...";
	### get file name
	my @fileName = split(/\//,$fileIn);
	my $fileName = pop(@fileName);

	open(IN,"$fileIn") || die "cannot open $fileIn!\n";
	open(OUT,">$folder/new/$fileName") || die "cannot create $folder/new/$fileName";

	foreach my $line(<IN>){
			chomp($line);
			# print($line);<>;
			my @tmp = split(/\t/,$line);
			my @id = split(/#/,$tmp[0]);

			if($geneID{$id[1]}){
				print OUT $line,"\n";
			}
	}
	close (IN);
	close (OUT);

	print " finished!\n";
}


exit;
