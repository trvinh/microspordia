#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

sub usage {
    my $msg = shift;
    print "example: perl filterFasta.pl -i <folder contains fasta files> -m phyloprofileInput.long\n";
    print "-i\tFasta folder\n";
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
	my @f = <IN>;
	close (IN);

	my $newFa ="";
	for(my $i=0; $i < scalar(@f); $i+=2){
			chomp(my $line = $f[$i]);
			# print($line);<>;
			my $id = $line; $id =~ s/>//;
			# print $id;<>;
			if($geneID{$id}){
				# print $id;<>;
				# print "$line\n$f[$i+1]";<>;
				$newFa .= "$line\n$f[$i+1]";

			}
	}
	if(length($newFa) > 2){
		open(OUT,">$folder/new/$fileName") || die "cannot create $folder/new/$fileName";
		print OUT "$newFa";
		close (OUT);
	}


	print " finished!\n";
}


exit;
