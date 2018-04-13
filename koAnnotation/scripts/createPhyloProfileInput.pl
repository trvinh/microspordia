#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
create PhyloProfile input from annotation result
2018.03.12
=cut

sub usage {
    my $msg = shift;
    print "example: perl createPhyloProfileInput.pl -i lca.list.koAnnotation.KO -t 4koAnnotataion_manually_taxonomy.list\n";
    print "-i\tinput annotated file\n";
    print "-t\tfile contains taxonomy IDs for all species\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_t,$opt_d);
getopts('i:t:d:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $tax = ($opt_t) ? $opt_t : usage("ERROR: No taxonomy file given\n");

### get taxonomy info
open(TAX,$tax) || die "Cannot open $tax!\n";
my %tax;  # $tax{specName} = taxonID
foreach my $line(<TAX>){
  chomp($line);
  my @tmp = split(/\t/,$line);
  $tax{$tmp[0]} = $tmp[2];
}

### read input
open(IN,$input);
my @in = <IN>;
close (IN);
my $in = join("",@in);
my @groups = split(/### /,$in);

open(OUT,">$input.phyloprofile");
print OUT "geneID\tncbiID\torthoID\tFAS\tpatristicDist\n";

foreach my $group(@groups){
  if(length($group) > 2){
    my @lines = split(/\n/,$group);
    my $groupID = shift(@lines); my @groupIDtmp = split(/\t/,$groupID); # group ID = $groupIDtmp[0]
    foreach my $line(@lines){
      if(length($line) > 2){
        chomp($line);
        my @tmp = split(/\t/,$line);
        my @specTMP = split(/:/,$tmp[0]);
        unless ($tax{$specTMP[0]}) {
          print "NO TAXONOMY ID FOR $specTMP[0]\n";<>;
        }

        print OUT $groupIDtmp[0],"\tncbi",$tax{$specTMP[0]},"\t",$tmp[0],"\t",$tmp[2],"\t",$tmp[3],"\n";
      }
    }
  }
}

print "FINISHED! Check output $input.phyloprofile\n";
exit;
