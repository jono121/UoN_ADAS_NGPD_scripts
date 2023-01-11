#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $version = "1.1";

my $help    = 0;
my $infile  = "";
my $outfile = "default.fasta";

GetOptions(
    "infile=s"  => \$infile,
    "outfile=s" => \$outfile,
    'help'      => \$help
) or die("Error in command line arguments\n");

my $usage = <<"USAGE";

Usage: reverse_complement.pl --infile <inputfilename.fasta>  --outfile <outputfilename.fasta> 
	  
 --infile <inputfilename.fasta>		Input filename. File in single line FASTA format
 --outfile <outputfilename.fasta>	Output filename e.g. output.rc.fasta
 

USAGE

unless ( -e $infile ) {
    print "\nFasta file not found?\n";
    exit;
}

if ( $help == 1 ) {
    print $usage;
    exit;
}

open( FASTA, $infile ) or die("Failed to open FASTA file: $!\n");

my @seqID        = ();
my @seqData      = ();
my @reverse_comp = ();

while ( my $line = <FASTA> ) {
    chomp $line;
    if   ( $line =~ /^>+/ ) { push @seqID,   $line; }
    else                    { push @seqData, $line; }
}
close FASTA;

foreach (@seqData) {
    my $complement = $_;
    $complement =~
      tr/ACGTURYMKSWBDHVNacgturymkswbdhvn/TGCAAYRKMSWVHDBNtgcaayrkmswvhdbn/;
    my $reverse = reverse($complement);
    push( @reverse_comp, $reverse );
}

my $arr_size = scalar @reverse_comp;

open( my $fh, '>', $outfile );

for ( my $i = 0 ; $i < $arr_size ; $i++ ) {
    print $fh "$seqID[$i]\n$reverse_comp[$i]\n";
}
close $fh;
