#!/usr/bin/perl

use warnings;
use strict;
use List::MoreUtils 'distinct';
use List::Util 'sum';
use List::UtilsBy 'rev_nsort_by';
use Getopt::Long;

my $version = 2.1;
my $help=0;

my $usage = <<"USAGE";

perl compareZ_2.1.pl

Script for comparing two fasta files and calculating Z score.

Sorts output based on Z score.

Can output each sequence's percent of total reads in the table using the 
optional '-p' switch.
i.e. perl compareZ_2.1.pl -p file1 file2 >output

Version: $version

USAGE



my $percent_flag = 0;

GetOptions
( 
	"percent" => \$percent_flag,
	"help"       => \$help
)
  or die("Error in command line arguments\n");

if ( $help == 1 ) {
    print $usage;
    exit;
}



unless ( exists $ARGV[1] ) { print "\nNeed to specify 2 fasta files to compare. \n\n"; exit; }

my %lookup_file1;
my %lookup_file2;
my @all_seqs = ();

# file1
open( FASTA1, "<$ARGV[0]" ) or die("Failed to open file: $!\n");
my @seqID1   = ();
my @seqData1 = ();

while ( my $line1 = <FASTA1> ) {

    if ( $line1 =~ /^>.+/ ) { chomp $line1; push @seqID1, $line1; }
    else {
        chomp $line1;
        push @seqData1, $line1;
        push @all_seqs, $line1;
        if ( exists $lookup_file1{$line1} ) {
            $lookup_file1{$line1} += 1;
        }   
        else { $lookup_file1{$line1} = 1; } 
    }
}

# file2
open( FASTA2, "<$ARGV[1]" ) or die("Failed to open file: $!\n");
my @seqID2   = ();
my @seqData2 = ();

while ( my $line2 = <FASTA2> ) {

    if ( $line2 =~ /^>.+/ ) { chomp $line2; push @seqID1, $line2; }
    else {
        chomp $line2;
        push @seqData2, $line2;
        push @all_seqs, $line2;
        if ( exists $lookup_file2{$line2} ) { $lookup_file2{$line2} += 1; }
        else                                { $lookup_file2{$line2} = 1; }
    }
}

my $totalcount1 = sum values %lookup_file1;
my $totalcount2 = sum values %lookup_file2;

my @uniq_sequences = distinct(@all_seqs);
my @results        = ();
foreach (@uniq_sequences) {

    #chomp $_;
    my $count1     = 0;
    my $count2     = 0;
    my $zfunction1 = 0;
    my $zfunction2 = 0;

    my $norm1   = 0;
    my $norm2   = 0;
    my $zresult = 0;

    my $percent1 = 0;
    my $percent2 = 0;

    if ( exists $lookup_file1{$_} ) { $count1 = $lookup_file1{$_}; }
    if ( exists $lookup_file2{$_} ) { $count2 = $lookup_file2{$_}; }

    $norm1      = $count1 / $totalcount1;
    $norm2      = $count2 / $totalcount2;
    $zfunction1 = $norm1 * ( 1 - $norm1 ) / $totalcount1;
    $zfunction2 = $norm2 * ( 1 - $norm2 ) / $totalcount2;
    $zresult    = sprintf( "%.2f",
        ( $norm1 - $norm2 ) / sqrt( $zfunction1 + $zfunction2 ) );

    $percent1 = sprintf( "%.4f", ( $norm1 * 100 ) )
      ;    # NB due to rounding, total % in columns may !=100
    $percent2 = sprintf( "%.4f", ( $norm2 * 100 ) )
      ;    # set to 5f or 6f to get more accurate figures

    if ( $percent_flag == 0 ) {
        my $line = "$_\t$count1\t$count2\t$zresult";
# full record:
#my $line = "$_\t$count1\t$count2\t$totalcount1\t$totalcount2\t$norm1\t$norm2\t$zfunction1\t$zfunction2\t$zresult";             
     push @results, ($line);
    }
    else {
        my $line = "$_\t$count1\t$count2\t$zresult\t$percent1\t$percent2";      
        push @results, ($line);
    }



}

# header
if ( $percent_flag == 0 ) { print "Sequence\t$ARGV[0]\t$ARGV[1]\tZ-score\n"; }
else {
    print "Sequence\t$ARGV[0]\t$ARGV[1]\tZ-score\t$ARGV[0]_%\t$ARGV[1]_%\n";
}

# use if not using List::UtilsBy
#my @sorted =  sort { (split(' ', $b))[3] <=> (split(' ', $a))[3] } @results;

my @sorted = rev_nsort_by { ( split( ' ', $_ ) )[3] } @results;

foreach (@sorted) { print "$_\n" }
