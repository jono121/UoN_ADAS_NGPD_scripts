#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $version = "1.2";
my $help=0;

my $usage = <<"USAGE";

perl iterate_motifs.pl -infile INPUT.fasta -outfile OUTPUT -peptide 5 -left AEGEF.txt -right DPAKA.txt

Script to identify peptide sequences between known flanking motifs.

-infile		FASTA file to process
-outfile	result filename	
-peptide	minimum number of amino acids between motifs to accept
-left		file containing N-terminal flanking sequence
-right		file containing C-terminal flanking sequence


Version: $version

USAGE



my $infile  = "";
my $outfile = "";
my $left    = "";
my $right   = "";
my $peptide = 1;

GetOptions(
    "infile=s"  => \$infile,
    "outfile=s" => \$outfile, 
    "left=s"    => \$left,
    "right=s"   => \$right,
    "peptide=i" => \$peptide,
    "help"       => \$help
);


if ( $help == 1 ) {
    print $usage;
    exit;
}

my $m_logfile = $outfile . "Motifs_tested.log";
unlink $m_logfile; # delete

my $l_motif       = (0);
my $r_motif   = (0);
my $l_array_total = (0);
my $r_array_total = (0);

#left list
open( L_LIST, "<$left" ) or die "can't open left list: $!";
my @left_list = <L_LIST>;
chomp @left_list;
close L_LIST;

#right list
open( R_LIST, "<$right" ) or die "can't open right list: $!";
my @right_list = <R_LIST>;
chomp @right_list;
close R_LIST;

open LOG, ">$outfile.Motif.log";
print LOG
  "File\tSeq_Count\tL_pat\tR_pat\tL+R_Motif\tL_only\tR_only\t%_L+R\n";

my $sequenceL  = 0;
my $sequenceR  = 0;
my $count      = 0;
my $LR_hits    = 0;
my $L_hits     = 0;
my $R_hits     = 0;
my $PC_L_and_R = 0;


my @seqID   = ();
my @seqData = ();

open FASTA, "<$infile";

while ( my $line = <FASTA> ) {
    chomp $line;
    if   ( $line =~ /^>+/ ) { push @seqID,   $line; }
    else                    { push @seqData, $line; }
}
close FASTA;

my $arr_size = scalar @seqData;

open OUT1, ">$outfile.LR.fasta";
open OUT2, ">$outfile.L.fasta";
open OUT3, ">$outfile.R.fasta";


foreach $sequenceL (@left_list) {
    foreach $sequenceR (@right_list) {
        $l_motif     = $sequenceL;
        $r_motif = $sequenceR;

        for ( my $i = 0 ; $i < $arr_size ; $i++ ) {

            if ( $seqData[$i] =~ /$l_motif(.{$peptide,}?)$r_motif/ ) {
                print OUT1 "$seqID[$i]\n$1\n";    # insert only
                $LR_hits++;
            }
            elsif ( $seqData[$i] =~ /$l_motif(.{$peptide,})/ ) {

                #	print OUT2 "$seqID[$i]\n$1\n";
                print OUT2 "$seqID[$i]\n$seqData[$i]\n";    # whole sequence
                $L_hits++;
            }
            elsif ( $seqData[$i] =~ /(.{$peptide,})$r_motif/ ) {

                #	print OUT3 "$seqID[$i]\n$1\n";
                print OUT3 "$seqID[$i]\n$seqData[$i]\n";    #  whole sequence
                $R_hits++;
            }
           
        }
        $PC_L_and_R = sprintf( "%.2f", ( ( $LR_hits / $arr_size ) * 100 ) );
        open LOG, ">>$outfile.Motif.log";
        print LOG
"$outfile\t$arr_size\t$l_motif\t$r_motif\t$LR_hits\t$L_hits\t$R_hits\t$PC_L_and_R\n";

        # reset counters
        $LR_hits  = 0;
        $L_hits   = 0;
        $R_hits   = 0;
    
    }
}


close OUT1;
close OUT2;
close OUT3;
close LOG;
