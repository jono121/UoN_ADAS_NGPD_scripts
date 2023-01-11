#!/usr/bin/env perl

# pipeline.pl

$|=1; # no write buffer


use strict;
use warnings;
use Getopt::Long;
use Text::Match::FastAlternatives;
use PerlIO::gzip;
use Time::localtime;
use List::Util 'sum';

my $version     = "1.12";
my $help        = 0;
my $full_output = 0;

my $outfile      = "";    
my $fastq_file   = "";  
my $barcode_file = "";    # barcode file in format: BC01	xxxxxxx

my @bc_ID;
my @bc_seq;
my @fh;
my @fh_rc;
my @fh_tr1;
my @fh_tr2;
my @fh_tr3;

# reverse complement settings

my $rev_comp = 0;         #default no reverse complementing required
my $rev_seq;
my $rc_seq;

#translate
my $dna2pep_return;

my $dnaframe1;
my $dnaframe2;
my $dnaframe3;

my $frame1;
my $frame2;
my $frame3;

my $codon;
my @codons;

#iterate
my $leftmotiflist  = "";
my $rightmotiflist = "";    #files containing motifs (must be specified)
my @fh_it;                  # filehandle array for iteration
my @left_list;
my @right_list;
my $left_motif;
my $right_motif;
my $min_insert     = 1;    
my $nomatch_motifs = 0;   
my $match_motifs   = 0;
my @result_motifs;
my $rank             = 0;
my $fastaoutput      = 0;
my $translate_output = 0;
my %pep_ranks        = ();
my $matchseq         = "";
my $seq_to_match     = "";
my %CODON_TABLE;
my $f1_match = 0;
my $f2_match = 0;
my $f3_match = 0;

# this is the custom "opq stops" aa table
%CODON_TABLE = (
   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
   TAT => 'Y',TAC => 'Y',TAA => 'o',TAG => 'q',
   TGT => 'C',TGC => 'C',TGA => 'p',TGG => 'W',	   
   CTA => 'L',CTG => 'L',CTC => 'L',CTT => 'L',
   CCA => 'P',CCG => 'P',CCC => 'P',CCT => 'P',
   CAT => 'H',CAC => 'H',CAA => 'Q',CAG => 'Q',
   CGA => 'R',CGG => 'R',CGC => 'R',CGT => 'R',
   ATT => 'I',ATC => 'I',ATA => 'I',ATG => 'M',
   ACA => 'T',ACG => 'T',ACC => 'T',ACT => 'T',
   AAT => 'N',AAC => 'N',AAA => 'K',AAG => 'K',
   AGT => 'S',AGC => 'S',AGA => 'R',AGG => 'R',
   GTA => 'V',GTG => 'V',GTC => 'V',GTT => 'V',
   GCA => 'A',GCG => 'A',GCC => 'A',GCT => 'A',
   GAT => 'D',GAC => 'D',GAA => 'E',GAG => 'E',
   GGA => 'G',GGG => 'G',GGC => 'G',GGT => 'G');
   
GetOptions(
    "infile=s"   => \$fastq_file,
    "barcodes=s" => \$barcode_file,
    "revcomp"    => \$rev_comp,
    "fulloutput" => \$full_output,
    "rank"       => \$rank,
    "translate"  => \$translate_output,
    "outfile=s"  => \$outfile,
    "left=s"     => \$leftmotiflist,
    "right=s"    => \$rightmotiflist,
    "minimum=i"  => \$min_insert,
    'help'       => \$help,
    "fasta"      => \$fastaoutput
) or die("Error in command line arguments\n");

my $requires = <<"REQUIRES";

pipeline.pl:

Script to process an Ion Torrent fastq ( or .gz compressed fastq file) 
for the UoN Phage Display NGS analysis pipeline.
N.B. fastq files are assumed to be in single line format;
it will not work with multiline fastq files.

1. Demultiplexes using a text file containing expected barcodes
2. optionally reverse complements (using the '--revcomp' switch)
3. translates in all 3 reading frames (outputs the resulting amino acid 
	sequences to a file with the '--translate' switch)
4. tests each frame for the presence of a pair of flanking motifs with a
	minimum (default = 1_) insert between
5. extracts inserts and optionally ranks by frequency for each barcode 
	(using the '--rank' switch)

Output from each stage is written to the appropriate files derived from 
the outfile name and appropriate barcode ID.

Barcode file should be supplied in the following format:

BCnn	GATCGATCGATC

Lines starting with # in the barcode file are excluded.

------------------------------------------------------------------------  

REQUIRES

my $usage = <<"USAGE";

Usage: pipeline.pl --infile <fastq_file>  --outfile <name> --barcodes <txt_file> --left <txt_file> --right <txt_file> 

--infile <fastq_file>		Input file in fastq or gz format
--outfile <filename>		Base filename used for output
--barcodes <txt_file>		List of barcodes and their corresponding sequences, tab separated
--fasta				Write demultiplexed sequences to fasta files
--revcomp			Reverse complement sequences before translation
--translate			Write each frame of DNA->protein translation to files (frame1.fasta, frame2.fasta, frame3.fasta)
--left <filename>		text file containing list of vector motifs N-terminal to the variable region
--right <filename>		text file containing list of vector motifs C-terminal to the variable region
--minimum <integer>		minimum number of amino acids between motifs to identify as an insert (default = 1)
--rank				Rank identified inserts by frequency (and print normalised values)
--fulloutput			Equivalent to '--fasta --translate --rank'

Version: $version

USAGE

# checks

if ( $help == 1 ) {
    print $requires;
    print $usage;

    exit;
}

if ( $outfile eq "" ) {
    print
"\nNo output filename specified - default of 'NGS_analysis' will be used\n";
    $outfile = "NGS_analysis";
}

unless ( -e $fastq_file ) {

    print "\nInput file not found in .fastq or .gz format. \n";
    print $usage;
    exit;
}

unless ( -e $barcode_file ) {
    print $usage;
    print "\nBarcode file not found.\n";
    exit;
}

# barcode file
open( BARCODES, "<$barcode_file" ) or die;

while ( my $line = <BARCODES> ) {
    chomp $line;
    my @separated = split( '\s+', $line );
    if ( $separated[0] !~ /^#/ )    # '#' a BC to exclude
    {
        push @bc_ID,  $separated[0];
        push @bc_seq, $separated[1];
    }
}

my $nomatch = "$outfile" . "_nomatch" . ".fasta";

open (NOMATCH, ">$nomatch" ) or die;

my $bc_total_num = scalar(@bc_ID);    # number of barcodes used

if ( $fastaoutput == 1 || $full_output == 1 ) {
    my $outfasta;    # demultiplexed fastas
    for ( my $n = 0 ; $n < $bc_total_num ; $n++ ) {
        $outfasta = $outfile . "_" . $bc_ID[$n] . ".fasta";
        open $fh[$n], '>', "$outfasta";
    }
}

# reverse complement fastas
    if ( $rev_comp == 1 ) {
		my $out_rc; 
        for ( my $m = 0 ; $m < $bc_total_num ; $m++ ) {
            $out_rc = $outfile . "_" . $bc_ID[$m] . ".revcomp.fasta";
            open $fh_rc[$m], '>', "$out_rc";
        }
    }

if ( $translate_output == 1 || $full_output == 1 ) {

    my $out_trans1;    # translated f1
    for ( my $o = 0 ; $o < $bc_total_num ; $o++ ) {
        $out_trans1 = $outfile . "_" . $bc_ID[$o] . ".frame1.fasta";
        open $fh_tr1[$o], '>', "$out_trans1";
    }

    my $out_trans2;    # translated f2
    for ( my $o = 0 ; $o < $bc_total_num ; $o++ ) {
        $out_trans2 = $outfile . "_" . $bc_ID[$o] . ".frame2.fasta";
        open $fh_tr2[$o], '>', "$out_trans2";
    }

    my $out_trans3;    # translated f3
    for ( my $o = 0 ; $o < $bc_total_num ; $o++ ) {
        $out_trans3 = $outfile . "_" . $bc_ID[$o] . ".frame3.fasta";
        open $fh_tr3[$o], '>', "$out_trans3";
    }

}

my $out_iterate;    # inserts fastas
for ( my $p = 0 ; $p < $bc_total_num ; $p++ ) {
    $out_iterate = $outfile . "_" . $bc_ID[$p] . ".LR.fasta";
    open $fh_it[$p], '>', "$out_iterate";
}

# motif flanking region files
open( LEFTMOTIFLIST, "<$leftmotiflist" )
  or die "can't open N-terminal motif list: $!";
@left_list = <LEFTMOTIFLIST>;
chomp @left_list;
close LEFTMOTIFLIST;

open( RIGHTMOTIFLIST, "<$rightmotiflist" )
  or die "can't open C-terminal motif list: $!";
@right_list = <RIGHTMOTIFLIST>;
chomp @right_list;
close RIGHTMOTIFLIST;

my $count = 0;
my $fh_file;

if ( $fastq_file =~ /.+\.gz$/ ) {
    open( $fh_file, '<:gzip', $fastq_file )
      or die("Failed to open fastq.gz file: $!");

}

if ( $fastq_file =~ /.+\.fastq$|.+\.fq$/i ) {
    open( $fh_file, "<$fastq_file" ) or die("Failed to open fastq file: $!");
}

my $i   = 1;
my $ID  = "";
my $seq = "";

my @totals;
my $total_seqs = 0;
my $quick_test = Text::Match::FastAlternatives->new(@bc_seq);

print "\n"; #new line before starting screen counter

while ( my $line = <$fh_file> ) {

   chomp $line;

 if ($total_seqs % 1000 ==0 || eof == 1) {print "Sequences processed: $total_seqs\r";} # quick counter

    if ( $i % 4 == 1 ) {
        $total_seqs++;
        $ID = $line;
        $ID =~ tr /@/>/;    # assumes one @ only in ID
    }

    if ( $i % 4 == 2 ) {
        $seq = $line;
  

    if ( $quick_test->match($seq) ) {   # check whether any barcodes are present

      BARCODE:

        for ( my $j = 0 ; $j < $bc_total_num ; $j++ ) {

            if ( $seq =~ m/^.{0,5}$bc_seq[$j]/ip )
            {    # check which barcode matches, p preserve pre/match/postmatch

                $seq = ${^POSTMATCH};

                if ( $full_output == 1 || $fastaoutput == 1 ) {
                    print { $fh[$j] } "$ID\n$seq\n";
                }

                #reverse complement
                if ( $rev_comp == 1 ) {
                    my $dna = $seq;

                    $seq = revcomp($dna);

                    if ( $full_output == 1 ) {
                        print { $fh_rc[$j] } "$ID\n$seq\n";
                    }
                }

                #translate
                if ( length $seq >= 3 )
                {  # need at least one DNA triplet before attempting translation

                    $dnaframe1 = $seq;
                    $dnaframe2 = substr( $seq, 1 );
                    $dnaframe3 = substr( $seq, 2 );

                    $frame1 = translate($dnaframe1);
                    $frame2 = translate($dnaframe2);
                    $frame3 = translate($dnaframe3);

                    if ( $full_output == 1 || $translate_output == 1 ) {
                        print { $fh_tr1[$j] } "$ID\n$frame1\n";
                        print { $fh_tr2[$j] } "$ID\n$frame2\n";
                        print { $fh_tr3[$j] } "$ID\n$frame3\n";
                    }
                }

                # iterate motifs

              LEFT: foreach $left_motif (@left_list) {
                  RIGHT: foreach $right_motif (@right_list) {

                        $matchseq = "";    # make sure is empty
                        $seq_to_match =
                          qr/$left_motif(.{$min_insert,}?)$right_motif/
                          ;               

                        if ( $frame1 =~ $seq_to_match ) {
                            $matchseq = $1;
                            $match_motifs++;
                            $f1_match++;

                            print { $fh_it[$j] }
                              "$ID\n$matchseq\n"
                              ; # write matching sequence out with its ID to the LR file
                            my $bc = $bc_ID[$j];    #  get this current barcode

                            if ( $rank == 1 || $full_output == 1 ) {
                                if ( exists $pep_ranks{$bc}{$matchseq} ) {
                                    $pep_ranks{$bc}{$matchseq}++;
                                }
                                else {
                                    $pep_ranks{$bc}{$matchseq} = 1;
                                }
                            }
                            next RIGHT;
                        }

                        elsif ( $frame2 =~ $seq_to_match ) {
                            $matchseq = $1;
                            $match_motifs++;
                            $f2_match++;

                            print { $fh_it[$j] } "$ID\n$matchseq\n";
                            my $bc = $bc_ID[$j];

                            if ( $rank == 1 || $full_output == 1 ) {
                                if ( exists $pep_ranks{$bc}{$matchseq} ) {
                                    $pep_ranks{$bc}{$matchseq}++;
                                }
                                else {
                                    $pep_ranks{$bc}{$matchseq} = 1;
                                }
                            }
                            next RIGHT;
                        }

                        elsif ( $frame3 =~ $seq_to_match ) {
                            $matchseq = $1;
                            $match_motifs++;
                            $f3_match++;

                            print { $fh_it[$j] } "$ID\n$matchseq\n";
                            my $bc = $bc_ID[$j];

                            if ( $rank == 1 || $full_output == 1 ) {
                                if ( exists $pep_ranks{$bc}{$matchseq} ) {
                                    $pep_ranks{$bc}{$matchseq}++;
                                }
                                else {
                                    $pep_ranks{$bc}{$matchseq} = 1;
                                }
                            }
                            next RIGHT;
                        }

                        else {
                            $nomatch_motifs++;
                        }
                    }
                }

                # reset and count
                $seq    = "";
                $ID     = "";
                $rc_seq = "";
                $count++;
                $totals[$j]++;    # match was found in a BC, so record a hit
                last BARCODE;     # abort BC search if one just found
            }
        }
  }
      elsif ( $quick_test->match($seq) ==0) 
    {print NOMATCH "$ID\n$seq\n";}
  }

     $i++; # all done, increment counter for next line of fastq file
       
       
}

# rankings based on frequency
if ( $rank == 1 || $full_output == 1 ) {
    my $fh;

    foreach my $barcode ( sort keys %pep_ranks ) {
        my $sum_count = 0;
        my $out_rank  = $outfile . "_" . $barcode . ".rank.table";
        open $fh, '>', "$out_rank";

        print $fh "Sequence\tFrequency\tNorm_freq\n";
        $sum_count = sum values %{ $pep_ranks{$barcode} };

        foreach my $peptide (
            sort { $pep_ranks{$barcode}{$b} <=> $pep_ranks{$barcode}{$a} }
            keys %{ $pep_ranks{$barcode} } )
        {
            my $freq =
              $pep_ranks{$barcode}{$peptide};    # get frequency value from hash
            my $norm = sprintf( "%.4f", ( ( $freq / $sum_count ) * 100 ) );
            print $fh "$peptide\t$freq\t$norm\n";
        }
    }
}

#################################################
sub revcomp {

    $_[0] =~ tr/ACGT/TGCA/;
    $rc_seq = reverse( $_[0] );
    return $rc_seq;
}

#################################################
sub translate {
    $dna2pep_return = "";
    my @dna2pep = ();
    @codons = split( /(.{3})/, $_[0] );
    $codon  = "";
    foreach $codon (@codons) {
        if ( length $codon == 3 ) {
            if ( exists $CODON_TABLE{ uc $codon } ) {
                push @dna2pep, $CODON_TABLE{ uc $codon };
            }
            else { push @dna2pep, "u"; }
        }
    }
    $dna2pep_return = join '', @dna2pep;
    return $dna2pep_return;
}

###### LOG ######

my $t         = localtime;
my $timeprint = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    $t->year + 1900,
    $t->mon + 1,
    $t->mday, $t->hour, $t->min, $t->sec
);    
my $output_log = $outfile . "_" . $barcode_file . ".log";
my $bc_percent = sprintf( "%.2f", ( ( $count / $total_seqs ) * 100 ) )
  ;    # % of reads with a barcode
my $match_percent = "";

if ( $match_motifs == 0 || $count == 0 ) { $match_percent = "0.00"; }
else {
    $match_percent = sprintf( "%.2f", ( ( $match_motifs / $count ) * 100 ) );
}

my $nomatch_percent = "";
if ( $nomatch_motifs == 0 || $count == 0 ) { $nomatch_percent = "0.00"; }
else {
    $nomatch_percent =
      sprintf( "%.2f", ( ( $nomatch_motifs / $count ) * 100 ) );
}

open( LOGFILE, ">$output_log" ) or die("Failed to open log file: $!");
print LOGFILE "$timeprint\n";
print LOGFILE"\n";
print LOGFILE "Total sequences processed: $total_seqs\n";
print LOGFILE"\n";
print LOGFILE
  "Total barcodes found: $count ($bc_percent% of total sequences)\n";
print LOGFILE"\n";
print LOGFILE "Barcode sequence:\treads:\treads as % of all barcodes found:\n";

for ( my $t = 0 ; $t < $bc_total_num ; $t++ ) {
    if ( exists( $totals[$t] ) == 0 ) { $totals[$t] = 0; }
    my $each_bc_percent = "";
    if ( $totals[$t] == 0 || $count == 0 ) { $each_bc_percent = "0.00"; }
    else {
        $each_bc_percent =
          sprintf( "%.2f", ( ( $totals[$t] / $count ) * 100 ) );
    }
    print LOGFILE "$bc_ID[$t]\t$totals[$t]\t$each_bc_percent\n";

}
print LOGFILE"\n";
print LOGFILE"\n";
print LOGFILE
"Translated sequences containing motifs and an insert of at least $min_insert amino acid : $match_motifs ($match_percent% of barcoded sequences)\n";
print LOGFILE
"Translated sequences containing incomplete/absent paired motifs: $nomatch_motifs ($nomatch_percent% of barcoded sequences)\n";

my $frame1_stats = "";
my $frame2_stats = "";
my $frame3_stats = "";

if ( $f1_match == 0 || $match_motifs == 0 ) { $frame1_stats = "0.00"; }
else {
    $frame1_stats = sprintf( "%.2f", ( ( $f1_match / $match_motifs ) * 100 ) );
}

if ( $f2_match == 0 || $match_motifs == 0 ) { $frame2_stats = "0.00"; }
else {
    $frame2_stats = sprintf( "%.2f", ( ( $f2_match / $match_motifs ) * 100 ) );
}

if ( $f3_match == 0 || $match_motifs == 0 ) { $frame3_stats = "0.00"; }
else {
    $frame3_stats = sprintf( "%.2f", ( ( $f3_match / $match_motifs ) * 100 ) );
}

print LOGFILE"\n";
print LOGFILE
"frequency of motifs in frame 1: $f1_match ($frame1_stats % of total motifs found)\n";
print LOGFILE
"frequency of motifs in frame 2: $f2_match ($frame2_stats % of total motifs found)\n";
print LOGFILE
"frequency of motifs in frame 3: $f3_match ($frame3_stats % of total motifs found)\n";
print LOGFILE"\n";

close LOGFILE;
