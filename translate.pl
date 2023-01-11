#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;

my $version = "2.0";
my $help=0;

my $usage = <<"USAGE";

perl translate.pl -i INPUT.fasta -s opq -f 1 >OUTPUT.fasta

-infile		input fasta file to process
-frame		frame 1,2,3 = frame +1,+2,+3
-stop		stop codon treatment one of 'norm' 'amber' 'q' or 'opq' where 
		norm = codons TAA, TGA, TAG translated to '*', 
		amber = TAG stop translated as Q, 
		q = codons TAA, TGA, TAG translated to 'Q'
		opq = codons TAA, TGA, TAG, translated to o, p, and q respectively


Version: $version

USAGE

my ($infile, $stop, $frame);

GetOptions(
    "infile=s"   => \$infile, 
	"stop=s"     => \$stop,
	"frame=i"     => \$frame,
	"help"       => \$help
	);

if ( $help == 1 ) {
    print $usage;
    exit;
}
	
if( ! defined $infile) {
print "\nWARNING: Cannot proceed without input fasta file\n\n$usage\n\n"; exit;
}
if( ! defined $stop) {
print "\nWARNING: Cannot proceed without stop definition\n\n$usage\n\n"; exit;
}
if( ! defined $frame) {
print "\nWARNING: Cannot proceed without frame to translate\n\n$usage\n\n"; exit;
}


my %CODON_TABLE;

if ($stop eq "norm")
	{
%CODON_TABLE = (
   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
   TAT => 'Y',TAC => 'Y',TAA => '*',TAG => '*',
   TGT => 'C',TGC => 'C',TGA => '*',TGG => 'W',	   
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
}

elsif ($stop eq "q")
	{
	%CODON_TABLE = (
   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
   TAT => 'Y',TAC => 'Y',TAA => 'Q',TAG => 'Q',
   TGT => 'C',TGC => 'C',TGA => 'Q',TGG => 'W',	   
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
	}
	
elsif ($stop eq "amber")
	{
%CODON_TABLE = (
   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
   TAT => 'Y',TAC => 'Y',TAA => '*',TAG => 'Q',
   TGT => 'C',TGC => 'C',TGA => '*',TGG => 'W',	   
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
}

elsif ($stop eq "opq")
	{
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
}


my @seqID=();
my @seqData=();

my ($dnaframe1, $dnaframe2, $dnaframe3, $translation);

my $dna2pep_return="";
my @codons =();

open FASTA, "<$infile";

while ( my $line = <FASTA> ) {
    chomp $line;
    if   ( $line =~ /^>+/ ) { push @seqID,   $line; }
    else                    { push @seqData, $line; }
}
close FASTA;

my $arr_size = scalar @seqData;

for (my $i = 0; $i<$arr_size; $i++){
	
if ( length $seqData[$i] >= 3 ) # need at least one DNA triplet before attempting translation
{ 
$dnaframe1 = $seqData[$i];
$dnaframe2 = substr( $seqData[$i], 1 );
$dnaframe3 = substr( $seqData[$i], 2 );
	
	if ($frame == 1) {$translation = translate($dnaframe1);}
	elsif ($frame == 2){$translation = translate($dnaframe2);}
	elsif ($frame == 3){$translation = translate($dnaframe3);}
		
print "$seqID[$i]\n$translation\n";
}
}


sub translate {
    $dna2pep_return = "";
    my @dna2pep = ();
    @codons = split( /(.{3})/, $_[0] );
    my $codon  = "";
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
