A selection of Perl scripts to process Ion Torrent DNA sequencing data to identify peptides enriched through panning of phage display libraries.

demultiplex.sh can be used to identify and remove the 96 Ion Torrent
DNA barcodes (also listed in barcodes.txt) from a DNA FASTA file.

The scripts reverse_complement.pl, translate.pl, iterate_motifs.pl, and compareZ_2.1.pl
are designed to be run on demultiplexed DNA fasta files. Instructions on 
using each file can be obtained by typing "perl translate.pl -help" at
the terminal for each of the .pl files.

The file pipeline.pl is designed to process an Ion Torrent fastq file, or a .gz compressed fastq file, 
and will perform the demultiplexing, reverse complementing, translating and motif iteration.
Resultant .LR.fasta files can then be compared using compareZ_2.1.pl to calculate Z scores.

Jan 2022.

