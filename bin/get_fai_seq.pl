#!/usr/bin/perl -w

use Bio::DB::Sam; 

use strict;
use Getopt::Long;

my $ref = '/nfs/cancer_ref01/human/37/genome.fa';
my $coord = '';
my $help = 0;

GetOptions( 'ref:s' => \$ref,
            'coord:s'  => \$coord,
            'help:s'   => \$help);

# check inputs
if ($help) { usage(); }

# load an indexed fasta file
print "ref $ref\n";
my $fai = Bio::DB::Sam::Fai->load( $ref );

if ($coord =~ /^(\S+?):(\d+)$/) {
    my $chr = $1;
    my $pos = $2;
    $coord = $chr . ':' . $pos . '-' . $pos;
}

my $seq = $fai->fetch("$coord");
print "$seq\n";

#------------------------------------------------------------------------------------------------------------#
sub usage {

    print <<HERE;

get_fai_seq.pl

Description - gets a sequence from a fasta reference file (and associated fai index file) for supplied coordinates


options...

   -coord         : coordinate string format chr:start-end eg 1:123-456
   -ref           : fasta reference file (with associated fai file). For vcf out file generation.
   -help          : Print this message

examples...

get_fai_seq.pl -ref /nfs/cancer_ref01/human/37/genome.fa -coord 7:12350-12400 



Author : las


HERE

exit(0);
}
