#!/usr/bin/perl

# for testing Sanger::CGP::Grass::GenomeData::Transcript class

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::GenomeData::Exon;
use Sanger::CGP::Grass::GenomeData::Transcript;

use Test::More 'no_plan';

# existing entry

my $stable_id = 'ENST00000473358';
my $display_id = 'ENST00000473358';
my $length = 712;
my $strand = 1;
my $start = 29554;
my $end = 31097;
my $coding_region_start = 5;
my $coding_region_end = 11;
my $ccds_id = 'CCDS30547.1';
my $entrez_id;
my $gene_name = 'MIR1302-11';
my $gene_stable_id = 'MIR1302-11';
my $gene_biotype = 'protein_coding';
my $translation_length = 6;
my $accession = 'ENSP00000334393';
my $exons = [];


# make an exon new object
my $phase = 0;
my $end_phase = -1;
my $estart = 50;
my $eend = 100;
my $ecoding_region_start = 60;
my $ecoding_region_end = 99;
my $seq = 'acgtACGTNn';
my $Exon = new Sanger::CGP::Grass::GenomeData::Exon(-phase => $phase,
				       -end_phase => $end_phase,
				       -start => $estart,
				       -end => $eend,
				       -coding_region_start => $ecoding_region_start,
				       -coding_region_end => $ecoding_region_end,
				       -seq => $seq );
push $exons, $Exon;

# make a new transcript object
my $Transcript = new Sanger::CGP::Grass::GenomeData::Transcript(-stable_id => $stable_id,
						   -display_id => $display_id,
						   -length => $length,
						   -strand => $strand,
						   -start => $start,
						   -end => $end,
						   -coding_region_start => $coding_region_start,
						   -coding_region_end => $coding_region_end,
						   -ccds_id => $ccds_id,
						   -entrez_id => $entrez_id,
						   -gene_name => $gene_name,
						   -gene_stable_id => $gene_stable_id,
						   -gene_biotype => $gene_biotype,
						   -translation_length => $translation_length,
						   -accession => $accession,
				       -exons => $exons );

ok($Transcript, 'object defined');

is (($Transcript->stable_id()), $stable_id , "get stable_id");
is (($Transcript->display_id()), $display_id , "get display_id");
is (($Transcript->length()), $length , "get length");
is (($Transcript->strand()), $strand , "get strand");
is (($Transcript->start()), $start , "get start");
is (($Transcript->end()), $end , "get end");
is (($Transcript->coding_region_start()), $coding_region_start , "get coding_region_start");
is (($Transcript->coding_region_end()), $coding_region_end , "get coding_region_end");
is (($Transcript->ccds_id()), $ccds_id , "get ccds_id");
is (($Transcript->entrez_id()), $entrez_id , "get entrez_id");
is (($Transcript->gene_name()), $gene_name , "get gene_name");
is (($Transcript->gene_stable_id()), $gene_stable_id , "get gene_stable_id");
is (($Transcript->gene_biotype()), $gene_biotype , "get gene_biotype");
is (($Transcript->translation_length()), $translation_length , "get translation_length");
is (($Transcript->accession()), $accession , "get accession");
is_deeply (($Transcript->exons()), $exons , "get exons");
