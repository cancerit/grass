#!/usr/bin/perl

# for testing Sanger::CGP::Grass::AnnoPoint class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::AnnoPoint;

use Test::More 'no_plan';

# existing entry
my $gene = 'TEST';
my $gene_id = 'TEST_ID';
my $biotype = 'coding';
my $transcript = 'TEST_TRANS';
my $transcript_id = 'TEST_TRANS_ID';
my $entrez_id = 'TEST_ENTREZ_ID';
my $strand = '-1';
my $phase = 0;
my $up2 = 'TG';
my $down2 = 'AA';
my $region = 'exon';
my $region_number = 3;
my $cds_pos = 57;
my $trans_pos = 186;
my $start_base = 3;
my $end_base = 7;
my $trans_region_count = 2;
my $trans_length = 1000;


# make a new object
my $Anno = new Sanger::CGP::Grass::AnnoPoint(-gene => $gene,
				-gene_id => $gene_id,
				-biotype => $biotype,
				-transcript => $transcript,
				-transcript_id => $transcript_id,
				-entrez_id => $entrez_id,
				-strand => $strand,
				-phase => $phase,
				-up2 => $up2,
				-down2 => $down2,
				-region => $region,
				-region_number => $region_number,
				-cds_pos => $cds_pos,
				-trans_pos => $trans_pos,
				-start_base => $start_base,
				-end_base => $end_base,
				-trans_region_count => $trans_region_count,
				-trans_length => $trans_length);

ok($Anno, 'object defined');

is (($Anno->gene()), $gene , "get gene");
is (($Anno->gene_id()), $gene_id , "get gene_id");
is (($Anno->biotype()), $biotype , "get biotype");
is (($Anno->transcript()), $transcript , "get transcript");
is (($Anno->transcript_id()), $transcript_id , "get transcript_id");
is (($Anno->entrez_id()), $entrez_id , "get entrez_id");
is (($Anno->strand()), $strand , "get strand");
is (($Anno->phase()), $phase , "get phase");
is (($Anno->up2()), $up2 , "get up2");
is (($Anno->down2()), $down2 , "get down2");
is (($Anno->region()), $region , "get region");
is (($Anno->region_number()), $region_number , "get region_number");
is (($Anno->cds_pos()), $cds_pos , "get cds_pos");
is (($Anno->trans_pos()), $trans_pos , "get trans_pos");
is (($Anno->start_base()), $start_base , "get start_base");
is (($Anno->end_base()), $end_base , "get end_base");
is (($Anno->trans_region_count()), $trans_region_count , "get trans_region_count");
is (($Anno->trans_length()), $trans_length , "get trans_length");
