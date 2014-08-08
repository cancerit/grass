#!/usr/bin/perl

# for testing Sanger::CGP::Grass::GenomeData::Exon class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::GenomeData::Exon;

use Test::More 'no_plan';

# existing entry


my $phase = 0;
my $end_phase = -1;
my $start = 50;
my $end = 100;
my $coding_region_start = 60;
my $coding_region_end = 99;
my $seq = 'acgtACGTNn';

# make a new object
my $Exon = new Sanger::CGP::Grass::GenomeData::Exon(-phase => $phase,
				       -end_phase => $end_phase,
				       -start => $start,
				       -end => $end,
				       -coding_region_start => $coding_region_start,
				       -coding_region_end => $coding_region_end,
				       -seq => $seq );

ok($Exon, 'object defined');

is (($Exon->phase()), $phase , "get phase");
is (($Exon->end_phase()), $end_phase , "get end_phase");
is (($Exon->start()), $start , "get start");
is (($Exon->end()), $end , "get end");
is (($Exon->coding_region_start()), $coding_region_start , "get coding_region_start");
is (($Exon->coding_region_end()), $coding_region_end , "get coding_region_end");
is (($Exon->seq()), $seq , "get seq");
