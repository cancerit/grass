#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Lucy Stebbings <cgpit@sanger.ac.uk>
#
# This file is part of grass.
#
# grass is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########


# for testing Sanger::CGP::Grass::RGannotator class

use strict;
use warnings FATAL => 'all';
use FindBin qw($Bin);
use Test::More 'no_plan';
use Const::Fast qw(const);
use Data::Dumper;

# results should be...
const my $within => 0;
const my $res1 => "TMCC1\tENSG00000172765\tENST00000393238\t-1\t1\texon\t4\t6\t_\tTMCC1\tENSG00000172765\tENST00000393238\t-1\t1\texon\t4\t6\t_\t715\n";
const my $res2 => "TFPI\tENSG00000003436\tENST00000233156\t-1\t1\tintron\t2\t8\t_\tprotein_coding\tTFPI\tENSG00000003436\tENST00000233156\t-1\t_\t5UTRintron\t1\t8\t_\tprotein_coding\t200\tENSG00000224063,TFPI\n";
const my $res3 => "TFPI\tTFPI\tENST00000233156\t-1\t1\tintron\t2\t8\t_\tprotein_coding\tTFPI\tTFPI\tENST00000233156\t-1\t_\t5UTRintron\t1\t8\t_\tprotein_coding\t200\tTFPI\n";

for my $module (qw(Sanger::CGP::Grass::GenomeData Sanger::CGP::Grass::DataEntry Sanger::CGP::Grass::Annotation::RGannotator)) {
  require_ok $module or BAIL_OUT "Can't load $module";
}

my $species = 'HUMAN';
#my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_74_1';
my $ensembl_api= '';

if(exists $ENV{GRASS_ENS_API}) {
  $ensembl_api = $ENV{GRASS_ENS_API};
}

#my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $genome_cache = "$Bin/../testData/vagrent.cache.gz";

# make a dataset
my $name = 'TEST';
my $chr1 = '3';
my $strand1 = '-';
my $pos1_start = 129389225;
my $pos1_end = 129389225;
my $chr2 = '3';
my $strand2 = '-';
my $pos2_start = 129390103;
my $pos2_end = 129390103;
my $shard = 'AA';
my $count = 6;
my $entry = new_ok('Sanger::CGP::Grass::DataEntry',
        [-name       => $name,
				 -chr1       => $chr1,
				 -strand1    => $strand1,
				 -pos1_start => $pos1_start,
				 -pos1_end   => $pos1_end,
				 -chr2       => $chr2,
				 -strand2    => $strand2,
				 -pos2_start => $pos2_start,
				 -pos2_end   => $pos2_end,
				 -shard      => $shard,
				 -count      => $count]);
my $dataset1 = [];
push @$dataset1, $entry;

my $name_2 = 37464190;
my $chr1_2 = '2';
my $strand1_2 = '+';
my $pos1_start_2 = 188365485;
my $pos1_end_2 = 188365837;
my $chr2_2 = '2';
my $strand2_2 = '+';
my $pos2_start_2 = 188417763;
my $pos2_end_2 = 188418155;
my $shard_2 = '';
my $count_2 = 4;
my $entry_2 = new_ok('Sanger::CGP::Grass::DataEntry',
        [-name       => $name_2,
				 -chr1       => $chr1_2,
				 -strand1    => $strand1_2,
				 -pos1_start => $pos1_start_2,
				 -pos1_end   => $pos1_end_2,
				 -chr2       => $chr2_2,
				 -strand2    => $strand2_2,
				 -pos2_start => $pos2_start_2,
				 -pos2_end   => $pos2_end_2,
				 -shard      => $shard_2,
				 -count      => $count_2]);
my $dataset2 = [];
push @$dataset2, $entry_2;

SKIP: {
  unless(-d $ensembl_api) {
    my $message = "SKIPPING: Root of Ensembl API not found at '$ensembl_api', set ENV: GRASS_ENS_API to enable these tests";
    warn "$message\n";
    skip $message, '1';
  }
  subtest 'Ensembl API approach' => sub {
    # set up access to genome data  from EnsemblDB (if species/ensembl_api supplied) or a cached flat file version (if genome_cache supplied)
    my $genome_data_ensembl = new_ok('Sanger::CGP::Grass::GenomeData',
               [-species     => $species,
                -ensembl_api => $ensembl_api,
                -gene_id_required => 1]);

    # make a new object
    my $rgann1 = new_ok('Sanger::CGP::Grass::Annotation::RGannotator',
               [-dataset         => $dataset1,
                -within          => $within,
                -genome_data     => $genome_data_ensembl,
                -list_between    => 0,
                -multi_annos     => 1,
                -show_biotype    => 0]);
    ok(!$rgann1->getAnnotation(), 'getAnnotation Ensembl should fail');

    is (($rgann1->within()), $within , "get within. ENSEMBL");
    is (($rgann1->list_between()), 0 , "get list_between. ENSEMBL");
    is (($rgann1->show_biotype()), 0 , "get show_biotype. ENSEMBL");
    is (($rgann1->format_for_printing()), $res1 , "get result. ENSEMBL");

    # make a new object
    my $rgann2 = new_ok('Sanger::CGP::Grass::Annotation::RGannotator',
               [-dataset         => $dataset2,
                -within          => $within,
                -genome_data     => $genome_data_ensembl,
                -list_between    => 1,
                -show_biotype    => 1]);
    ok($rgann2->getAnnotation(), '2getAnnotation Ensembl');

    is (($rgann2->within()), $within , "2get within. ENSEMBL");
    is (($rgann2->list_between()), 1 , "2get list_between. ENSEMBL");
    is (($rgann2->show_biotype()), 1 , "2get show_biotype. ENSEMBL");
    is (($rgann2->format_for_printing()), $res2 , "2get result. ENSEMBL");
  };
};

subtest 'Cache approach' => sub {
  my $genome_data_cache = new_ok('Sanger::CGP::Grass::GenomeData',
                 [-genome_cache => $genome_cache,
                  -gene_id_required => 0]);
  my $rgann3 = new_ok('Sanger::CGP::Grass::Annotation::RGannotator',
             [-dataset         => $dataset2,
              -within          => $within,
              -genome_data    => $genome_data_cache,
              -list_between    => 1,
              -show_biotype    => 1]);
  ok($rgann3->getAnnotation(), '3getAnnotation Cache');

  is (($rgann3->within()), $within , "3get within. Cache");
  is (($rgann3->list_between()), 1 , "3get list_between. Cache");
  is (($rgann3->show_biotype()), 1 , "3get show_biotype. Cache");
  is (($rgann3->format_for_printing()), $res3 , "3get result. Cache");
};

#------------------------------------------------------------------------------------------------#
