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


# for testing Sanger::CGP::Grass::Annotation::RGcombinationAnnotator class


use strict;
use warnings FATAL => 'all';
use FindBin qw($Bin);
use Data::Dumper;

use Test::More 'no_plan';

for my $module (qw(Sanger::CGP::Grass::GenomeData Sanger::CGP::Grass::DataEntry Sanger::CGP::Grass::Annotation::RGendAnnotator Sanger::CGP::Grass::Annotation::RGcombinationAnnotator)) {
  require_ok $module or BAIL_OUT "Can't load $module";
}

#my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $genome_cache = "$Bin/../testData/vagrent.cache.gz";

my $species = 'HUMAN';
#my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_58_1';
#my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_74_1';
my $ensembl_api= '';

if(exists $ENV{GRASS_ENS_API}) {
  $ensembl_api = $ENV{GRASS_ENS_API};
}

my $res = get_result_58();
my $res_74 = get_result_74();

my $name = 'TEST';
my $chr1 = '3';
my $strand1 = '-';
my $pos1_start = 129389225;
my $pos1_end = 129389225;
my $chr2 = '3';
my $strand2 = '+';
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
				 -count      => $count ]);
my $within = 0;

subtest 'Cache approach' => sub {
  my $genome_data_cache = new_ok('Sanger::CGP::Grass::GenomeData',
                 [-genome_cache => $genome_cache,
                  -gene_id_required => 0]);

  my $rgend1_cache = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
                [-entry    => $entry,
                 -end      => 1,
                 -genome_data => $genome_data_cache,
                 -within   => $within]);
  my $ann1_cache = $rgend1_cache->annotate();

  my $rgend2_cache = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
                [-entry    => $entry,
                 -end      => 2,
                 -genome_data => $genome_data_cache,
                 -within   => $within]);
  my $ann2_cache = $rgend2_cache->annotate();

  my $combi_cache = new_ok('Sanger::CGP::Grass::Annotation::RGcombinationAnnotator',
                 [-anno1   => $ann1_cache->[0],
                  -anno2   => $ann2_cache->[0],
                  -strand1 => $strand1,
                  -strand2 => $strand2,
                  -shard   => $shard]);
  ok($combi_cache->combine(), 'Successful combine - Cache');
  my $rganno_cache = $combi_cache->anno();
  is (($combi_cache->strand1()), $strand1 , "get strand1. Cache");
  is (($combi_cache->strand2()), $strand2 , "get strand2. Cache");
  is (($combi_cache->shard()), $shard , "get shard. Cache");
  is_deeply ($rganno_cache, $res_74 , "get result. Cache");
};

SKIP: {
  unless(-d $ensembl_api) {
    my $message = "SKIPPING: Root of Ensembl API not found at '$ensembl_api', set ENV: GRASS_ENS_API to enable these tests";
    warn "$message\n";
    skip $message, '1';
  }
  subtest 'Ensembl API approach' => sub {
    my $genome_data_ensembl = new_ok('Sanger::CGP::Grass::GenomeData',
               [-species     => $species,
                -ensembl_api => $ensembl_api,
                -gene_id_required => 0]);
    my $rgend1 = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
                  [-entry    => $entry,
                   -end      => 1,
                   -genome_data => $genome_data_ensembl,
                   -within   => $within]);
    my $ann1 = $rgend1->annotate();

    my $rgend2 = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
                  [-entry    => $entry,
                   -end      => 2,
                   -genome_data => $genome_data_ensembl,
                   -within   => $within]);
    my $ann2 = $rgend2->annotate();

    my $combi = new_ok('Sanger::CGP::Grass::Annotation::RGcombinationAnnotator',
                   [-anno1   => $ann1->[0],
                    -anno2   => $ann2->[0],
                    -strand1 => $strand1,
                    -strand2 => $strand2,
                    -shard   => $shard]);
    ok($combi->combine(), 'Successful combine - Ensembl');
    my $rganno = $combi->anno();

    is (($combi->strand1()), $strand1 , "get strand1. Ensembl");
    is (($combi->strand2()), $strand2 , "get strand2. Ensembl");
    is (($combi->shard()), $shard , "get shard. Ensembl");
    is_deeply ($rganno, $res_74 , "get result. Ensembl");
  };
};

#------------------------------------------------------------------------------------------------#
sub get_result_58 {
  my $res;
    eval <<'END';
$res = bless( {
                 'L5' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000329333',
                                  'trans_length' => 6152,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => '1',
                                  'region_number' => 5,
                                  'strand' => '-1',
                                  'trans_region_count' => 7,
                                  'transcript_id' => 'ENST00000329333'
                                }, 'Sanger::CGP::Grass::AnnoPoint' ),
                 'Ltype' => 'exon',
                 'H5' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'CG',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000329333',
                                  'trans_length' => 6152,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AA',
                                  'phase' => '1',
                                  'region_number' => 5,
                                  'strand' => '-1',
                                  'trans_region_count' => 7,
                                  'transcript_id' => 'ENST00000329333'
                                }, 'Sanger::CGP::Grass::AnnoPoint' ),
                 'Htype' => 'exon',
                 'id_fusion_flag' => 510,
                 'id_rg' => 'TEST',
                 'L3' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000329333',
                                  'trans_length' => 6152,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => '1',
                                  'region_number' => 5,
                                  'strand' => '-1',
                                  'trans_region_count' => 7,
                                  'transcript_id' => 'ENST00000329333'
                                }, 'Sanger::CGP::Grass::AnnoPoint' ),
                 'Llength' => 18901,
                 'Hlength' => 157118,
                 'H3' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'CG',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000329333',
                                  'trans_length' => 6152,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AA',
                                  'phase' => '1',
                                  'region_number' => 5,
                                  'strand' => '-1',
                                  'trans_region_count' => 7,
                                  'transcript_id' => 'ENST00000329333'
                                }, 'Sanger::CGP::Grass::AnnoPoint' )
               }, 'Sanger::CGP::Grass::Anno' );
END
    return($res);
}
#------------------------------------------------------------------------------------------------#
sub get_result_74 {
  my $res;
    eval <<'END';
$res = bless( {
                 'L5' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => '1',
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Sanger::CGP::Grass::AnnoPoint' ),
                 'Ltype' => 'exon',
                 'H5' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'CG',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AA',
                                  'phase' => '1',
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Sanger::CGP::Grass::AnnoPoint' ),
                 'Htype' => 'exon',
                 'id_fusion_flag' => 510,
                 'id_rg' => 'TEST',
                 'L3' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => '1',
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Sanger::CGP::Grass::AnnoPoint' ),
                 'Llength' => 18901,
                 'Hlength' => 157118,
                 'H3' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'CG',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'translation_length' => 653,
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AA',
                                  'phase' => '1',
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Sanger::CGP::Grass::AnnoPoint' )
               }, 'Sanger::CGP::Grass::Anno' );
END
    return($res);
}
#------------------------------------------------------------------------------------------------#
