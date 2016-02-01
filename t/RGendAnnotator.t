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


# for testing Sanger::CGP::Grass::Annotation::RGendAnnotator class

use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use FindBin qw($Bin);
use Test::More 'no_plan';

my $genome_cache = "$Bin/../testData/vagrent.cache.gz";

for my $module (qw(Sanger::CGP::Grass::GenomeData Sanger::CGP::Grass::DataEntry Sanger::CGP::Grass::Annotation::RGendAnnotator)) {
  require_ok $module or BAIL_OUT "Can't load $module";
}

#/lustre/scratch104/sanger/am3/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz

my $ensembl_api= '';
if(exists $ENV{GRASS_ENS_API}) {
  $ensembl_api = $ENV{GRASS_ENS_API};
}

# make an entry object
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
my $end = 2;
my $within = 0;

# make an entry object where the gene is on the plus strand
my $name_2 = 'TEST_plus';
my $chr1_2 = '1';
my $strand1_2 = '+';
my $pos1_start_2 = 1552041;
my $pos1_end_2 = 1552041;
my $chr2_2 = '1';
my $strand2_2 = '+';
my $pos2_start_2 = 1558398;
my $pos2_end_2 = 1558410;
my $shard_2 = '';
my $count_2 = 4;
my $end_2 = 1;
my $within_2 = 100;


my $species = 'HUMAN';


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
         -count      => $count_2 ]);

my $res = get_result_74();
my $res_2 = get_result_74_2();

subtest 'Cache approach' => sub {
  my $genome_data_cache = new_ok('Sanger::CGP::Grass::GenomeData',
             [-genome_cache => $genome_cache,
              -gene_id_required => 0,
              -use_putative_coding => 0]);

  my $End2 = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
             [-entry   => $entry,
              -end     => $end,
              -within  => $within,
              -genome_data => $genome_data_cache]);
  my $anns2 = $End2->annotate();
  is (($End2->entry()), $entry , "get entry");
  is (($End2->end()), $end , "get end");
  is (($End2->within()), $within , "get within");
  is (($End2->ccds_only()), 1 , "get ccds_only");
  is_deeply ($anns2, $res , "get result cache");

  # make a new object - plus strand
  my $End2_2 = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
             [-entry   => $entry_2,
              -end     => $end_2,
              -within  => $within_2,
              -genome_data => $genome_data_cache]);
  my $anns2_2 = $End2_2->annotate();
  is (($End2_2->entry()), $entry_2 , "get entry (plus strand gene)");
  is (($End2_2->end()), $end_2 , "get end (plus strand gene)");
  is (($End2_2->within()), $within_2 , "get within (plus strand gene)");
  is (($End2_2->ccds_only()), 1 , "get ccds_only (plus strand gene)");
  is_deeply ($anns2_2, $res_2 , "get result cache (plus strand gene)");
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
              -gene_id_required => 0,
              -use_putative_coding => 0]);

    my $End = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
               [-entry   => $entry,
                -end     => $end,
                -within  => $within,
                -genome_data => $genome_data_ensembl]);
    my $anns = $End->annotate();
    is (($End->entry()), $entry , "get entry");
    is (($End->end()), $end , "get end");
    is (($End->within()), $within , "get within");
    is (($End->ccds_only()), 1 , "get ccds_only");
    is_deeply ($anns, $res , "get result ensembl");

    my $End_2 = new_ok('Sanger::CGP::Grass::Annotation::RGendAnnotator',
               [-entry   => $entry_2,
                -end     => $end_2,
                -within  => $within_2,
                -genome_data => $genome_data_ensembl]);
    my $anns_2 = $End_2->annotate();
    is (($End_2->entry()), $entry_2 , "get entry (plus strand gene)");
    is (($End_2->end()), $end_2 , "get end (plus strand gene)");
    is (($End_2->within()), $within_2 , "get within (plus strand gene)");
    is (($End_2->ccds_only()), 1 , "get ccds_only (plus strand gene)");
    is_deeply ($anns_2, $res_2 , "get result ensembl (plus strand gene)"); # has a putative coding transcript as well that is missing from the cache version
  };
};

#------------------------------------------------------------------------------------------------#
sub get_result_58 {
    my $res;
    eval <<'END';
$res = [
          bless( {
                   'id_rg' => 'TEST',
                   'H5' => bless( {
                                    'gene_id' => 'ENSG00000172765',
                                    'region' => '5UTRexon',
                                    'gene' => 'TMCC1',
                                    'transcript' => 'ENST00000432054',
                                    'trans_length' => 5631,
                                    'biotype' => 'protein_coding',
                                    'strand' => '-1',
                                    'region_number' => 2,
                                    'trans_region_count' => 4,
                                    'transcript_id' => 'ENST00000432054'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'Htype' => '5UTRexon',
                   'Hlength' => 17472,
                   'H3' => bless( {
                                    'gene_id' => 'ENSG00000172765',
                                    'region' => '5UTRexon',
                                    'gene' => 'TMCC1',
                                    'transcript' => 'ENST00000432054',
                                    'trans_length' => 5631,
                                    'biotype' => 'protein_coding',
                                    'strand' => '-1',
                                    'region_number' => 2,
                                    'trans_region_count' => 4,
                                    'transcript_id' => 'ENST00000432054'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' )
                 }, 'Sanger::CGP::Grass::Anno' ),
          bless( {
                   'id_rg' => 'TEST',
                   'H5' => bless( {
                                    'gene_id' => 'ENSG00000172765',
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
                   'Hlength' => 157118,
                   'H3' => bless( {
                                    'gene_id' => 'ENSG00000172765',
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
                 }, 'Sanger::CGP::Grass::Anno' ),
          bless( {
                   'id_rg' => 'TEST',
                   'H5' => bless( {
                                    'gene_id' => 'ENSG00000172765',
                                    'up2' => 'CG',
                                    'region' => 'exon',
                                    'gene' => 'TMCC1',
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
                   'Hlength' => 157118,
                   'H3' => bless( {
                                    'gene_id' => 'ENSG00000172765',
                                    'up2' => 'CG',
                                    'region' => 'exon',
                                    'gene' => 'TMCC1',
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
                 }, 'Sanger::CGP::Grass::Anno' ),
          bless( {
                   'id_rg' => 'TEST',
                   'H5' => bless( {
                                    'gene_id' => 'ENSG00000172765',
                                    'region' => '5UTRexon',
                                    'gene' => 'TMCC1',
                                    'transcript' => 'ENST00000307824',
                                    'trans_length' => 6152,
                                    'biotype' => 'protein_coding',
                                    'strand' => '-1',
                                    'region_number' => 5,
                                    'trans_region_count' => 7,
                                    'transcript_id' => 'ENST00000307824'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'Htype' => '5UTRexon',
                   'Hlength' => 222300,
                   'H3' => bless( {
                                    'gene_id' => 'ENSG00000172765',
                                    'region' => '5UTRexon',
                                    'gene' => 'TMCC1',
                                    'transcript' => 'ENST00000307824',
                                    'trans_length' => 6152,
                                    'biotype' => 'protein_coding',
                                    'strand' => '-1',
                                    'region_number' => 5,
                                    'trans_region_count' => 7,
                                    'transcript_id' => 'ENST00000307824'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' )
                 }, 'Sanger::CGP::Grass::Anno' )
        ];
END
    return($res);
}

#------------------------------------------------------------------------------------------------#
sub get_result_74 {
    my $res;
    eval <<'END';
$res = [
          bless( {
                   'id_rg' => 'TEST',
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
                 }, 'Sanger::CGP::Grass::Anno' )
        ];
END
    return($res);
}
#------------------------------------------------------------------------------------------------#
sub get_result_74_2 {
  my $res;
    eval <<'END';
$res = [
          bless( {
                   'L5' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1070,
                                    'transcript' => 'ENST00000505820',
                                    'trans_length' => 3305,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 20,
                                    'transcript_id' => 'ENST00000505820'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'Ltype' => 'intron',
                   'L3' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1070,
                                    'transcript' => 'ENST00000505820',
                                    'trans_length' => 3305,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 20,
                                    'transcript_id' => 'ENST00000505820'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'id_rg' => 'TEST_plus',
                   'Llength' => 1401
                 }, 'Sanger::CGP::Grass::Anno' ),
          bless( {
                   'L5' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1066,
                                    'transcript' => 'ENST00000520777',
                                    'trans_length' => 3321,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 20,
                                    'transcript_id' => 'ENST00000520777'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'Ltype' => 'intron',
                   'L3' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1066,
                                    'transcript' => 'ENST00000520777',
                                    'trans_length' => 3321,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 20,
                                    'transcript_id' => 'ENST00000520777'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'id_rg' => 'TEST_plus',
                   'Llength' => 1401
                 }, 'Sanger::CGP::Grass::Anno' ),
          bless( {
                   'L5' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1056,
                                    'transcript' => 'ENST00000355826',
                                    'trans_length' => 3290,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 20,
                                    'transcript_id' => 'ENST00000355826'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'Ltype' => 'intron',
                   'L3' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1056,
                                    'transcript' => 'ENST00000355826',
                                    'trans_length' => 3290,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 20,
                                    'transcript_id' => 'ENST00000355826'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'id_rg' => 'TEST_plus',
                   'Llength' => 1401
                 }, 'Sanger::CGP::Grass::Anno' ),
          bless( {
                   'L5' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1005,
                                    'transcript' => 'ENST00000518681',
                                    'trans_length' => 3116,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 19,
                                    'transcript_id' => 'ENST00000518681'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'Ltype' => 'intron',
                   'L3' => bless( {
                                    'gene_id' => 'MIB2',
                                    'region' => 'intron',
                                    'gene' => 'MIB2',
                                    'translation_length' => 1005,
                                    'transcript' => 'ENST00000518681',
                                    'trans_length' => 3116,
                                    'biotype' => 'protein_coding',
                                    'phase' => '2',
                                    'region_number' => 2,
                                    'strand' => '1',
                                    'trans_region_count' => 19,
                                    'transcript_id' => 'ENST00000518681'
                                  }, 'Sanger::CGP::Grass::AnnoPoint' ),
                   'id_rg' => 'TEST_plus',
                   'Llength' => 1401
                 }, 'Sanger::CGP::Grass::Anno' )
        ];
END
    return($res);
}

#------------------------------------------------------------------------------------------------#
