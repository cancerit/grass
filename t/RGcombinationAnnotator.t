#!/usr/bin/perl

# for testing Grass::Annotation::RGcombinationAnnotator class

##########LICENCE##########
# PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
# Copyright (C) 2014 ICGC PanCancer Project
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not see:
#   http://www.gnu.org/licenses/gpl-2.0.html
##########LICENCE##########

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Data::Dumper;

use Grass::GenomeData;
use Grass::DataEntry;
use Grass::Annotation::RGendAnnotator;
use Grass::Annotation::RGcombinationAnnotator;

use Test::More 'no_plan';

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
my $entry = new Grass::DataEntry(-name       => $name,
				 -chr1       => $chr1,
				 -strand1    => $strand1,
				 -pos1_start => $pos1_start,
				 -pos1_end   => $pos1_end,
				 -chr2       => $chr2,
				 -strand2    => $strand2,
				 -pos2_start => $pos2_start,
				 -pos2_end   => $pos2_end,
				 -shard      => $shard,
				 -count      => $count );
my $within = 0;
my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $species = 'HUMAN';
#my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_58_1';
my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_74_1';

my $genome_data_ensembl = new Grass::GenomeData(-species     => $species,
						-ensembl_api => $ensembl_api,
						-gene_id_required => 0);

my $genome_data_cache = new Grass::GenomeData(-genome_cache => $genome_cache,
					      -gene_id_required => 0);



my $rgend1 = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
						   -end      => 1,
						   -genome_data => $genome_data_ensembl,
						   -within   => $within );
my $ann1 = $rgend1->annotate();

my $rgend2 = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
						   -end      => 2,
						   -genome_data => $genome_data_ensembl,
						   -within   => $within );
my $ann2 = $rgend2->annotate();

my $combi = new Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1->[0],
							  -anno2   => $ann2->[0],
							  -strand1 => $strand1,
							  -strand2 => $strand2,
							  -shard   => $shard);
$combi->combine(); 
my $rganno = $combi->anno();



my $rgend1_cache = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
							 -end      => 1,
							 -genome_data => $genome_data_cache,
							 -within   => $within );
my $ann1_cache = $rgend1_cache->annotate();

my $rgend2_cache = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
							 -end      => 2,
							 -genome_data => $genome_data_cache,
							 -within   => $within );
my $ann2_cache = $rgend2_cache->annotate();

my $combi_cache = new Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1_cache->[0],
								-anno2   => $ann2_cache->[0],
								-strand1 => $strand1,
								-strand2 => $strand2,
								-shard   => $shard);
$combi_cache->combine(); 
my $rganno_cache = $combi_cache->anno();

#my $res = get_result_58();
my $res = get_result_74();

ok($combi, 'object defined. ENSEMBL');
is (($combi->strand1()), $strand1 , "get strand1. ENSEMBL");
is (($combi->strand2()), $strand2 , "get strand2. ENSEMBL");
is (($combi->shard()), $shard , "get shard. ENSEMBL");

is ((Dumper($rganno)), $res , "get result. ENSEMBL");

ok($combi_cache, 'object defined. CACHE');
is (($combi_cache->strand1()), $strand1 , "get strand1. CACHE");
is (($combi_cache->strand2()), $strand2 , "get strand2. CACHE");
is (($combi_cache->shard()), $shard , "get shard. CACHE");

is ((Dumper($rganno_cache)), $res , "get result. CACHE");

#------------------------------------------------------------------------------------------------#
sub get_result_58 {
    my $res = <<'END';
$VAR1 = bless( {
                 'L5' => bless( {
                                  'gene_id' => 'ENSG00000172765',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'transcript' => 'ENST00000432054',
                                  'trans_length' => 5631,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => 1,
                                  'region_number' => 2,
                                  'strand' => '-1',
                                  'trans_region_count' => 4,
                                  'transcript_id' => 'ENST00000432054'
                                }, 'Grass::AnnoPoint' ),
                 'Ltype' => 'exon',
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
                                }, 'Grass::AnnoPoint' ),
                 'Htype' => '5UTRexon',
                 'id_fusion_flag' => '100',
                 'id_rg' => 'TEST',
                 'L3' => bless( {
                                  'gene_id' => 'ENSG00000172765',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'transcript' => 'ENST00000432054',
                                  'trans_length' => 5631,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => 1,
                                  'region_number' => 2,
                                  'strand' => '-1',
                                  'trans_region_count' => 4,
                                  'transcript_id' => 'ENST00000432054'
                                }, 'Grass::AnnoPoint' ),
                 'Llength' => 18901,
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
                                }, 'Grass::AnnoPoint' )
               }, 'Grass::Anno' );
END
    return($res);
}
#------------------------------------------------------------------------------------------------#
sub get_result_74 {
    my $res = <<'END';
$VAR1 = bless( {
                 'L5' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => 1,
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Grass::AnnoPoint' ),
                 'Ltype' => 'exon',
                 'H5' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'CG',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AA',
                                  'phase' => 1,
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Grass::AnnoPoint' ),
                 'Htype' => 'exon',
                 'id_fusion_flag' => 510,
                 'id_rg' => 'TEST',
                 'L3' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'TT',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AT',
                                  'phase' => 1,
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Grass::AnnoPoint' ),
                 'Llength' => 18901,
                 'Hlength' => 157118,
                 'H3' => bless( {
                                  'gene_id' => 'TMCC1',
                                  'up2' => 'CG',
                                  'region' => 'exon',
                                  'gene' => 'TMCC1',
                                  'transcript' => 'ENST00000393238',
                                  'trans_length' => 5992,
                                  'biotype' => 'protein_coding',
                                  'down2' => 'AA',
                                  'phase' => 1,
                                  'region_number' => 4,
                                  'strand' => '-1',
                                  'trans_region_count' => 6,
                                  'transcript_id' => 'ENST00000393238'
                                }, 'Grass::AnnoPoint' )
               }, 'Grass::Anno' );
END
    return($res);
}
#------------------------------------------------------------------------------------------------#
