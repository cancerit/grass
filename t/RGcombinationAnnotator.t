#!/usr/bin/perl

# for testing Sanger::CGP::Grass::Annotation::RGcombinationAnnotator class

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Data::Dumper;

use Sanger::CGP::Grass::GenomeData;
use Sanger::CGP::Grass::DataEntry;
use Sanger::CGP::Grass::Annotation::RGendAnnotator;
use Sanger::CGP::Grass::Annotation::RGcombinationAnnotator;

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
my $entry = new Sanger::CGP::Grass::DataEntry(-name       => $name,
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

my $genome_data_ensembl = new Sanger::CGP::Grass::GenomeData(-species     => $species,
						-ensembl_api => $ensembl_api,
						-gene_id_required => 0);

my $genome_data_cache = new Sanger::CGP::Grass::GenomeData(-genome_cache => $genome_cache,
					      -gene_id_required => 0);



my $rgend1 = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry    => $entry,
						   -end      => 1,
						   -genome_data => $genome_data_ensembl,
						   -within   => $within );
my $ann1 = $rgend1->annotate();

my $rgend2 = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry    => $entry,
						   -end      => 2,
						   -genome_data => $genome_data_ensembl,
						   -within   => $within );
my $ann2 = $rgend2->annotate();

my $combi = new Sanger::CGP::Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1->[0],
							  -anno2   => $ann2->[0],
							  -strand1 => $strand1,
							  -strand2 => $strand2,
							  -shard   => $shard);
$combi->combine(); 
my $rganno = $combi->anno();



my $rgend1_cache = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry    => $entry,
							 -end      => 1,
							 -genome_data => $genome_data_cache,
							 -within   => $within );
my $ann1_cache = $rgend1_cache->annotate();

my $rgend2_cache = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry    => $entry,
							 -end      => 2,
							 -genome_data => $genome_data_cache,
							 -within   => $within );
my $ann2_cache = $rgend2_cache->annotate();

my $combi_cache = new Sanger::CGP::Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1_cache->[0],
								-anno2   => $ann2_cache->[0],
								-strand1 => $strand1,
								-strand2 => $strand2,
								-shard   => $shard);
$combi_cache->combine(); 
my $rganno_cache = $combi_cache->anno();

my $res = get_result_58();
my $res_74 = get_result_74();

ok($combi, 'object defined. ENSEMBL');
is (($combi->strand1()), $strand1 , "get strand1. ENSEMBL");
is (($combi->strand2()), $strand2 , "get strand2. ENSEMBL");
is (($combi->shard()), $shard , "get shard. ENSEMBL");

is ((Dumper($rganno)), $res_74 , "get result. ENSEMBL");

ok($combi_cache, 'object defined. CACHE');
is (($combi_cache->strand1()), $strand1 , "get strand1. CACHE");
is (($combi_cache->strand2()), $strand2 , "get strand2. CACHE");
is (($combi_cache->shard()), $shard , "get shard. CACHE");

is ((Dumper($rganno_cache)), $res_74 , "get result. CACHE");

#------------------------------------------------------------------------------------------------#
sub get_result_58 {
    my $res = <<'END';
$VAR1 = bless( {
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
    my $res = <<'END';
$VAR1 = bless( {
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
