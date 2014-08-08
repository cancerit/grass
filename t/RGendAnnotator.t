#!/usr/bin/perl

# for testing Sanger::CGP::Grass::Annotation::RGendAnnotator class

use strict;
use warnings FATAL => 'all';

use Data::Dumper;

use Sanger::CGP::Grass::GenomeData;
use Sanger::CGP::Grass::DataEntry;
use Sanger::CGP::Grass::Annotation::RGendAnnotator;

use Test::More 'no_plan';

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
my $entry_2 = new Sanger::CGP::Grass::DataEntry(-name       => $name_2,
				 -chr1       => $chr1_2,
				 -strand1    => $strand1_2,
				 -pos1_start => $pos1_start_2,
				 -pos1_end   => $pos1_end_2,
				 -chr2       => $chr2_2,
				 -strand2    => $strand2_2,
				 -pos2_start => $pos2_start_2,
				 -pos2_end   => $pos2_end_2,
				 -shard      => $shard_2,
				 -count      => $count_2 );

my $end_2 = 1;
my $within_2 = 100;

my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $species = 'HUMAN';
#my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_58_1';
my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_74_1';

my $genome_data_ensembl = new Sanger::CGP::Grass::GenomeData(-species     => $species,
						-ensembl_api => $ensembl_api,
						-gene_id_required => 0,
                                                -use_putative_coding => 0);

my $genome_data_cache = new Sanger::CGP::Grass::GenomeData(-genome_cache => $genome_cache,
						-gene_id_required => 0,
                                                -use_putative_coding => 0);

# make a new object
my $End = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry   => $entry,
						-end     => $end,
						-within  => $within,
						-genome_data => $genome_data_ensembl);
my $anns = $End->annotate();
my $End2 = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry   => $entry,
						 -end     => $end,
						 -within  => $within,
						 -genome_data => $genome_data_cache);
my $anns2 = $End2->annotate();

# make a new object - plus strand
my $End_2 = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry   => $entry_2,
						-end     => $end_2,
						-within  => $within_2,
						-genome_data => $genome_data_ensembl);
my $anns_2 = $End_2->annotate();
my $End2_2 = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry   => $entry_2,
						 -end     => $end_2,
						 -within  => $within_2,
						 -genome_data => $genome_data_cache);
my $anns2_2 = $End2_2->annotate();

#my $res = get_result_58();
my $res = get_result_74();

#my $res_2 = get_result_58_2();
my $res_2 = get_result_74_2();

ok($End, 'object defined');
is (($End->entry()), $entry , "get entry");
is (($End->end()), $end , "get end");
is (($End->within()), $within , "get within");
is (($End->ccds_only()), 1 , "get ccds_only");
is ((Dumper($anns)), $res , "get result ensembl");

is ((Dumper($anns2)), $res , "get result cache");


# check on the plus strand
ok($End_2, 'object defined _2');
is (($End_2->entry()), $entry_2 , "get entry (plus strand gene)");
is (($End_2->end()), $end_2 , "get end (plus strand gene)");
is (($End_2->within()), $within_2 , "get within (plus strand gene)");
is (($End_2->ccds_only()), 1 , "get ccds_only (plus strand gene)");
is ((Dumper($anns_2)), $res_2 , "get result ensembl (plus strand gene)"); # has a putative coding transcript as well that is missing from the cache version
is ((Dumper($anns2_2)), $res_2 , "get result cache (plus strand gene)");
#print Dumper($anns2_2);

#------------------------------------------------------------------------------------------------#
sub get_result_58 {
    my $res = <<'END';
$VAR1 = [
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
    my $res = <<'END';
$VAR1 = [
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
    my $res = <<'END';
$VAR1 = [
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
