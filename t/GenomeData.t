#!/usr/bin/perl

# for testing Sanger::CGP::Grass::GenomeData class

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Data::Dumper;

use Sanger::CGP::Grass::GenomeData;

use Test::More 'no_plan';

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


my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $species = 'HUMAN';
my $ensembl_api_58 = '/software/pubseq/PerlModules/Ensembl/www_58_1';
my $ensembl_api_74 = '/software/pubseq/PerlModules/Ensembl/www_74_1';

my $genome_data_ensembl = new Sanger::CGP::Grass::GenomeData(-species     => $species,
						-ensembl_api => $ensembl_api_74,
						-gene_id_required => 0,
						-use_all_biotypes => 0);

my $genome_data_cache = new Sanger::CGP::Grass::GenomeData(-genome_cache => $genome_cache,
					      -gene_id_required => 0);

ok($genome_data_ensembl, 'object defined');
is (($genome_data_ensembl->species()), $species , "get species");
is (($genome_data_ensembl->ensembl_api()), $ensembl_api_74 , "get ensembl_api");
ok($genome_data_ensembl->registry(), 'registry object defined');

#test_ensembl($genome_data_ensembl, '74');
test_cache($genome_data_cache);

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
sub test_ensembl {
    my ($genome_data, $version) = @_;
    test_between_ensembl($genome_data, $version);
    test_fetch_transcript_by_region_ensembl($genome_data, $version);
    test_thin_out_translist_ensembl($genome_data, $version);

}
#------------------------------------------------------------------------------------------------#
sub test_cache {
    my ($genome_data) = @_;

    test_between_cache($genome_data);
    test_fetch_transcript_by_region_cache($genome_data);
    test_thin_out_translist_cache($genome_data);
    #test_5prime_truncated_cache($genome_data);
}
#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
sub test_between_ensembl {
    my ($genome_data) = @_;

    my $chr = '2';
    my $pos_start = 188365837;
    my $pos_end = 188417763;
    
    my $between = 'ENSG00000224063,TFPI';
    
    my $between_out = $genome_data->get_gene_list($chr, $pos_start, $pos_end);
    
    ok($genome_data->slice_ad(), 'slice adaptor object defined. ENSEMBL');
    is ($between_out, $between , "get between genes. ENSEMBL");
}
#------------------------------------------------------------------------------------------------#
sub test_between_cache {
    my ($genome_data) = @_;

    my $chr = '2';
    my $pos_start = 188365837;
    my $pos_end = 188417763;
    
    my $between = 'TFPI';
    
    my $between_out = $genome_data->get_gene_list($chr, $pos_start, $pos_end);
    
    is (($between), $between_out , "get between genes. CACHE");
}
#------------------------------------------------------------------------------------------------#
sub test_fetch_transcript_by_region_ensembl {
    my ($genome_data, $version) = @_;

    my $chr = '3';
    my $pos_start = 129389225;
    my $pos_end = 129389230;

    my $fetch_region = 'Sanger::CGP::Grass::GenomeData::Transcript';

    my $all_count = 5;
    if ($version eq '74') {
	$all_count = 3;
   }
    my $fetch_region_out = $genome_data->fetch_transcripts_by_region($chr, $pos_start, $pos_end);
    my $translist_count = scalar(@$fetch_region_out);
    ok($genome_data->trans_ad(), 'transcript adaptor object defined. ENSEMBL');
    is (ref($fetch_region_out->[0]), $fetch_region , "get fetch_region transcript object. ENSEMBL");
    is ($translist_count , $all_count , "get fetch_region transcript count. ENSEMBL");
}
#------------------------------------------------------------------------------------------------#
sub test_fetch_transcript_by_region_cache {
    my ($genome_data) = @_;

    my $chr = '3';
    my $pos_start = 129389225;
    my $pos_end = 129389230;

    my $fetch_region = 'Sanger::CGP::Grass::GenomeData::Transcript';
    my $all_count = 3;

    my $fetch_region_out = $genome_data->fetch_transcripts_by_region($chr, $pos_start, $pos_end);
    my $translist_count = scalar(@$fetch_region_out);

    is (ref($fetch_region_out->[0]), $fetch_region , "get fetch_region cache transcript object. CACHE");
    is ($translist_count , $all_count , "get fetch_region transcript count. CACHE");
}
#------------------------------------------------------------------------------------------------#
sub test_thin_out_translist_ensembl {
    my ($genome_data, $version) = @_;

    my $chr = '3';
    my $pos_start = 129389225;
    my $pos_end = 129389225;

    my $all_count = 5;
    my $thinned_count = 4;
    my $is_ccds = 1;

    if ($version eq '74') {
	$all_count = 3;
	$thinned_count = 1;
   }

    my $full_translist = $genome_data->fetch_transcripts_by_region($chr, $pos_start, $pos_end);
    my $all_count_out = scalar(@$full_translist);

    my ($thinned_translist, $ccds_only) = $genome_data->thin_out_translist($full_translist);
    my $thinned_count_out = scalar(@$thinned_translist);

    is (($all_count_out), $all_count , "get thin_out_transcript all count. ENSEMBL");
    is (($thinned_count_out), $thinned_count , "get thin_out_transcript thinned count. ENSEMBL");
    is (($ccds_only), $is_ccds , "get thin_out_transcript ccds_only. ENSEMBL");
}
#------------------------------------------------------------------------------------------------#
sub test_thin_out_translist_cache {
    my ($genome_data) = @_;

    my $chr = '3';
    my $pos_start = 129389225;
    my $pos_end = 129389225;

    my $all_count = 3; # 3 not 4 - cache doesn't contain any putative_coding transcripts - ensembl version does
    my $thinned_count = 1;
    my $is_ccds = 1;

    my $full_translist = $genome_data->fetch_transcripts_by_region($chr, $pos_start, $pos_end);
    my $all_count_out = scalar(@$full_translist);

    my ($thinned_translist, $ccds_only) = $genome_data->thin_out_translist($full_translist);
    my $thinned_count_out = scalar(@$thinned_translist);

    is (($all_count_out), $all_count , "get thin_out_transcript all count. CACHE");
    is (($thinned_count_out), $thinned_count , "get thin_out_transcript thinned count. CACHE");
    is (($ccds_only), $is_ccds , "get thin_out_transcript ccds_only. CACHE");
}
#------------------------------------------------------------------------------------------------#
sub test_5prime_truncated_cache {
    my ($genome_data) = @_;

    my $chr = '19';
    my $pos_start = 47011780;
    my $pos_end = 47012370;

    my $all_count = 1; # 3 not 4 - cache doesn't contain any putative_coding transcripts - ensembl version does
    my $phase = 1;

    my $full_translist = $genome_data->fetch_transcripts_by_region($chr, $pos_start, $pos_end);
    my $all_count_out = (scalar(@$full_translist));
#    print "TEST: " . $full_translist->[0]->stable_id . " " . $full_translist->[0]->exons->[0]->phase . "\n";


    my ($thinned_translist, $ccds_only) = $genome_data->thin_out_translist($full_translist);
    my $thinned_count_out = scalar(@$thinned_translist);

    is (($all_count_out), $all_count , "get 5prime_truncated all count. CACHE");

    # this one should work but there might be an issue with start phase on 5' truncated transcripts
    is (($full_translist->[0]->exons->[0]->phase), $phase , "get start phase. CACHE");
}
#------------------------------------------------------------------------------------------------#
sub get_result {
    my $res = <<'END';
END
    return($res);
}
