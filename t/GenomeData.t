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


# for testing Sanger::CGP::Grass::GenomeData class

use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use FindBin qw($Bin);

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

for my $module (qw(Sanger::CGP::Grass::GenomeData)) {
  require_ok $module or BAIL_OUT "Can't load $module";
}

my $genome_cache = "$Bin/../testData/vagrent.cache.gz";
my $species = 'HUMAN';
#my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_74_1';
my $ensembl_api = '';
if(exists $ENV{GRASS_ENS_API}) {
  $ensembl_api = $ENV{GRASS_ENS_API};
}

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
                -use_all_biotypes => 0]);
    is (($genome_data_ensembl->species()), $species , "get species");
    is (($genome_data_ensembl->ensembl_api()), $ensembl_api , "get ensembl_api");
    ok($genome_data_ensembl->registry(), 'registry object defined');
    test_ensembl($genome_data_ensembl, '74');
  };
};

my $genome_data_cache = new_ok('Sanger::CGP::Grass::GenomeData',
              [-genome_cache => $genome_cache,
					      -gene_id_required => 0]);
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
