#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
# 
# Author: Lucy Stebbings <cgpit@sanger.ac.uk>
# 
# This file is part of cgpPindel.
# 
# cgpPindel is free software: you can redistribute it and/or modify it under
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

use Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl;

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


my $species = 'HUMAN';
my $ensembl_api_58 = '/software/pubseq/PerlModules/Ensembl/www_58_1';
my $ensembl_api_74 = '/software/pubseq/PerlModules/Ensembl/www_74_1';

my $version = '74';
my $ensembl_api = $ensembl_api_74;

my $genome_data_ensembl = new Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl(-species     => $species,
								   -ensembl_api => $ensembl_api,
								   -gene_id_required => 1,
								   -use_all_biotypes => 1);

ok($genome_data_ensembl, 'object defined');
is (($genome_data_ensembl->gene_id_required()), 1 , "get gene_id_required");
is (($genome_data_ensembl->use_all_biotypes()), 1 , "get use_all_biotypes");
is (($genome_data_ensembl->species()), $species , "get species");
is (($genome_data_ensembl->ensembl_api()), $ensembl_api , "get ensembl_api");
ok($genome_data_ensembl->registry(), 'registry object defined');

test_between_ensembl($genome_data_ensembl, $version);
test_fetch_transcript_by_region_ensembl($genome_data_ensembl, $version);

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
sub test_fetch_transcript_by_region_ensembl {
    my ($genome_data, $version) = @_;

    my $chr = '3';
    my $pos_start = 129389225;
    my $pos_end = 129389230;

    my $fetch_region = 'Sanger::CGP::Grass::GenomeData::Transcript';

    my $all_count = 0;
    if ($version eq '58') {  $all_count = 5; }
    if ($version eq '74') { $all_count = 4; } # one of these is a putative coding transcript so you would only get 3 if use_all_biotypes = 0

    my $fetch_region_out = $genome_data->fetch_transcripts_by_region($chr, $pos_start, $pos_end);
    my $translist_count = scalar(@$fetch_region_out);
    ok($genome_data->trans_ad(), 'transcript adaptor object defined. ENSEMBL');
    is (ref($fetch_region_out->[0]), $fetch_region , "get fetch_region transcript object. ENSEMBL");
    is ($translist_count , $all_count , "get fetch_region transcript count. ENSEMBL");
}
#------------------------------------------------------------------------------------------------#
sub get_result {
    my $res = <<'END';
END
    return($res);
}
