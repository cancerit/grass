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


# for testing Sanger::CGP::Grass::GenomeDataCache class

use strict;
use warnings FATAL => 'all';

use Data::Dumper;

use Sanger::CGP::Grass::GenomeData::GenomeDataCache;

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


my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz';

my $genome_data_cache = new Sanger::CGP::Grass::GenomeData::GenomeDataCache(-genome_cache => $genome_cache,
							       -gene_id_required => 0);

ok($genome_data_cache, 'object defined');
is (($genome_data_cache->genome_cache()), $genome_cache , "get genome_cache");
is (($genome_data_cache->gene_id_required()), 0 , "get gene_id_required");

test_between_cache($genome_data_cache);
test_fetch_transcript_by_region_cache($genome_data_cache);
test_5prime_truncated_cache($genome_data_cache);

#------------------------------------------------------------------------------------------------#
sub test_between_cache {
    my ($genome_data) = @_;

    my $chr = '2';
    my $pos_start = 188365837;
    my $pos_end = 188417763;

    my $between = 'TFPI';

    my $between_out = $genome_data->get_gene_list($chr, $pos_start, $pos_end);

    is (($between_out), $between , "get between genes. CACHE");
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

    is (($all_count_out), $all_count , "get 5prime_truncated all count. CACHE");

    # this one *should* work but there might be an issue with start phase on 5' truncated transcripts INVESTIGATE
#    is (($full_translist->[0]->exons->[0]->phase), $phase , "get start phase. CACHE");
}
#------------------------------------------------------------------------------------------------#
sub get_result {
    my $res = <<'END';
END
    return($res);
}
