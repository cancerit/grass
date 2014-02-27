#!/usr/bin/env perl

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
my $species = 'HUMAN';
my $ensembl_api = 58;
my $registry = set_ensembl($ensembl_api, $species);


my $rgend1 = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
						   -end      => 1,
						   -species  => $species,
						   -registry => $registry,
						   -within   => $within );
my $ann1 = $rgend1->annotate();

my $rgend2 = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
						   -end      => 2,
						   -species  => $species,
						   -registry => $registry,
						   -within   => $within );
my $ann2 = $rgend2->annotate();

my $combi = new Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1->[0],
							  -anno2   => $ann2->[0],
							  -strand1 => $strand1,
							  -strand2 => $strand2,
							  -shard   => $shard);
$combi->combine(); 
my $rganno = $combi->anno();

my $res = get_result();

ok($combi, 'object defined');
is (($combi->strand1()), $strand1 , "get strand1");
is (($combi->strand2()), $strand2 , "get strand2");
is (($combi->shard()), $shard , "get shard");

is ((Dumper($rganno)), $res , "get result");

#------------------------------------------------------------------------------------------------#
sub set_ensembl {
    my $ensembl_api = shift; # either passed in or from the CGP DB
    my $species = shift;

    my $ensembl = '';
    if    ($ensembl_api =~ /^\d+$/)         { $ensembl = "www_" . $ensembl_api . "_1"; }
    elsif ($ensembl_api =~ /^\d+_\d+$/)     { $ensembl = "www_" . $ensembl_api; }
    elsif ($ensembl_api =~ /^www_\d+_\d+$/) { $ensembl =          $ensembl_api; }
    else                                    { print "ensembl_api format $ensembl_api not recognised\n"; exit; }
    print "$ensembl\n";

# put this here because may need to pass in which ensembl version to use
    my $sent = 'use lib qw(/software/pubseq/PerlModules/Ensembl/' . $ensembl . '/ensembl/modules/  );
                use Bio::EnsEMBL::Registry;';
    eval $sent; # delays it until run time so that $ensembl is set first

    #print "@INC\n";
    # team ENSEMBL install (fast, build 37)
    my $registry = 'Bio::EnsEMBL::Registry';

    # local sanger CGP ensembl install
    $registry->clear();
    $registry->load_registry_from_db( -host => 'cgp-ensembl-db.internal.sanger.ac.uk',
				      -user => 'cgpense-ro');

    unless ($registry->get_adaptor($species,'core','slice')) {
	if ($species eq 'CEREVISIAE') { $species = 'SACCHAROMYCES'; }
	$registry->load_registry_from_db( -host => 'cgp-ensembl-db.internal.sanger.ac.uk',
					  -user => 'cgpense-ro');
    }
    unless ($registry->get_adaptor($species,'core','slice')) {
	print "use remote ensembl server\n";
	if ($species eq 'CEREVISIAE') { $species = 'Saccharomyces cerevisiae'; }
	$registry->clear();
	$registry->load_registry_from_db( -host => 'ensembldb.ensembl.org',
				          -user => 'anonymous');
    }
    unless ($registry->get_adaptor($species,'core','slice')) {
	print "could not get connection to local or remote ensembl registry\n"; 
	exit;
    }
    return($registry);
}
#------------------------------------------------------------------------------------------------#
sub get_result {
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
