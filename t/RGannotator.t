#!/usr/bin/env perl

# for testing Grass::RGannotator class

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

use Grass::DataEntry;
use Grass::Annotation::RGannotator;

use Test::More 'no_plan';

# make a dataset
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
my $dataset1 = [];
push @$dataset1, $entry;

my $name_2 = 37464190;
my $chr1_2 = '2';
my $strand1_2 = '+';
my $pos1_start_2 = 188365485;
my $pos1_end_2 = 188365837;
my $chr2_2 = '2';
my $strand2_2 = '+';
my $pos2_start_2 = 188417763;
my $pos2_end_2 = 188418155;
my $shard_2 = '';
my $count_2 = 4;
my $entry_2 = new Grass::DataEntry(-name       => $name_2,
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
my $dataset2 = [];
push @$dataset2, $entry_2;

# set values
my $within = 0;
my $species = 'HUMAN';
my $ensembl_api = 58;
my $registry = set_ensembl($ensembl_api, $species);

# results should be...
my $res1 = "TMCC1\tENSG00000172765\tENST00000329333\t-1\t1\texon\t5\t7\t_\tTMCC1\tENSG00000172765\tENST00000329333\t-1\t1\texon\t5\t7\t_\t715\nTMCC1\tENSG00000172765\tENST00000307824\t-1\t1\texon\t5\t7\t_\tTMCC1\tENSG00000172765\tENST00000307824\t-1\t_\t5UTRexon\t5\t7\t_\t100\n";
my $res2 = "TFPI\tENSG00000003436\tENST00000233156\t-1\t1\tintron\t2\t8\t_\tprotein_coding\tTFPI\tENSG00000003436\tENST00000233156\t-1\t_\t5UTRintron\t1\t8\t_\tprotein_coding\t200\tENSG00000224063,TFPI\n";

# make a new object
my $rgann1 = new Grass::Annotation::RGannotator(-dataset         => $dataset1,
						-within          => $within,
						-species         => $species,
						-registry        => $registry,
						-list_between    => 0,
						-show_biotype    => 0 );
$rgann1->getAnnotation();

# make a new object
my $rgann2 = new Grass::Annotation::RGannotator(-dataset         => $dataset2,
						-within          => $within,
						-species         => $species,
						-registry        => $registry,
						-list_between    => 1,
						-show_biotype    => 1 );
$rgann2->getAnnotation();

ok($rgann1, 'object defined');
ok($rgann2, 'object defined');

is (($rgann1->within()), $within , "get within");
is (($rgann1->species()), $species , "get species");
is (($rgann1->registry()), $registry , "get registry");
is (($rgann1->list_between()), 0 , "get list_between");
is (($rgann1->show_biotype()), 0 , "get show_biotype");
is (($rgann1->format_for_printing()), $res1 , "get result");

is (($rgann2->within()), $within , "2get within");
is (($rgann2->species()), $species , "2get species");
is (($rgann2->registry()), $registry , "2get registry");
is (($rgann2->list_between()), 1 , "2get list_between");
is (($rgann2->show_biotype()), 1 , "2get show_biotype");
is (($rgann2->format_for_printing()), $res2 , "2get result");

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
