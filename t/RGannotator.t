#!/usr/bin/perl

# for testing Sanger::CGP::Grass::RGannotator class

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::GenomeData;
use Sanger::CGP::Grass::DataEntry;
use Sanger::CGP::Grass::Annotation::RGannotator;

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
my $dataset2 = [];
push @$dataset2, $entry_2;

# set values
my $within = 0;
my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $species = 'HUMAN';
my $ensembl_api = '/software/pubseq/PerlModules/Ensembl/www_58_1';

# set up access to genome data  from EnsemblDB (if species/ensembl_api supplied) or a cached flat file version (if genome_cache supplied)
my $genome_data_ensembl = new Sanger::CGP::Grass::GenomeData(-species     => $species,
						-ensembl_api => $ensembl_api,
                                                -gene_id_required => 1);

my $genome_data_cache = new Sanger::CGP::Grass::GenomeData(-genome_cache => $genome_cache,
					      -gene_id_required => 0);

# results should be...
my $res1 = "TMCC1\tENSG00000172765\tENST00000329333\t-1\t1\texon\t5\t7\t_\tTMCC1\tENSG00000172765\tENST00000329333\t-1\t1\texon\t5\t7\t_\t715\nTMCC1\tENSG00000172765\tENST00000307824\t-1\t1\texon\t5\t7\t_\tTMCC1\tENSG00000172765\tENST00000307824\t-1\t_\t5UTRexon\t5\t7\t_\t100\n";
my $res2 = "TFPI\tENSG00000003436\tENST00000233156\t-1\t1\tintron\t2\t8\t_\tprotein_coding\tTFPI\tENSG00000003436\tENST00000233156\t-1\t_\t5UTRintron\t1\t8\t_\tprotein_coding\t200\tENSG00000224063,TFPI\n";
my $res3 = "TFPI\tTFPI\tENST00000233156\t-1\t1\tintron\t2\t8\t_\tprotein_coding\tTFPI\tTFPI\tENST00000233156\t-1\t_\t5UTRintron\t1\t8\t_\tprotein_coding\t200\tTFPI\n";

# make a new object
my $rgann1 = new Sanger::CGP::Grass::Annotation::RGannotator(-dataset         => $dataset1,
						-within          => $within,
						-genome_data     => $genome_data_ensembl,
						-list_between    => 0,
						-show_biotype    => 0 );
$rgann1->getAnnotation();

# make a new object
my $rgann2 = new Sanger::CGP::Grass::Annotation::RGannotator(-dataset         => $dataset2,
						-within          => $within,
						-genome_data     => $genome_data_ensembl,
						-list_between    => 1,
						-show_biotype    => 1 );
$rgann2->getAnnotation();

my $rgann3 = new Sanger::CGP::Grass::Annotation::RGannotator(-dataset         => $dataset2,
						-within          => $within,
						-genome_data    => $genome_data_cache,
						-list_between    => 1,
						-show_biotype    => 1 );
$rgann3->getAnnotation();

ok($rgann1, 'object defined. ENSEMBL');
ok($rgann2, '2object defined. ENSEMBL');
ok($rgann3, 'object defined');

is (($rgann1->within()), $within , "get within. ENSEMBL");
is (($rgann1->list_between()), 0 , "get list_between. ENSEMBL");
is (($rgann1->show_biotype()), 0 , "get show_biotype. ENSEMBL");
is (($rgann1->format_for_printing()), $res1 , "get result. ENSEMBL");

is (($rgann2->within()), $within , "2get within. ENSEMBL");
is (($rgann2->list_between()), 1 , "2get list_between. ENSEMBL");
is (($rgann2->show_biotype()), 1 , "2get show_biotype. ENSEMBL");
is (($rgann2->format_for_printing()), $res2 , "2get result. ENSEMBL");

is (($rgann3->within()), $within , "3get within. CACHE");
is (($rgann3->list_between()), 1 , "3get list_between. CACHE");
is (($rgann3->show_biotype()), 1 , "3get show_biotype. CACHE");
is (($rgann3->format_for_printing()), $res3 , "3get result. CACHE");

#------------------------------------------------------------------------------------------------#
