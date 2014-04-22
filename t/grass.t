#!/usr/bin/perl

# for testing grass.pl script

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use File::Copy qw(copy);
use Text::Diff;
use Capture::Tiny qw(capture);

use strict;
use warnings FATAL => 'all';

use Test::More 'no_plan';

# got to do something about this...
my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $vagrent = '/software/CGP/projects/vagrent_vcf/perl/modules/';

my $script = dirname(abs_path($0)).'/../bin/' . 'grass.pl';

test_coord_input();
test_file_input();

#----------------------------------------------------------------------------------------#
sub test_coord_input {

    my $test_coord = 'X:-:84500101-84500632,X:+:84510265-84510265';
    my $test_coord_res = 'gene	gene_id	transcript_id	strand	end_phase	region	region_number	total_region_count	first/last	gene	gene_id	transcript_id	strand	phase	region	region_number	total_region_count	first/last	fusion_flag
ZNF711	ZNF711	ENST00000373165	1	_	5UTRintron	1	9	_	ZNF711	ZNF711	ENST00000373165	1	1	exon	4	9	first_base	570
';

    my ($test_coord_out, @result) = capture { system("perl -I $vagrent $script -genome_cache $genome_cache -coord $test_coord"); };
    #my $test_coord_out = `perl -I $vagrent $script -genome_cache $genome_cache -coord $test_coord`;

    is ($test_coord_out, $test_coord_res, "check supply coord");

}
#----------------------------------------------------------------------------------------#
sub test_file_input {

    my $infile = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I_and_II';
    my $outfile = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I_and_II_ann';

    my $testfile = 'testout_Brass_I_and_II';
    my $testfile_out = 'testout_Brass_I_and_II_ann';

    copy $infile, $testfile;

    my ($test_file_res, @result) = capture { system("perl -I $vagrent $script -genome_cache $genome_cache -file $testfile"); };
    #my $test_file_res = `perl -I $vagrent $script -genome_cache $genome_cache -file $testfile`;

    my $diff = diff "$testfile_out", "$outfile";

    ok(!$diff, 'correct file created');
    
    unless ($diff) { unlink $testfile; unlink $testfile_out; }

}
#----------------------------------------------------------------------------------------#
