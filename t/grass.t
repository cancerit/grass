#!/usr/bin/perl

# for testing grass.pl script

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use Time::localtime;
use File::Copy qw(copy);
use Text::Diff;
use Capture::Tiny qw(capture);

use strict;
use warnings FATAL => 'all';

use Test::More 'no_plan';

# got to do something about this...
my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $ref = '/nfs/cancer_ref01/human/37/genome.fa';

my $script = dirname(abs_path($0)).'/../bin/' . 'grass.pl';

test_coord_input();
test_file_input();
test_file_input_bedpe();

#----------------------------------------------------------------------------------------#
sub test_coord_input {

    my $test_coord = 'X:-:84500101-84500632,X:+:84510265-84510265';
    my $test_coord_res = 'gene	gene_id	transcript_id	strand	end_phase	region	region_number	total_region_count	first/last	gene	gene_id	transcript_id	strand	phase	region	region_number	total_region_count	first/last	fusion_flag
ZNF711	ZNF711	ENST00000373165	1	_	5UTRintron	1	9	_	ZNF711	ZNF711	ENST00000373165	1	1	exon	4	9	first_base	570
';

    my ($test_coord_out, @result) = capture { system("perl $script -genome_cache $genome_cache -coord $test_coord -ref $ref "); };

    is ($test_coord_out, $test_coord_res, "check supply coord");

}
#----------------------------------------------------------------------------------------#
sub test_file_input {

    my $infile = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I_and_II';
    my $outfile = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I_and_II_ann';

    my $testfile = 'testout_Brass_I_and_II';
    my $testfile_out = 'testout_Brass_I_and_II_ann';

    copy $infile, $testfile;
    my ($test_file_res, @result) = capture { system("perl $script -genome_cache $genome_cache -file $testfile"); };
    if ($test_file_res) { print "OUTPUT AND ERRORS: " . $test_file_res . "@result\n"; }

    my $diff = diff "$testfile_out", "$outfile";

    ok(!$diff, 'correct file created');
    
    unless ($diff) { unlink $testfile; unlink $testfile_out; }
}
#----------------------------------------------------------------------------------------#
sub test_file_input_bedpe {

    my $infile = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I.bedpe';
    my $outfile = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I_ann.bedpe';
    my $outfile_vcf = dirname(abs_path($0)).'/../testData/' . 'testout_Brass_I_ann.vcf';

    my $testfile = 'testout_Brass_I.bedpe';
    my $testfile_out = 'testout_Brass_I_ann.bedpe';
    my $testfile_out_vcf = 'testout_Brass_I_ann.vcf';

    copy $infile, $testfile;
    my ($test_file_res, @result) = capture { system("perl $script -genome_cache $genome_cache -file $testfile -ref $ref "); };
    if ($test_file_res) { print "OUTPUT AND ERRORS: " . $test_file_res . "@result\n"; }

    # the current file date is put in the output vcf file so it will always differ to the 'correct' file in testData unless updated
    my $date_in_correct_testfile = '20140512';
    my $tm = localtime;
    my $today = ($tm->year() + 1900) . sprintf("%02d",($tm->mon() + 1)) . sprintf("%02d",($tm->mday()));
    # do an inline substitution
    my $command1 = "perl -pi -e \'s/##fileDate=$today/##fileDate=$date_in_correct_testfile/g\' $testfile_out_vcf";
    my $command2 = "perl -pi -e \'s/##source_$today/##source_$date_in_correct_testfile/g\' $testfile_out_vcf";
    my $command3 = "perl -pi -e \'s/##vcfProcessLog_$today/##vcfProcessLog_$date_in_correct_testfile/g\' $testfile_out_vcf";
    `$command1`;
    `$command2`;
    `$command3`;

    my $diff1 = diff "$testfile_out", "$outfile";
    my $diff2 = diff "$testfile_out_vcf", "$outfile_vcf";

    ok(!$diff1, 'correct bedpe file created');
    ok(!$diff2, 'correct vcf file created');
    
    unless ($diff1 || $diff2) { unlink $testfile; unlink $testfile_out; unlink $testfile_out_vcf; }
}
#----------------------------------------------------------------------------------------#
