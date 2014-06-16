#!/usr/bin/perl

# for testing Sanger::CGP::Grass::DataFlankingBases class

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::FlankingBases;

use File::Copy qw(copy);
use Text::Diff;
use Capture::Tiny qw(capture);

use Test::More 'no_plan';

my $ref = '/nfs/cancer_ref01/human/37/genome.fa';

test_file_input($ref);

#----------------------------------------------------------------------------------------#

sub test_file_input {
    my ($ref) = @_;

    my $infile = dirname(abs_path($0)).'/../testData/' . 'test_FlankingBases_ann.in';
    my $outfile = dirname(abs_path($0)).'/../testData/' . 'test_FlankingBases_ann.out';

    my $testfile = 'test_FlankingBases_ann.bedpe';

    copy $infile, $testfile;

# make a new object
    my $FlankingBases = new Sanger::CGP::Grass::FlankingBases(-infile => $testfile,
							      -ref    => $ref );
    ok($FlankingBases, 'object defined');

    $FlankingBases->process();

    my $diff = diff "$testfile", "$outfile";

    ok(!$diff, 'correct file created');
    
    unless ($diff) { unlink $testfile; }
}
#----------------------------------------------------------------------------------------#
