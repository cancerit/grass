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


# for testing Sanger::CGP::Grass::DataFlankingBases class

use strict;
use warnings FATAL => 'all';

use File::Copy qw(copy);
use Capture::Tiny qw(capture);
use FindBin qw($Bin);
use File::Temp qw(tempdir);

use Test::More 'no_plan';

for my $module (qw(Sanger::CGP::Grass::FlankingBases)) {
  require_ok $module or BAIL_OUT "Can't load $module";
}

my $ref = q{};

if(exists $ENV{GRASS_GRCH37_FA}) {
  $ref = $ENV{GRASS_GRCH37_FA};
}

# to look at the created file, set CLEANUP to 0 and uncomment the note line to see its location in the output tmp space
my $tmpdir = tempdir( CLEANUP => 1 );
#note $tmpdir;

SKIP: {
  unless(-e $ref) {
    my $message = "SKIPPING: Reference *.fa not found at '$ref', set ENV: GRASS_GRCH37_FA to enable these tests";
    warn "$message\n";
    skip $message, '1';
  }
  subtest 'flanking bases' => sub {
    my $infile = "$Bin/../testData/test_FlankingBases_ann.in";
    my $outfile = "$Bin/../testData/test_FlankingBases_ann.out";

    my $testfile = "$tmpdir/test_FlankingBases_ann.bedpe";

    copy $infile, $testfile;

    # make a new object
    my $FlankingBases = new_ok('Sanger::CGP::Grass::FlankingBases',
                  [-infile => $testfile,
                   -ref    => $ref]);

    $FlankingBases->process();

    my $slurp_tf = slurp_file($testfile);
    my $slurp_of = slurp_file($outfile);

    is_deeply($slurp_tf, $slurp_of, 'correct file created');
  };
};

sub slurp_file {
  my $in = shift;
  open my $fh, '<', $in or die $!;
  my @data = <$fh>;
  close $fh;
  return \@data;
}
