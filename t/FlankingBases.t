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


# for testing Sanger::CGP::Grass::DataFlankingBases class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::FlankingBases;

use File::Copy qw(copy);
use Text::Diff;
use Capture::Tiny qw(capture);
use FindBin qw($Bin);

use Test::More 'no_plan';

my $ref = '/nfs/cancer_ref01/human/37/genome.fa';

test_file_input($ref);

#----------------------------------------------------------------------------------------#

sub test_file_input {
    my ($ref) = @_;

    my $infile = $Bin.'/../testData/' . 'test_FlankingBases_ann.in';
    my $outfile = $Bin.'/../testData/' . 'test_FlankingBases_ann.out';

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
