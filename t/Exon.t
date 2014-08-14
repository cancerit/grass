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


# for testing Sanger::CGP::Grass::GenomeData::Exon class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::GenomeData::Exon;

use Test::More 'no_plan';

# existing entry


my $phase = 0;
my $end_phase = -1;
my $start = 50;
my $end = 100;
my $coding_region_start = 60;
my $coding_region_end = 99;
my $seq = 'acgtACGTNn';

# make a new object
my $Exon = new Sanger::CGP::Grass::GenomeData::Exon(-phase => $phase,
				       -end_phase => $end_phase,
				       -start => $start,
				       -end => $end,
				       -coding_region_start => $coding_region_start,
				       -coding_region_end => $coding_region_end,
				       -seq => $seq );

ok($Exon, 'object defined');

is (($Exon->phase()), $phase , "get phase");
is (($Exon->end_phase()), $end_phase , "get end_phase");
is (($Exon->start()), $start , "get start");
is (($Exon->end()), $end , "get end");
is (($Exon->coding_region_start()), $coding_region_start , "get coding_region_start");
is (($Exon->coding_region_end()), $coding_region_end , "get coding_region_end");
is (($Exon->seq()), $seq , "get seq");
