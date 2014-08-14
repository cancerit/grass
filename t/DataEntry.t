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


# for testing Sanger::CGP::Grass::DataEntry class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::DataEntry;

use Test::More 'no_plan';

my $name = 'TEST';
my $chr1 = '3';
my $strand1 = '+';
my $pos1_start = 23;
my $pos1_end = 45;
my $chr2 = 'X';
my $strand2 = '-';
my $pos2_start = 67;
my $pos2_end = 89;
my $shard = 'AT';
my $count = 6;

my $coord = $chr1 . ':' . $strand1 . ':' . $pos1_start . '-' . $pos1_end . ',' . $chr2 . ':' . $strand2 . ':' . $pos2_start . '-' . $pos2_end . ',' . $shard;

# make a new object
my $Entry = new Sanger::CGP::Grass::DataEntry(-name       => $name,
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

my $Entry_coord = new Sanger::CGP::Grass::DataEntry(-coord => $coord );

ok($Entry, 'object defined');

is (($Entry->name()), $name , "get name");
is (($Entry->chr1()), $chr1 , "get chr1");
is (($Entry->strand1()), $strand1 , "get strand1");
is (($Entry->pos1_start()), $pos1_start , "get pos1_start");
is (($Entry->pos1_end()), $pos1_end , "get pos1_end");
is (($Entry->chr2()), $chr2 , "get chr2");
is (($Entry->strand2()), $strand2 , "get strand2");
is (($Entry->pos2_start()), $pos2_start , "get pos2_start");
is (($Entry->pos2_end()), $pos2_end , "get pos2_end");
is (($Entry->shard()), $shard , "get shard");
is (($Entry->count()), $count , "get count");


is (($Entry_coord->coord()), $coord , "get coord");
is (($Entry_coord->name()), $coord , "get name coord");
is (($Entry_coord->chr1()), $chr1 , "get chr1 coord");
is (($Entry_coord->strand1()), $strand1 , "get strand1 coord");
is (($Entry_coord->pos1_start()), $pos1_start , "get pos1_start coord");
is (($Entry_coord->pos1_end()), $pos1_end , "get pos1_end coord");
is (($Entry_coord->chr2()), $chr2 , "get chr2 coord");
is (($Entry_coord->strand2()), $strand2 , "get strand2 coord");
is (($Entry_coord->pos2_start()), $pos2_start , "get pos2_start coord");
is (($Entry_coord->pos2_end()), $pos2_end , "get pos2_end coord");
is (($Entry_coord->shard()), $shard , "get shard coord");
is (($Entry_coord->count()), undef , "get count coord");
