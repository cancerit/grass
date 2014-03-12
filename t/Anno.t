#!/usr/bin/perl

# for testing Grass::Anno class

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

use Grass::Anno;
use Grass::AnnoPoint;

use Test::More 'no_plan';

# existing entry

my $id_rg = 1234;
my $id_fusion_flag = 720;
my $Ltype = 'intron';
my $Htype = 'exon';
my $Llength = 100;
my $Hlength = 150;

my $apl5 = new Grass::AnnoPoint();
my $apl3 = new Grass::AnnoPoint();
my $aph5 = new Grass::AnnoPoint();
my $aph3 = new Grass::AnnoPoint();


# make a new object
my $Anno = new Grass::Anno(-id_rg => $id_rg,
                           -id_fusion_flag => $id_fusion_flag,
			   -L5 => $apl5,
			   -L3 => $apl3,
			   -H5 => $aph5,
			   -H3 => $aph3,
			   -Ltype => $Ltype,
			   -Htype => $Htype,
			   -Llength => $Llength,
			   -Hlength => $Hlength );

ok($Anno, 'object defined');

is (($Anno->id_rg()), $id_rg , "get id_rg");
is (($Anno->id_fusion_flag()), $id_fusion_flag , "get id_fusion_flag");
is (($Anno->L5()), $apl5 , "get L5");
is (($Anno->L3()), $apl3 , "get L3");
is (($Anno->H5()), $aph5 , "get H5");
is (($Anno->H3()), $aph3 , "get H3");
is (($Anno->Ltype()), $Ltype , "get Ltype");
is (($Anno->Htype()), $Htype , "get Htype");
is (($Anno->Llength()), $Llength , "get Llength");
is (($Anno->Hlength()), $Hlength , "get Hlength");
