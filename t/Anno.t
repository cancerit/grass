#!/usr/bin/perl

# for testing Sanger::CGP::Grass::Anno class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Grass::Anno;
use Sanger::CGP::Grass::AnnoPoint;

use Test::More 'no_plan';

# existing entry

my $id_rg = 1234;
my $id_fusion_flag = 720;
my $Ltype = 'intron';
my $Htype = 'exon';
my $Llength = 100;
my $Hlength = 150;

my $apl5 = new Sanger::CGP::Grass::AnnoPoint();
my $apl3 = new Sanger::CGP::Grass::AnnoPoint();
my $aph5 = new Sanger::CGP::Grass::AnnoPoint();
my $aph3 = new Sanger::CGP::Grass::AnnoPoint();


# make a new object
my $Anno = new Sanger::CGP::Grass::Anno(-id_rg => $id_rg,
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
