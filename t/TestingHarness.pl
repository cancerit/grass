#!/usr/bin/perl
# Author: las (copied from wholegenomepipeline)

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use TAP::Harness;
use Cwd 'abs_path';

## locate the absolute path to the root of the project
my $abs_exe_path = abs_path($0);
my $abs_project_root_path = $abs_exe_path;
$abs_project_root_path =~ s/\/TestingHarness\.pl$//;

## setup the environment to point to local code base

my $adapt_path = $abs_project_root_path.'/../lib';
unshift @INC, $adapt_path;

## locate all the test files in the project.
my $dh;
opendir($dh,$abs_project_root_path) or die "|$abs_project_root_path| $!";
my @test_files = map{"$abs_project_root_path/$_"} grep {/\.t$/} readdir($dh);
close $dh;

## we need to send the lib locations of other cgp modules to the harness.
my $args = {
        verbosity => 1,
        lib => \@INC
};


my $harness = new TAP::Harness($args);
$harness->runtests(@test_files);
