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


# for testing grass.pl script

use Time::localtime;
use File::Copy qw(copy);
use Text::Diff;
use Capture::Tiny qw(capture);
use File::Temp qw(tempdir);

use strict;
use warnings FATAL => 'all';

use Test::More 'no_plan';
use FindBin qw($Bin);

# got to do something about this...
#my $genome_cache = '/lustre/scratch104/sanger/am3/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz';
my $genome_cache = "$Bin/../testData/vagrent.cache.gz";
my $ref = "$Bin/../testData/genome.fa";
if(exists $ENV{GRASS_GRCH37_FA}) {
  $ref = $ENV{GRASS_GRCH37_FA};
}

my $script = $Bin.'/../bin/' . 'grass.pl';
my $tmpdir = tempdir( CLEANUP => 1 );
note $tmpdir;

test_coord_input();
test_file_input();
test_file_input_bedpe();

#----------------------------------------------------------------------------------------#
sub test_coord_input {

    my $test_coord = 'X:-:84500101-84500632,X:+:84510265-84510265';
    my $test_coord_res = 'gene	gene_id	transcript_id	strand	end_phase	region	region_number	total_region_count	first/last	gene	gene_id	transcript_id	strand	phase	region	region_number	total_region_count	first/last	fusion_flag
ZNF711	ZNF711	ENST00000373165	1	_	5UTRintron	1	9	_	ZNF711	ZNF711	ENST00000373165	1	1	exon	4	9	first_base	570
';
    my $command = "$^X $script -genome_cache $genome_cache -coord $test_coord -ref $ref";
    note $command;
    my ($test_coord_out, $err, @result) = capture { system($command); };
    die "$err\n\t" if ($result[0]);
    is ($test_coord_out, $test_coord_res, "check supply coord");

}
#----------------------------------------------------------------------------------------#
sub test_file_input {

    my $infile = "$Bin/../testData/testout_Brass_I_and_II";
    my $outfile = "$Bin/../testData/testout_Brass_I_and_II_ann";


    my $testfile = "$tmpdir/testout_Brass_I_and_II";
    my $testfile_out = "$tmpdir/testout_Brass_I_and_II_ann";

    ok(copy($infile, $testfile), 'setup files');
    my $command = "$^X $script -genome_cache $genome_cache -file $testfile";
    note $command;
    my ($out, $err, @result) = capture { system($command); };
    die "$err\n\t" if ($result[0]);

    my $tf_out = slurp_file($testfile_out);
    my $of = slurp_file($outfile);
#use Data::Dumper;
#warn Dumper($tf_out);
#warn Dumper($of);
    is_deeply($tf_out, $of, 'correct file created');
}
#----------------------------------------------------------------------------------------#
sub test_file_input_bedpe {

    my $infile = "$Bin/../testData/testout_Brass.bedpe";
    my $outfile = "$Bin/../testData/testout_Brass_ann.bedpe";
    my $outfile_vcf = "$Bin/../testData/testout_Brass_ann.vcf";

    my $testfile = "$tmpdir/testout_Brass.bedpe";
    my $testfile_out = "$tmpdir/testout_Brass_ann.bedpe";
    my $testfile_out_vcf = "$tmpdir/testout_Brass_ann.vcf";

    ok(copy($infile, $testfile), 'setup files');
    unless (-e $ref) { die "ERROR: ref file $ref does not exist. exiting\n"; }
    my $command = "$^X $script -genome_cache $genome_cache -file $testfile -ref $ref";
    note "$command\n";
    my ($out, $err, @result) = capture { system($command); };
    die "$err\n\t" if ($result[0]);

    # the current file date is put in the output vcf file so it will always differ to the 'correct' file in testData unless updated
    my $date_in_correct_testfile = '20140512';
    my $tm = localtime;
    my $today = ($tm->year() + 1900) . sprintf("%02d",($tm->mon() + 1)) . sprintf("%02d",($tm->mday()));
    # do an inline substitution
    my $command1 = "$^X -pi -e \'s/##fileDate=$today/##fileDate=$date_in_correct_testfile/g\' $testfile_out_vcf";
    my $command2 = "$^X -pi -e \'s/##source_$today/##source_$date_in_correct_testfile/g\' $testfile_out_vcf";
    my $command3 = "$^X -pi -e \'s/##vcfProcessLog_$today/##vcfProcessLog_$date_in_correct_testfile/g\' $testfile_out_vcf";
    `$command1`;
    `$command2`;
    `$command3`;

    my $slurp_testfile_out = slurp_file($testfile_out);
    my $slurp_outfile = slurp_file($outfile);
    is_deeply($slurp_testfile_out, $slurp_outfile, 'correct bedpe file created');

    my $slurp_testfile_out_vcf = slurp_file($testfile_out_vcf);
    my $slurp_outfile_vcf = slurp_file($outfile_vcf);

    # path of reference file will differ so remove from comparison
    splice(@{$slurp_testfile_out_vcf}, 3,1);
    splice(@{$slurp_outfile_vcf}, 3,1);

    is_deeply($slurp_testfile_out_vcf, $slurp_outfile_vcf, 'correct vcf file created');
}
#----------------------------------------------------------------------------------------#

sub slurp_file {
  my $in = shift;
  open my $fh, '<', $in or die $!;
  my @data = <$fh>;
  close $fh;
  return \@data;
}
