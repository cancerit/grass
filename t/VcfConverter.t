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


# for testing Sanger::CGP::Grass::VcfConverter class

use strict;
use warnings FATAL => 'all';

use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Vcf::VcfProcessLog;
use Sanger::CGP::Grass::VcfConverter;
use Sanger::CGP::Grass::VcfContigs;

use Time::localtime;
use File::Copy qw(copy);
use Text::Diff;
use Capture::Tiny qw(capture);
use File::Basename qw(basename);

use Test::More 'no_plan';
use FindBin qw($Bin);

my $ref = "$Bin/../testData/genome.fa";

test_brassI_bedpe_file_input($ref);
test_brassII_bedpe_file_input($ref);

#------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

sub test_brassI_bedpe_file_input {
    my ($ref) = @_;

    my $infile = $Bin.'/../testData/test_VcfConverter_ann.bedpe';
    my $outfile = $Bin.'/../testData/test_VcfConverter_ann.vcf';

    my $testfile = 'test_VcfConverter_ann.bedpe';
    my $testoutfile = 'test_VcfConverter_ann.vcf';

    copy $infile, $testfile;

    my ($wt_sample, $mt_sample, $contigs, $process_logs, $source) = get_support_objects($ref);

    # make a new object
    my $VcfConverter = new Sanger::CGP::Grass::VcfConverter(-infile  => $testfile,
                                                            -contigs => $contigs );
    ok($VcfConverter, 'object defined');

    $VcfConverter->convert($wt_sample, $mt_sample, $process_logs, $ref, $source);

    # the current file date is put in the output vcf file so it will always differ to the 'correct' file in testData unless updated
    my $date_in_correct_testfile = '20140512';
    my $tm = localtime;
    my $today = ($tm->year() + 1900) . sprintf("%02d",($tm->mon() + 1)) . sprintf("%02d",($tm->mday()));
    # do an inline substitution
    my $command1 = "perl -pi -e \'s/##fileDate=$today/##fileDate=$date_in_correct_testfile/g\' $testoutfile";
    my $command2 = "perl -pi -e \'s/##source_$today/##source_$date_in_correct_testfile/g\' $testoutfile";
    my $command3 = "perl -pi -e \'s/##vcfProcessLog_$today/##vcfProcessLog_$date_in_correct_testfile/g\' $testoutfile";
    `$command1`;
    `$command2`;
    `$command3`;

    my $testoutdata = slurp_file($testoutfile);
    my $outdata = slurp_file($outfile);

    # versions and path of reference file will differ so remove from comparison
    splice(@{$testoutdata}, 2,2);
    splice(@{$outdata}, 2,2);

    # this will move if more headers added
    splice(@{$testoutdata}, 59,1);
    splice(@{$outdata}, 59,1);

    is_deeply($testoutdata, $outdata, 'correct file created');

	unlink $testfile;
	unlink $testoutfile;
}
#----------------------------------------------------------------------------------------#

sub test_brassII_bedpe_file_input {
    my ($ref) = @_;

    my $infile = $Bin.'/../testData/test_VcfConverterII_ann.bedpe';
    my $outfile = $Bin.'/../testData/test_VcfConverterII_ann.vcf';

    my $testfile = 'test_VcfConverterII_ann.bedpe';
    my $testoutfile = 'test_VcfConverterII_ann.vcf';

    copy $infile, $testfile;

    my ($wt_sample, $mt_sample, $contigs, $process_logs, $source) = get_support_objects($ref);

    # make a new object
    my $VcfConverter = new Sanger::CGP::Grass::VcfConverter(-infile  => $testfile,
                                                            -contigs => $contigs );
    ok($VcfConverter, 'object defined');

    $VcfConverter->convert($wt_sample, $mt_sample, $process_logs, $ref, $source);

    # the current file date is put in the output vcf file so it will always differ to the 'correct' file in testData unless updated
    my $date_in_correct_testfile = '20140512';
    my $tm = localtime;
    my $today = ($tm->year() + 1900) . sprintf("%02d",($tm->mon() + 1)) . sprintf("%02d",($tm->mday()));
    # do an inline substitution
    my $command1 = "perl -pi -e \'s/##fileDate=$today/##fileDate=$date_in_correct_testfile/g\' $testoutfile";
    my $command2 = "perl -pi -e \'s/##source_$today/##source_$date_in_correct_testfile/g\' $testoutfile";
    my $command3 = "perl -pi -e \'s/##vcfProcessLog_$today/##vcfProcessLog_$date_in_correct_testfile/g\' $testoutfile";
    `$command1`;
    `$command2`;
    `$command3`;

    my $testoutdata = slurp_file($testoutfile);
    my $outdata = slurp_file($outfile);

    # versions and path of reference file will differ so remove from comparison
    splice(@{$testoutdata}, 2,2);
    splice(@{$outdata}, 2,2);

    # this will move if more headers added
    splice(@{$testoutdata}, 59,1);
    splice(@{$outdata}, 59,1);

    is_deeply($testoutdata, $outdata, 'correct file created');

	unlink $testfile;
	unlink $testoutfile;
}
#----------------------------------------------------------------------------------------#
sub get_support_objects {

    my ($ref) = @_;

    my $tumour = 'PD1234a';
    my $acc_tumour = 123;
    my $acc_source_tumour = 'ID_SAMPLE_COSMIC';
    my $study_tumour = 12345;
    my $normal = 'PD1234b';
    my $acc_normal = 45;
    my $acc_source_normal = 'ID_SAMPLE_COSMIC';
    my $study_normal = 12346;
    my $assembly = 'GRC37';
    my $species = 'HUMAN';
    my $platform = 'HiSeq';
    my $protocol = 'genomic';

    # could get out the bam header and get the contig and sample vcf objects
    #$header = `samtools view -H /nfs/cancer_trk0003/00000009/191035.v1.brm.bam`;

    # ...or just create the contig objects from the fai file and the wt and tumour sample objects here...

    my $contig_o = new Sanger::CGP::Grass::VcfContigs(-fai => ($ref . '.fai'),
						     -species => $species,
						     -assembly => $assembly);
    my $contigs = $contig_o->generate();

    my $mt_sample = new Sanger::CGP::Vcf::Sample( -name => $tumour,
						  -study => $study_tumour,
						  -platform => $platform,
						  -seq_protocol => $protocol,
						  -accession => $acc_tumour,
						  -accession_source => $acc_source_tumour,
						  -description => 'Mutant' );

    my $wt_sample = new Sanger::CGP::Vcf::Sample( -name => $normal,
						  -study => $study_normal,
						  -platform => $platform,
						  -seq_protocol => $protocol,
						  -accession => $acc_normal,
						  -accession_source => $acc_source_normal,
						  -description => 'Normal' );

    # make the process logs to put in the vcf header
    my $opts = join " ", $0, @ARGV;
    my @process_logs = ();
    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog( -input_vcf_source => basename($0),
							     -input_vcf_ver => Sanger::CGP::Grass->VERSION,
							     -input_vcf_param => $opts );

    my $source = basename($0). '_v'. Sanger::CGP::Grass->VERSION;

    return($wt_sample, $mt_sample, $contigs, \@process_logs, $source);
}
#----------------------------------------------------------------------------------------#


sub slurp_file {
  my $file = shift;
  open my $fh, '<', $file || die $!;
  my @lines = <$fh>;
  close $fh;
  return \@lines;
}
