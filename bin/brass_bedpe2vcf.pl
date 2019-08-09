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


# takes brassI output
# gets flanking bases
# converts to vcf

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use Sanger::CGP::Grass::FlankingBases;
use Sanger::CGP::Grass::VcfContigs;
use Sanger::CGP::Grass::VcfConverter;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Vcf::VcfProcessLog;

use Carp;
use Pod::Usage qw(pod2usage);
use Data::Dumper;
use Getopt::Long;

$| = 1;

my $file = '';
my $file2 = '';
my $outfile = '';
my $is_bedpe = 0;
my $list_between = 0;
my $show_biotype = 0;
my $use_all_biotypes = 0;
my $gene_id_required = 0;

# if using ensembl api and for vcf generation
my $species = ''; # eg HUMAN

# for generating vcf format
my $ref = '';
my $assembly = '';
my $platform = '';
my $protocol = '';
my $tumour = '';
my $acc_tumour = '';
my $acc_source_tumour = '';
my $study_tumour = '';
my $normal = '';
my $acc_normal = '';
my $acc_source_normal = '';
my $study_normal = '';

my $help = 0;

GetOptions( 'file:s'        => \$file,
	    'species:s'     => \$species,
	    'ref:s'         => \$ref,
	    'assembly:s'    => \$assembly,
	    'platform:s'    => \$platform,
	    'protocol:s'    => \$protocol,
        'version'       => \$version,
	    'tumour:s'            => \$tumour,
	    'acc_tumour:s'        => \$acc_tumour,
	    'acc_source_tumour:s' => \$acc_source_tumour,
	    'study_tumour:s'      => \$study_tumour,
	    'normal:s'           => \$normal,
	    'acc_normal:s'       => \$acc_normal,
	    'acc_source_normal:s'=> \$acc_source_normal,
	    'study_normal:s'     => \$study_normal,
	    'help'          => \$help,
         );

# check inputs
if ($help) { usage(); }

if($version){
  print 'Version: '.Sanger::CGP::Grass->VERSION."\n";
  exit;
}

# get the flanking bases for use in vcf conversion (appends the up and downstream bases to the end of each line)
if ($file) {
    my $FlankingBases = new Sanger::CGP::Grass::FlankingBases(-infile => $file,
							      -ref    => $ref,
							      -flip_strand => 1);
    $FlankingBases->process();
}

make_vcf_file($file, $ref, $species, $assembly, $platform, $protocol, $tumour, $acc_tumour, $acc_source_tumour, $study_tumour, $normal, $acc_normal, $acc_source_normal, $study_normal);

#------------------------------------------------------------------------------------------------------------#
sub make_vcf_file {
    my ($file, $ref, $species, $assembly, $platform, $protocol, $tumour, $acc_tumour, $acc_source_tumour, $study_tumour, $normal, $acc_normal, $acc_source_normal, $study_normal) = @_;

    # could get out the bam header and get the contig and sample vcf objects
    #$header = `samtools view -H /nfs/cancer_trk0003/00000009/191035.v1.brm.bam`;

    # ...or just create the contig objects from the fai file and the wt and tumour sample objects here...

    my $contigs = [];
    if ($ref && $species && $assembly) {
	my $contig_o = new Sanger::CGP::Grass::VcfContigs(-fai => ($ref . '.fai'),
							  -species => $species,
							  -assembly => $assembly);
	$contigs = $contig_o->generate();
    }

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

    if ($file) {
	my $VcfConverter = new Sanger::CGP::Grass::VcfConverter(-infile  => $file,
								-contigs => $contigs,
								-flip_strand => 1);
	$VcfConverter->convert($wt_sample, $mt_sample, \@process_logs, $ref, $source);
    }
}
#----------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
sub usage {

    print <<HERE;

brass_bedpe2vcf.pl

Description - converts brass I and II bedpe files to vcf output

options...

   -species       : species (HUMAN, MOUSE, DOG etc)
   -file          : input file - format type: bedpe)
   -outfile       : what file to print output to. Default is input_file.bedpe
   -ref           : fasta reference file (with associated fai file). For vcf out file generation.
   -assembly      : sequence assembly used (eg GRCh37). For vcf out file generation.
   -platform      : sequencing platform used (eg HiSeq). For vcf out file generation.
   -protocol      : sequencing experimental design (eg genomic, pulldown). For vcf out file generation.
   -tumour         : name of tumour sample. For vcf out file generation.
   -acc_tumour     : name of tumour sample accession id. For vcf out file generation.
   -acc_source_tumour : source of tumour sample accession id. For vcf out file generation.
   -study_tumour   : study id associated with tumour sample. For vcf out file generation.
   -normal         : name of normal sample. For vcf out file generation.
   -acc_normal     : name of normal sample accession id. For vcf out file generation.
   -acc_source_normal : source of normal sample accession id. For vcf out file generation.
   -study_normal   : study id associated with normal sample. For vcf out file generation.
   -help           : Print this message

examples...

brass_bedpe2vcf.pl -species HUMAN -ref /nfs/cancer_ref01/human/37/genome.fa -assembly GRCh37 -platform HiSeq -protocol genomic -tumour HCC1395 -acc_tumour 1234 -acc_source_tumour COSMIC_SAMPLE_ID -study_tumour 111 -normal 1395BL -acc_normal 1235 -acc_source_normal COSMIC_SAMPLE_ID -study_normal 222 -file HCC1395_191535.v1.bedpe



Author : las


HERE

exit(0);
}
