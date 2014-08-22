#!/usr/bin/perl

use ExtUtils::MakeMaker;

WriteMakefile(
  NAME          => 'Grass',
  VERSION_FROM  => 'lib/Sanger/CGP/Grass.pm',
  EXE_FILES     => [qw( bin/grass.pl
                        bin/brass_bedpe2vcf.pl)],
  PREREQ_PM     => {
                     'Const::Fast' => 0.014,
                     'Try::Tiny' => 0.19,
                     'List::Util' => 1.38,
                     'Devel::Cover' => 1.09,
                     'Pod::Coverage' => 0.23,
                     'Attribute::Abstract' => 1.07,
                     'Sub::Exporter::Progressive' => 0.001011,
                     'Bio::DB::Sam' => 1.39,
                     'Bio::Root::Version' => 1.006923, # this is the only way to get Bio::Perl to behave
                                                       # and return the version correctly
                   }
);
