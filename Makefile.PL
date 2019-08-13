#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
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

use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share module => 'Sanger::CGP::Grass', 'share';

WriteMakefile(
  NAME          => 'Grass',
  LICENSE       => 'agpl_3', # http://search.cpan.org/~dagolden/CPAN-Meta-2.142690/lib/CPAN/Meta/Spec.pm#license
  VERSION_FROM  => 'lib/Sanger/CGP/Grass.pm',
  EXE_FILES     => [qw( bin/grass.pl
                        bin/brass_bedpe2vcf.pl)],
  PREREQ_PM     => {
                     'Const::Fast' => 0.014,
                     'Try::Tiny' => 0.19,
                     'List::Util' => 1.38,
                     'Pod::Usage' => 1.64,
                     'Attribute::Abstract' => 1.02,
                     'File::ShareDir' => 1.102,
                     'Bio::Root::Version' => 1.006924,# this is the only way to get Bio::Perl to behave
                                                       # and return the version correctly
                   }
);

package MY;
use File::ShareDir::Install qw( postamble );
