#!/usr/bin/perl -w

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


use Bio::DB::Sam; 

use strict;
use Getopt::Long;

my $ref = '/nfs/cancer_ref01/human/37/genome.fa';
my $coord = '';
my $help = 0;

GetOptions( 'ref:s' => \$ref,
            'coord:s'  => \$coord,
            'help:s'   => \$help);

# check inputs
if ($help) { usage(); }

# load an indexed fasta file
print "ref $ref\n";
my $fai = Bio::DB::Sam::Fai->load( $ref );

if ($coord =~ /^(\S+?):(\d+)$/) {
    my $chr = $1;
    my $pos = $2;
    $coord = $chr . ':' . $pos . '-' . $pos;
}

my $seq = $fai->fetch("$coord");
print "$seq\n";

#------------------------------------------------------------------------------------------------------------#
sub usage {

    print <<HERE;

get_fai_seq.pl

Description - gets a sequence from a fasta reference file (and associated fai index file) for supplied coordinates


options...

   -coord         : coordinate string format chr:start-end eg 1:123-456
   -ref           : fasta reference file (with associated fai file). For vcf out file generation.
   -help          : Print this message

examples...

get_fai_seq.pl -ref /nfs/cancer_ref01/human/37/genome.fa -coord 7:12350-12400 



Author : las


HERE

exit(0);
}
