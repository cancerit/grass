package Sanger::CGP::Grass::GenomeData::Exon;

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


use strict;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

=head1 NAME

Exon

=head1 SYNOPSIS

use Sanger::CGP::Grass::GenomeData::Exon;

my $grass_exon = new Sanger::CGP::Grass::GenomeData::Exon(-phase  => $exon_phase,
    -end_phase  => $exon_end_phase,
    -start  => $exon_start,
    -end    => $exon_end,
    -coding_region_start => $exon_cds_start,
    -coding_region_end   => $exon_cds_end,
    -seq    => $exon_seq );


=head1 DESCRIPTION

Class to hold exon details

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::GenomeData::Exon();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    if (defined($args{-phase}))       { $self->phase($args{-phase}); }
    if (defined($args{-end_phase}))   { $self->end_phase($args{-end_phase}); }
    if ($args{-start})                { $self->start($args{-start}); }
    if ($args{-end})                  { $self->end($args{-end}); }
    if ($args{-coding_region_start})  { $self->coding_region_start($args{-coding_region_start}); }
    if ($args{-coding_region_end})    { $self->coding_region_end($args{-coding_region_end}); }
    if ($args{-seq})                  { $self->seq($args{-seq}); }

    return $self;
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 phase

  Arg (1)    : $phase
  Example    : $phase = $Object->phase($phase);
  Description: which phase to get from Ensembl
  Return     : phase object

=cut

sub phase {
    my $self = shift;
    $self->{phase} = shift if @_;
    return $self->{phase};
}
#-----------------------------------------------------------------------#

=head2 end_phase

  Arg (1)    : $end_phase
  Example    : $end_phase = $Object->end_phase($end_phase);
  Description: which end_phase to get from Ensembl
  Return     : end_phase object

=cut

sub end_phase {
    my $self = shift;
    $self->{end_phase} = shift if @_;
    return $self->{end_phase};
}
#-----------------------------------------------------------------------#

=head2 start

  Arg (1)    : $start
  Example    : $start = $Object->start($start);
  Description: start
  Return     : start object

=cut

sub start {
    my $self = shift;
    $self->{start} = shift if @_;
    return $self->{start};
}
#-----------------------------------------------------------------------#

=head2 end

  Arg (1)    : $end
  Example    : $end = $Object->end($end);
  Description: end
  Return     : end object

=cut

sub end {
    my $self = shift;
    $self->{end} = shift if @_;
    return $self->{end};
}
#-----------------------------------------------------------------------#

=head2 coding_region_start

  Arg (1)    : $coding_region_start
  Example    : $coding_region_start = $Object->coding_region_start($coding_region_start);
  Description: coding_region_start
  Return     : coding_region_start object

=cut

sub coding_region_start {
    my $self = shift;
    $self->{coding_region_start} = shift if @_;
    return $self->{coding_region_start};
}
#-----------------------------------------------------------------------#

=head2 coding_region_end

  Arg (1)    : $coding_region_end
  Example    : $coding_region_end = $Object->coding_region_end($coding_region_end);
  Description: coding_region_end
  Return     : coding_region_end object

=cut

sub coding_region_end {
    my $self = shift;
    $self->{coding_region_end} = shift if @_;
    return $self->{coding_region_end};
}
#-----------------------------------------------------------------------#

=head2 seq

  Arg (1)    : $seq
  Example    : $seq = $Object->seq($seq);
  Description: seq
  Return     : seq object

=cut

sub seq {
    my $self = shift;
    $self->{seq} = shift if @_;
    return $self->{seq};
}
#-----------------------------------------------------------------------#
1;
