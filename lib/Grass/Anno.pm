## Grass::Anno

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
#
# Author las
#
=head1 NAME

Project

=head1 SYNOPSIS

=head1 DESCRIPTION

Class holding fusion analysis annotation results.
Lists annotations for a pair of coordinates/coordinate ranges.

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Grass::Anno;
use Grass::AnnoPoint;

use strict;

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    : 
  Example    : $object = new Grass::Anno();
  Description: make a new Anno object
  Return     : Anno object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    # set the input parameters

    if ($args{-id_rg})          { $self->id_rg($args{-id_rg}); }
    if ($args{-id_fusion_flag}) { $self->id_fusion_flag($args{-id_fusion_flag}); }
    if ($args{-L5})             { $self->L5($args{-L5}); } # anno point objects
    if ($args{-L3})             { $self->L3($args{-L3}); }
    if ($args{-H5})             { $self->H5($args{-H5}); }
    if ($args{-H3})             { $self->H3($args{-H3}); }
    if ($args{-Ltype})          { $self->Ltype($args{-Ltype}); } # type of annotation over this end
    if ($args{-Htype})          { $self->Htype($args{-Htype}); }
    if ($args{-Llength})        { $self->Llength($args{-Llength}); }
    if ($args{-Hlength})        { $self->Hlength($args{-Hlength}); }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 id_rg

  Arg (0)    :
  Example    : $id_rg = $pair->id_rg($id_rg);
  Description: id_rearrangement_group
  Return     : $id_rg

=cut

sub id_rg {
    my $self = shift;
    $self->{id_rg} = shift if @_;
    return $self->{id_rg};
}
#-----------------------------------------------------------------------#

=head2 id_fusion_flag

  Arg (0)    :
  Example    : $id_fusion_flag = $pair->id_fusion_flag($id_fusion_flag);
  Description: id_fusion_flag
  Return     : $id_fusion_flag

=cut

sub id_fusion_flag {
    my $self = shift;
    $self->{id_fusion_flag} = shift if @_;
    return $self->{id_fusion_flag};
}
#-----------------------------------------------------------------------#

=head2 L5

  Arg (0)    :
  Example    : $L5 = $pair->L5($L5);
  Description: object representing the 5 prime end of the Lower end of the rearrangement (L5) 
  Return     : $L5

=cut

sub L5 {
    my $self = shift;
    $self->{L5} = shift if @_;
    return $self->{L5};
}
#-----------------------------------------------------------------------#

=head2 L3

  Arg (0)    :
  Example    : $L3 = $pair->L3($L3);
  Description: object representing the 3 prime end of the Lower end of the rearrangement (L3) 
  Return     : $L3

=cut

sub L3 {
    my $self = shift;
    $self->{L3} = shift if @_;
    return $self->{L3};
}
#-----------------------------------------------------------------------#

=head2 Ltype

  Arg (0)    :
  Example    : $Ltype = $pair->Ltype($Ltype);
  Description: type of the lower end of the rearrangement
  Return     : $Ltype

=cut

sub Ltype {
    my $self = shift;
    $self->{Ltype} = shift if @_;
    return $self->{Ltype};
}
#-----------------------------------------------------------------------#

=head2 H5

  Arg (0)    :
  Example    : $H5 = $pair->H5($H5);
  Description: object representing the 5 prime end of the Higher end of the rearrangement (H5) 
  Return     : $H5

=cut

sub H5 {
    my $self = shift;
    $self->{H5} = shift if @_;
    return $self->{H5};
}
#-----------------------------------------------------------------------#

=head2 H3

  Arg (0)    :
  Example    : $H3 = $pair->H3($H3);
  Description: object representing the 3 prime end of the Higher end of the rearrangement (H3) 
  Return     : $H3

=cut

sub H3 {
    my $self = shift;
    $self->{H3} = shift if @_;
    return $self->{H3};
}
#-----------------------------------------------------------------------#

=head2 Htype

  Arg (0)    :
  Example    : $Htype = $pair->Htype($Htype);
  Description: type of the higher end of the rearrangement
  Return     : $Ltype

=cut

sub Htype {
    my $self = shift;
    $self->{Htype} = shift if @_;
    return $self->{Htype};
}
#-----------------------------------------------------------------------#

=head2 Llength

  Arg (0)    :
  Example    : $Llength = $pair->Llength($Llength);
  Description: length of cds to L side of breakpoint
  Return     : $Llength

=cut

sub Llength {
    my $self = shift;
    $self->{Llength} = shift if @_;
    return $self->{Llength};
}
#-----------------------------------------------------------------------#

=head2 Hlength

  Arg (0)    :
  Example    : $H_length = $pair->Hlength($Hlength);
  Description: length of cds to H side of breakpoint
  Return     : $Hlength

=cut

sub Hlength {
    my $self = shift;
    $self->{Hlength} = shift if @_;
    return $self->{Hlength};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

1;
