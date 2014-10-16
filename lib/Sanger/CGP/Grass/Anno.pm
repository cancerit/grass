package Sanger::CGP::Grass::Anno;

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


use Sanger::CGP::Grass::AnnoPoint;

use strict;
use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

=head1 NAME

Anno

=head1 SYNOPSIS

use Sanger::CGP::Grass::Anno;

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

=head1 DESCRIPTION

Class holding grass fusion analysis annotation results.

Lists annotations for a pair of coordinates/coordinate ranges.

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::Anno();
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

  Arg (1)    : $id_rg
  Example    : $id_rg = $object->id_rg($id_rg);
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

  Arg (1)    : $id_fusion_flag
  Example    : $id_fusion_flag = $object->id_fusion_flag($id_fusion_flag);
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

  Arg (1)    : $L5
  Example    : $L5 = $object->L5($L5);
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

  Arg (1)    : $L3
  Example    : $L3 = $object->L3($L3);
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

  Arg (1)    : $Ltype
  Example    : $Ltype = $object->Ltype($Ltype);
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

  Arg (1)    : $H5
  Example    : $H5 = $object->H5($H5);
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

  Arg (1)    : $H3
  Example    : $H3 = $object->H3($H3);
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

  Arg (1)    : $Htype
  Example    : $Htype = $object->Htype($Htype);
  Description: type of the higher end of the rearrangement
  Return     : $Htype

=cut

sub Htype {
    my $self = shift;
    $self->{Htype} = shift if @_;
    return $self->{Htype};
}
#-----------------------------------------------------------------------#

=head2 Llength

  Arg (1)    : $Llength
  Example    : $Llength = $object->Llength($Llength);
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

  Arg (1)    : $Hlength
  Example    : $H_length = $object->Hlength($Hlength);
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
