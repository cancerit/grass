### Grass::DataEntry

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

Params

=head1 SYNOPSIS

=head1 DESCRIPTION

    A holder for a data entry

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Grass::DataEntry;

use strict;

#-----------------------------------------------------------------------#

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless $self,$class;

    if ($args{-name}) { $self->name($args{-name}) }
    if ($args{-chr1}) { $self->chr1($args{-chr1}) }
    if ($args{-strand1}) { $self->strand1($args{-strand1}) }
    if ($args{-pos1_start}) { $self->pos1_start($args{-pos1_start}) }
    if ($args{-pos1_end}) { $self->pos1_end($args{-pos1_end}) }
    if ($args{-chr2}) { $self->chr2($args{-chr2}) }
    if ($args{-strand2}) { $self->strand2($args{-strand2}) }
    if ($args{-pos2_start}) { $self->pos2_start($args{-pos2_start}) }
    if ($args{-pos2_end}) { $self->pos2_end($args{-pos2_end}) }
    if ($args{-shard}) { $self->shard($args{-shard}) }
    if ($args{-count}) { $self->count($args{-count}) }

    return $self;
}

#-----------------------------------------------------------------------#
######### code ############
#-----------------------------------------------------------------------#

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 name 
 
    Arg (1)    : name
    Example    : $name = $params->name($name);
    Description: name of the entry
    Return     : 

=cut

sub name {
    my $self = shift;
    $self->{name} = shift if @_;
    return($self->{name});
}
#-----------------------------------------------------------------------#
=head2 chr1 
 
    Arg (1)    : chr1
    Example    : $chr1 = $params->chr1($chr1);
    Description: read 1 chromosome
    Return     : 

=cut

sub chr1 {
    my $self = shift;
    $self->{chr1} = shift if @_;
    return($self->{chr1});
}
#-----------------------------------------------------------------------#
=head2 strand1
 
    Arg (1)    : strand1
    Example    : $strand1 = $params->strand1($strand1);
    Description: read 1 strand
    Return     : 

=cut

sub strand1 {
    my $self = shift;
    $self->{strand1} = shift if @_;
    return($self->{strand1});
}
#-----------------------------------------------------------------------#
=head2 pos1_start
 
    Arg (1)    : pos1_start
    Example    : $pos1_start = $params->pos1_start($pos1_start);
    Description: read 1 - start of location range
    Return     : pos1_start

=cut

sub pos1_start {
    my $self = shift;
    $self->{pos1_start} = shift if @_;
    return($self->{pos1_start});
}
#-----------------------------------------------------------------------#
=head2 pos1_end
 
    Arg (1)    : pos1_end
    Example    : $pos1_end = $params->pos1_end($pos1_end);
    Description: read 1 - end of location range
    Return     : pos1_end

=cut

sub pos1_end {
    my $self = shift;
    $self->{pos1_end} = shift if @_;
    return($self->{pos1_end});
}
#-----------------------------------------------------------------------#
=head2 chr2
 
    Arg (1)    : chr2
    Example    : $chr2 = $params->chr2($chr2);
    Description: read 2 chromosome
    Return     : 

=cut

sub chr2 {
    my $self = shift;
    $self->{chr2} = shift if @_;
    return($self->{chr2});
}
#-----------------------------------------------------------------------#
=head2  strand2
 
    Arg (1)    : strand2
    Example    : $strand2 = $params->strand2($strand2);
    Description: read 2 strand
    Return     : 

=cut

sub strand2 {
    my $self = shift;
    $self->{strand2} = shift if @_;
    return($self->{strand2});
}
#-----------------------------------------------------------------------#
=head2  pos2_start
 
    Arg (2)    : pos2_start
    Example    : $pos2_start = $params->pos2_start($pos2_start);
    Description: read 2 - start of location range
    Return     : pos2_start

=cut

sub pos2_start {
    my $self = shift;
    $self->{pos2_start} = shift if @_;
    return($self->{pos2_start});
}
#-----------------------------------------------------------------------#
=head2 pos2_end
 
    Arg (1)    : pos2_end
    Example    : $pos2_end = $params->pos2_end($pos2_end);
    Description: read 2 - end of location range
    Return     : pos2_end

=cut

sub pos2_end {
    my $self = shift;
    $self->{pos2_end} = shift if @_;
    return($self->{pos2_end});
}
#-----------------------------------------------------------------------#
=head2 shard
 
    Arg (1)    : shard
    Example    : $shard = $params->shard($shard);
    Description: shard sequence (Non-templated sequence)
                 If sequence is mappable, should treat as 2 breakpoints (link in cosmic)
    Return     : shard

=cut

sub shard {
    my $self = shift;
    $self->{shard} = shift if @_;
    if ($self->{shard}) { $self->{shard} =~ s/ //g; }
    return($self->{shard});
}
#-----------------------------------------------------------------------#
=head2 count
 
    Arg (1)    : count
    Example    : $count = $params->count($count);
    Description: number of read pairs that contributed to this group
    Return     : count

=cut

sub count {
    my $self = shift;
    $self->{count} = shift if @_;
    return($self->{count});
}
#-----------------------------------------------------------------------#

1;
