package Sanger::CGP::Grass::DataEntry;

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


use strict;
use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

=head1 NAME

DataEntry

=head1 SYNOPSIS

use Sanger::CGP::Grass::DataEntry;

my $Entry = new Sanger::CGP::Grass::DataEntry(-name       => $name,
				 -chr1       => $chr1,
				 -strand1    => $strand1,
				 -pos1_start => $pos1_start,
				 -pos1_end   => $pos1_end,
				 -chr2       => $chr2,
				 -strand2    => $strand2,
				 -pos2_start => $pos2_start,
				 -pos2_end   => $pos2_end,
				 -shard      => $shard,
				 -count      => $count );

=head1 DESCRIPTION

    A holder for a pair of breakpoint entries that constitute a potential rearrangement

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::DataEntry();
  Description: make a new DataEntry object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless $self,$class;

    $self->{debug} = 0;

    if ($args{-name}) { $self->name($args{-name}) }
    if ($args{-coord}) { $self->coord($args{-coord}) }
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

    $self->name; # get it to generate a name/coord string if not provided
    if ($self->{debug}) { print "Use coord: " . $self->coord . "\n"; }

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
    Example    : $name = $object->name($name);
    Description: name of the entry
    Return     : name

=cut

sub name {
    my $self = shift;
    $self->{name} = shift if @_;
    if ($self->{chr1} && !($self->{name})) { $self->_create_name(); }
    return($self->{name});
}
#-----------------------------------------------------------------------#

=head2 coord

    Arg (1)    : coord
    Example    : $coord = $object->coord($coord);
    Description: coord of the entry
    Return     : coord

=cut

sub coord {
    my $self = shift;
    my $coord =  shift if @_;
    if ($coord) {
	if ($self->_parse_coords($coord)) { $self->{coord} = $coord; }
	else {
	    print "ERROR parsing coordinate $coord\n";
	    return;
	}
    }
    return($self->{coord});
}
#-----------------------------------------------------------------------#

=head2 chr1

    Arg (1)    : chr1
    Example    : $chr1 = $object->chr1($chr1);
    Description: breakpoint 1 chromosome
    Return     : chr1

=cut

sub chr1 {
    my $self = shift;
    $self->{chr1} = shift if @_;
    return($self->{chr1});
}
#-----------------------------------------------------------------------#

=head2 strand1

    Arg (1)    : strand1
    Example    : $strand1 = $object->strand1($strand1);
    Description: breakpoint 1 strand
    Return     : strand1

=cut

sub strand1 {
    my $self = shift;
    $self->{strand1} = shift if @_;
    return($self->{strand1});
}
#-----------------------------------------------------------------------#

=head2 pos1_start

    Arg (1)    : pos1_start
    Example    : $pos1_start = $object->pos1_start($pos1_start);
    Description: breakpoint 1 - start of location range
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
    Example    : $pos1_end = $object->pos1_end($pos1_end);
    Description: breakpoint 1 - end of location range
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
    Example    : $chr2 = $object->chr2($chr2);
    Description: breakpoint 2 chromosome
    Return     : chr2

=cut

sub chr2 {
    my $self = shift;
    $self->{chr2} = shift if @_;
    return($self->{chr2});
}
#-----------------------------------------------------------------------#

=head2  strand2

    Arg (1)    : strand2
    Example    : $strand2 = $object->strand2($strand2);
    Description: breakpoint 2 strand
    Return     : strand2

=cut

sub strand2 {
    my $self = shift;
    $self->{strand2} = shift if @_;
    return($self->{strand2});
}
#-----------------------------------------------------------------------#

=head2  pos2_start

    Arg (2)    : pos2_start
    Example    : $pos2_start = $object->pos2_start($pos2_start);
    Description: breakpoint 2 - start of location range
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
    Example    : $pos2_end = $object->pos2_end($pos2_end);
    Description: breakpoint 2 - end of location range
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
    Example    : $shard = $object->shard($shard);
    Description: shard sequence (Non-templated sequence)
                 If sequence is mappable, should treat as 2 breakpoint pairs
    Return     : shard

=cut

sub shard {
    my $self = shift;
    $self->{shard} = shift if @_;
    if ($self->{shard}) { $self->{shard} =~ s/[ \.]+//g; }
    return($self->{shard});
}
#-----------------------------------------------------------------------#

=head2 count

    Arg (1)    : count
    Example    : $count = $object->count($count);
    Description: number of read pairs that contributed to this group
    Return     : count

=cut

sub count {
    my $self = shift;
    $self->{count} = shift if @_;
    return($self->{count});
}
#-----------------------------------------------------------------------#
sub _create_name {
    my $self = shift;

    $self->{name} =  $self->{chr1} . ':' . $self->{strand1} . ':' . $self->{pos1_start} . '-' . $self->{pos1_end} . ',' . $self->{chr2} . ':' . $self->{strand2} . ':' . $self->{pos2_start} . '-' . $self->{pos2_end};

    if ($self->{shard}) { $self->{name} .= ',' . $self->{shard}; }

    unless ($self->{coord}) { $self->{coord} = $self->{name}; }
}
#-----------------------------------------------------------------------#

sub _parse_coords {
    my ($self, $coord_string) = @_;

    my ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard);
    my ($coord1, $coord2);

    if ($coord_string =~ /;/) {
	# parse coordinate - formats 1-:123-456;2+:345-678 or 1-:123;2+:345
	($coord1, $coord2, $shard) = split ';', $coord_string;
	return(0) unless ($coord1 && $coord2);

	my ($chr_strand1,$pos_string1) = split ':', $coord1;
	my ($chr_strand2,$pos_string2) = split ':', $coord2;
	($chr1, $strand1) = ($chr_strand1 =~ /(\S+?)(\S)$/);
	($chr2, $strand2) = ($chr_strand2 =~ /(\S+?)(\S)$/);
	return(0) unless ($chr1 && $strand1 && $pos_string1 && $chr2 && $strand2 && $pos_string2);

	$pos_string1 =~ s/[,\s]//g; # take out commas and spaces in the coordinates
	$pos1_start = $pos_string1;
	$pos1_end = $pos_string1;
	if ($pos_string1 =~ /(\d+)-(\d+)/) { $pos1_start = $1; $pos1_end = $2; }

	$pos_string2 =~ s/[,\s]//g;
	$pos2_start = $pos_string2;
	$pos2_end = $pos_string2;
	if ($pos_string2 =~ /(\d+)-(\d+)/) { $pos2_start = $1; $pos2_end = $2; }

	unless ($chr1 && $strand1 && $pos1_start && $pos1_end && $chr2 && $strand2 && $pos2_start && $pos2_end) {
	    print "$coord_string not parsed correctly. skip\n";
	    return(0);
	}
    }
    elsif ($coord_string =~ /,/) {
	# parse coordinate - formats 1:-:123-456,2:+:345-678 or 1:-:123,2:+:345
	($coord1, $coord2, $shard) = split ',', $coord_string;

	my ($pos_string1,$pos_string2);
	($chr1,$strand1,$pos_string1) = split ':', $coord1;
	($chr2,$strand2,$pos_string2) = split ':', $coord2;

	$pos_string1 =~ s/[,\s]//g; # take out commas and spaces in the coordinates
	$pos1_start = $pos_string1;
	$pos1_end = $pos_string1;
	if ($pos_string1 =~ /(\d+)-(\d+)/) { $pos1_start = $1; $pos1_end = $2; }

	$pos_string2 =~ s/[,\s]//g;
	$pos2_start = $pos_string2;
	$pos2_end = $pos_string2;
	if ($pos_string2 =~ /(\d+)-(\d+)/) { $pos2_start = $1; $pos2_end = $2; }
    }

    if ($shard) {
	$shard =~ s/ //g;
	unless ($shard =~ /^[aAtTgGcCnN]+$/) { print "WARN: shard sequence $shard is not valid. Can not use it\n"; $shard = undef; }
    }

    $self->{chr1} = $chr1;
    $self->{strand1} = $strand1;
    $self->{pos1_start} = $pos1_start;
    $self->{pos1_end} = $pos1_end;
    $self->{chr2} = $chr2;
    $self->{strand2} = $strand2;
    $self->{pos2_start} = $pos2_start;
    $self->{pos2_end} = $pos2_end;
    $self->{shard} = $shard;

    return(1);
}
#------------------------------------------------------------------------------------------------#

1;
