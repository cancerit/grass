package Sanger::CGP::Grass::AnnoPoint;

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

AnnoPoint

=head1 SYNOPSIS

use Sanger::CGP::Grass::AnnoPoint;

my $Anno = new Sanger::CGP::Grass::AnnoPoint(-gene => $gene,
				-gene_id => $gene_id,
				-biotype => $biotype,
				-transcript => $transcript,
				-transcript_id => $transcript_id,
				-entrez_id => $entrez_id,
				-strand => $strand,
				-phase => $phase,
				-up2 => $up2,
				-down2 => $down2,
				-region => $region,
				-region_number => $region_number,
				-cds_pos => $cds_pos,
				-trans_pos => $trans_pos,
				-start_base => $start_base,
				-end_base => $end_base,
				-trans_region_count => $trans_region_count,
				-trans_length => $trans_length);

=head1 DESCRIPTION

Class holding individual annotation points (L5, L3, H5, H3)

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::AnnoPoint();
  Description: make a new AnnoPoint object
  Return     : AnnoPoint object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    # set the input parameters

    if ($args{-gene}) { $self->gene($args{-gene}); }
    if ($args{-gene_id}) { $self->gene_id($args{-gene_id}); }
    if ($args{-biotype}) { $self->biotype($args{-biotype}); }
    if ($args{-transcript}) { $self->transcript($args{-transcript}); }
    if ($args{-transcript_id}) { $self->transcript_id($args{-transcript_id}); }
    if ($args{-entrez_id}) { $self->entrez_id($args{-entrez_id}); }
    if ($args{-strand}) { $self->strand($args{-strand}); }
    if (defined($args{-phase})) { $self->phase($args{-phase}); }
    if ($args{-up2}) { $self->up2($args{-up2}); }
    if ($args{-down2}) { $self->down2($args{-down2}); }
    if ($args{-region}) { $self->region($args{-region}); }
    if ($args{-region_number}) { $self->region_number($args{-region_number}); }
    if ($args{-cds_pos}) { $self->cds_pos($args{-cds_pos}); }
    if ($args{-trans_pos}) { $self->trans_pos($args{-trans_pos}); }
    if ($args{-start_base}) { $self->start_base($args{-start_base}); }
    if ($args{-end_base}) { $self->end_base($args{-end_base}); }
    if ($args{-trans_region_count}) { $self->trans_region_count($args{-trans_region_count}); }
    if ($args{-trans_length}) { $self->trans_length($args{-trans_length}); }
    if ($args{-translation_length}) { $self->translation_length($args{-translation_length}); }

    return $self;
}

##-----------------------------------------------------------------------#

=head2 gene

  Arg (1)    : $gene
  Example    : $gene = $object->gene($gene);
  Description: gene name
  Return     : $gene

=cut

sub gene {
    my $self = shift;
    $self->{gene} = shift if @_;
    return $self->{gene};
}
#-----------------------------------------------------------------------#

=head2 gene_id

  Arg (1)    : $gene_id
  Example    : $gene_id = $object->gene_id($gene_id);
  Description: gene id
  Return     : $gene_id

=cut

sub gene_id {
    my $self = shift;
    $self->{gene_id} = shift if @_;
    return $self->{gene_id};
}
#-----------------------------------------------------------------------#

=head2 biotype

  Arg (1)    : $biotype
  Example    : $biotype = $object->biotype($biotype);
  Description: gene id
  Return     : $biotype

=cut

sub biotype {
    my $self = shift;
    $self->{biotype} = shift if @_;
    return $self->{biotype};
}
#-----------------------------------------------------------------------#

=head2 transcript

  Arg (1)    : $transcript
  Example    : $transcript = $object->transcript($transcript);
  Description: transcript name
  Return     : $transcript

=cut

sub transcript {
    my $self = shift;
    $self->{transcript} = shift if @_;
    return $self->{transcript};
}
#-----------------------------------------------------------------------#

=head2 transcript_id

  Arg (1)    : $transcript_id
  Example    : $transcript_id = $object->transcript_id($transcript_id);
  Description: transcript id
  Return     : $transcript_id

=cut

sub transcript_id {
    my $self = shift;
    $self->{transcript_id} = shift if @_;
    return $self->{transcript_id};
}
##-----------------------------------------------------------------------#

=head2 entrez_id

  Arg (1)    : $entrez_id
  Example    : $entrez_id = $object->entrez_id($entrez_id);
  Description: entrez id
  Return     : $entrez_id

=cut

sub entrez_id {
    my $self = shift;
    $self->{entrez_id} = shift if @_;
    return $self->{entrez_id};
}
#-----------------------------------------------------------------------#

=head2 strand

  Arg (1)    : $strand
  Example    : $strand = $object->strand($strand);
  Description: strand
  Return     : $strand

=cut

sub strand {
    my $self = shift;
    $self->{strand} = shift if @_;
    return $self->{strand};
}
#-----------------------------------------------------------------------#

=head2 phase

  Arg (1)    : $phase
  Example    : $phase = $object->phase($phase);
  Description: phase
  Return     : $phase

=cut

sub phase {
    my $self = shift;
    my $phase = shift if @_;
    if ((defined($phase)) && (($phase eq '0') || ($phase eq '1') || ($phase eq '2') || ($phase eq '-1'))) {
	$self->{phase} = $phase;
    }
    return $self->{phase};
}
#-----------------------------------------------------------------------#

=head2 up2

  Arg (1)    : $up2
  Example    : $up2 = $object->up2($up2);
  Description: 2 coding bases upstream of this breakpoint (to work out if stops are formed in fusion)
  Return     : $up2

=cut

sub up2 {
    my $self = shift;
    $self->{up2} = shift if @_;
    return $self->{up2};
}
#-----------------------------------------------------------------------#

=head2 down2

  Arg (1)    : $down2
  Example    : $down2 = $object->down2($down2);
  Description: 2 coding bases downstream of this breakpoint (to work out if stops are formed in fusion)
  Return     : $down2

=cut

sub down2 {
    my $self = shift;
    $self->{down2} = shift if @_;
    return $self->{down2};
}
#-----------------------------------------------------------------------#

=head2 region

  Arg (1)    : $region
  Example    : $region = $object->region($region);
  Description: region - whether its upstream, downstream, UTR5, UTR3, exon, intron or coding
  Return     : $region

=cut

sub region {
    my $self = shift;
    $self->{region} = shift if @_;
    return $self->{region};
}
#-----------------------------------------------------------------------#

=head2 region_number

  Arg (1)    : $region_number
  Example    : $region_number = $object->region_number($region_number);
  Description: region number. Which number exon or intron this is counting from the start of the transcript
  Return     : $region_number

=cut

sub region_number {
    my $self = shift;
    $self->{region_number} = shift if @_;
    return $self->{region_number};
}
#-----------------------------------------------------------------------#

=head2 cds_pos

  Arg (1)    : $cds_pos
  Example    : $cds_pos = $object->cds_pos($cds_pos);
  Description: cds position
  Return     : $cds_pos

=cut

sub cds_pos {
    my $self = shift;
    $self->{cds_pos} = shift if @_;
    return $self->{cds_pos};
}
#-----------------------------------------------------------------------#

=head2 trans_pos

  Arg (1)    : $trans_pos
  Example    : $trans_pos = $object->trans_pos($trans_pos);
  Description: transcript position
  Return     : $trans_pos

=cut

sub trans_pos {
    my $self = shift;
    $self->{trans_pos} = shift if @_;
    return $self->{trans_pos};
}
#-----------------------------------------------------------------------#

=head2 start_base

  Arg (1)    : $start_base
  Example    : $start_base = $object->start_base($start_base);
  Description: For exact breakpoints, whether this brakpoint is the first base of an exon
  Return     : $start_base

=cut

sub start_base {
    my $self = shift;
    $self->{start_base} = shift if @_;
    return $self->{start_base};
}
#-----------------------------------------------------------------------#

=head2 end_base

  Arg (1)    : $end_base
  Example    : $end_base = $object->end_base($end_base);
  Description: For exact breakpoints, whether this breakpoint is the last base of an exon
  Return     : $end_base

=cut

sub end_base {
    my $self = shift;
    $self->{end_base} = shift if @_;
    return $self->{end_base};
}
#-----------------------------------------------------------------------#

=head2 trans_region_count

  Arg (1)    : $trans_region_count
  Example    : $trans_region_count = $object->trans_region_count($trans_region_count);
  Description: number of exons in this transcript
  Return     : $trans_region_count

=cut

sub trans_region_count {
    my $self = shift;
    $self->{trans_region_count} = shift if @_;
    return $self->{trans_region_count};
}
#-----------------------------------------------------------------------#

=head2 trans_length

  Arg (1)    : $trans_length
  Example    : $trans_length = $object->trans_length($trans_length);
  Description: length of this transcript
  Return     : $trans_length

=cut

sub trans_length {
    my $self = shift;
    $self->{trans_length} = shift if @_;
    return $self->{trans_length};
}
#-----------------------------------------------------------------------#

=head2 translation_length

  Arg (1)    : $translation_length
  Example    : $translation_length = $object->translation_length($translation_length);
  Description: length of this transcript's translation product
  Return     : $translation_length

=cut

sub translation_length {
    my $self = shift;
    $self->{translation_length} = shift if @_;
    return $self->{translation_length};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

1;
