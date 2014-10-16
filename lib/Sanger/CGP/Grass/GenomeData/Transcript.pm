package Sanger::CGP::Grass::GenomeData::Transcript;

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

Transcript

=head1 SYNOPSIS

use Sanger::CGP::Grass::GenomeData::Transcript;

my $grass_transcript = new Sanger::CGP::Grass::GenomeData::Transcript(-display_id => $accession,
							 -stable_id  => $accession,
							 -length     => $length,
							 -strand     => $Strand,
							 -start      => $start,
							 -end        => $end,
							 -coding_region_start => $coding_region_start,
							 -coding_region_end   => $coding_region_end,
							 -ccds_id   => $ccds,
							 -entrez_id => '',
							 -gene_name => $GeneName,
							 -gene_biotype   => $biotype,
							 -gene_stable_id => $gene_id,
							 -translation_length => $translation_length,
							 -accession          => $ProteinAccession,
							 -exons              => $grass_exons);

=head1 DESCRIPTION

Class to hold transcript details

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::GenomeData::Transcript();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    if ($args{-stable_id})            { $self->stable_id($args{-stable_id}); }
    if ($args{-display_id})           { $self->display_id($args{-display_id}); }
    if ($args{-length})               { $self->length($args{-length}); }
    if ($args{-strand})               { $self->strand($args{-strand}); }
    if ($args{-start})                { $self->start($args{-start}); }
    if ($args{-end})                  { $self->end($args{-end}); }
    if ($args{-coding_region_start}) { $self->coding_region_start($args{-coding_region_start}); }
    if ($args{-coding_region_end})   { $self->coding_region_end($args{-coding_region_end}); }

    if ($args{-ccds_id})        { $self->ccds_id($args{-ccds_id}); }
    if ($args{-entrez_id})      { $self->entrez_id($args{-entrez_id}); }

    if ($args{-gene_name})      { $self->gene_name($args{-gene_name}); }
    if ($args{-gene_stable_id}) { $self->gene_stable_id($args{-gene_stable_id}); }
    if ($args{-gene_biotype})   { $self->gene_biotype($args{-gene_biotype}); }

    if ($args{-translation_length}) { $self->translation_length($args{-translation_length}); }
    if ($args{-accession})          { $self->accession($args{-accession}); }

    if ($args{-exons}) { $self->exons($args{-exons}); }

    return $self;
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 stable_id

  Arg (1)    : $stable_id
  Example    : $stable_id = $Object->stable_id($stable_id);
  Description: stable_id
  Return     : stable_id object

=cut

sub stable_id {
    my $self = shift;
    $self->{stable_id} = shift if @_;
    return $self->{stable_id};
}
#-----------------------------------------------------------------------#

=head2 display_id

  Arg (1)    : $display_id
  Example    : $display_id = $Object->display_id($display_id);
  Description: display_id
  Return     : display_id object

=cut

sub display_id {
    my $self = shift;
    $self->{display_id} = shift if @_;
    return $self->{display_id};
}
#-----------------------------------------------------------------------#

=head2 length

  Arg (1)    : $length
  Example    : $length = $Object->length($length);
  Description: length
  Return     : length object

=cut

sub length {
    my $self = shift;
    $self->{length} = shift if @_;
    return $self->{length};
}
#-----------------------------------------------------------------------#

=head2 strand

  Arg (1)    : $strand
  Example    : $strand = $Object->strand($strand);
  Description: strand
  Return     : strand object

=cut

sub strand {
    my $self = shift;
    $self->{strand} = shift if @_;
    return $self->{strand};
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

=head2 ccds_id

  Arg (1)    : $ccds_id
  Example    : $ccds_id = $Object->ccds_id($ccds_id);
  Description: ccds_id
  Return     : ccds_id object

=cut

sub ccds_id {
    my $self = shift;
    $self->{ccds_id} = shift if @_;
    return $self->{ccds_id};
}
#-----------------------------------------------------------------------#

=head2 entrez_id

  Arg (1)    : $entrez_id
  Example    : $entrez_id = $Object->entrez_id($entrez_id);
  Description: entrez_id
  Return     : entrez_id object

=cut

sub entrez_id {
    my $self = shift;
    $self->{entrez_id} = shift if @_;
    return $self->{entrez_id};
}
#-----------------------------------------------------------------------#

=head2 gene_name

  Arg (1)    : $gene_name
  Example    : $gene_name = $Object->gene_name($gene_name);
  Description: gene_name
  Return     : gene_name object

=cut

sub gene_name {
    my $self = shift;
    $self->{gene_name} = shift if @_;
    return $self->{gene_name};
}
#-----------------------------------------------------------------------#

=head2 gene_stable_id

  Arg (1)    : $gene_stable_id
  Example    : $gene_stable_id = $Object->gene_stable_id($gene_stable_id);
  Description: gene_stable_id
  Return     : gene_stable_id object

=cut

sub gene_stable_id {
    my $self = shift;
    $self->{gene_stable_id} = shift if @_;
    return $self->{gene_stable_id};
}
#-----------------------------------------------------------------------#

=head2 gene_biotype

  Arg (1)    : $gene_biotype
  Example    : $gene_biotype = $Object->gene_biotype($gene_biotype);
  Description: gene_biotype
  Return     : gene_biotype object

=cut

sub gene_biotype {
    my $self = shift;
    $self->{gene_biotype} = shift if @_;
    return $self->{gene_biotype};
}
#-----------------------------------------------------------------------#

=head2 translation_length

  Arg (1)    : $translation_length
  Example    : $translation_length = $Object->translation_length($translation_length);
  Description: translation_length
  Return     : translation_length object

=cut

sub translation_length {
    my $self = shift;
    $self->{translation_length} = shift if @_;
    return $self->{translation_length};
}
#-----------------------------------------------------------------------#

=head2 accession

  Arg (1)    : $accession
  Example    : $accession = $Object->accession($accession);
  Description: accession
  Return     : accession object

=cut

sub accession {
    my $self = shift;
    $self->{accession} = shift if @_;
    return $self->{accession};
}
#-----------------------------------------------------------------------#

=head2 exons

  Arg (1)    : $exons
  Example    : $exons = $Object->exons($exons);
  Description: exon objects
  Return     : reference to an array of exon objects

=cut

sub exons {
    my $self = shift;
    $self->{exons} = shift if @_;
    return $self->{exons};
}
#-----------------------------------------------------------------------#
1;
