package Sanger::CGP::Grass::GenomeData;

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

use Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl;
use Sanger::CGP::Grass::GenomeData::GenomeDataCache;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;


=head1 NAME

GenomeData

=head1 SYNOPSIS

use Sanger::CGP::Grass::GenomeData;

my $genome_data_cache = new Sanger::CGP::Grass::GenomeData(-genome_cache => $genome_cache,
					      -gene_id_required => 0);

or

my $genome_data_ensembl = new Sanger::CGP::Grass::GenomeData(-species     => $species,
						-ensembl_api => $ensembl_api_74,
						-gene_id_required => 0,
						-use_all_biotypes => 0);

depending on whether a genome_cache file (generated using vagrent from Ensembl files) or access to the ensembl server is available.

=head1 DESCRIPTION

Class to access genome data from either ensembl DB or a cached flat file version

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::GenomeData();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    if (defined($args{-gene_id_required})) { $self->gene_id_required($args{-gene_id_required}); } # set this before setting genome_cache file
    if (defined($args{-use_all_biotypes})) { $self->use_all_biotypes($args{-use_all_biotypes}); } # set this before setting genome_cache file
    if ($args{-genome_cache}) { $self->genome_cache($args{-genome_cache}); }
    if ($args{-species})      { $self->species($args{-species}); }
    if ($args{-ensembl_api})  { $self->ensembl_api($args{-ensembl_api}); }

    return $self;
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 gene_id_required

  Arg (1)    : 1/0
  Example    : $gene_id_required = $Object->gene_id_required($gene_id_required);
  Description: define whether gene_id (id, not just name) is required
  Return     : 1/0

=cut

sub gene_id_required {
    my $self = shift;
    $self->{gene_id_required} = shift if @_;
    return $self->{gene_id_required};
}
#-----------------------------------------------------------------------#

=head2 use_all_biotypes

  Arg (1)    : 1/0
  Example    : $use_all_biotypes = $Object->use_all_biotypes($use_all_biotypes);
  Description: define whether all transcript biotypes are required (1) or just protein_coding (0)
  Return     : 1/0

=cut

sub use_all_biotypes {
    my $self = shift;
    $self->{use_all_biotypes} = shift if @_;
    return $self->{use_all_biotypes};
}
#-----------------------------------------------------------------------#

=head2 genome_cache

  Arg (1)    : $genome_cache
  Example    : $genome_cache = $Object->genome_cache($genome_cache);
  Description: get/set the ensembl genome_cache object to use
  Return     : genome_cache object

=cut

sub genome_cache {
    my $self = shift;
    my $genome_cache = shift if @_;
    if ($genome_cache) {
	$self->{cache} = new Sanger::CGP::Grass::GenomeData::GenomeDataCache(-genome_cache => $genome_cache,
								-gene_id_required => $self->{gene_id_required});
    }
    return $self->{genome_cache};
}
#-----------------------------------------------------------------------#

=head2 species

  Arg (1)    : $species
  Example    : $species = $Object->species($species);
  Description: which species to get from Ensembl
  Return     : species object

=cut

sub species {
    my $self = shift;
    $self->{species} = shift if @_;
    return $self->{species};
}
#-----------------------------------------------------------------------#

=head2 ensembl_api

  Arg (1)    : $ensembl_api
  Example    : $ensembl_api = $Object->ensembl_api($ensembl_api);
  Description: get/set the ensembl ensembl_api object to use
  Return     : ensembl_api object

=cut

sub ensembl_api {
    my $self = shift;
    my $ensembl_api = shift if @_;
    if ($ensembl_api) {
	$self->{ensembl} = new Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl(-species     => $self->{species},
								    -ensembl_api => $ensembl_api,
								    -gene_id_required => $self->{gene_id_required},
								    -use_all_biotypes => $self->{use_all_biotypes});
    }
    return $self->{ensembl}->{ensembl_api};
}
#-----------------------------------------------------------------------#

=head2 registry

  Arg (1)    : $registry
  Example    : $registry = $Object->registry($registry);
  Description: Ensembl registry
  Return     : registry object

=cut

sub registry {
    my $self = shift;
    return $self->{ensembl}->{registry};
}
#-----------------------------------------------------------------------#

=head2 slice_ad

  Arg (1)    : $slice_ad
  Example    : $slice_ad = $Object->slice_ad($slice_ad);
  Description: Ensembl slice_adaptor
  Return     : slice_ad object

=cut

sub slice_ad {
    my $self = shift;
    return $self->{ensembl}->{slice_ad};
}
#-----------------------------------------------------------------------#

=head2 trans_ad

  Arg (1)    : $trans_ad
  Example    : $trans_ad = $Object->trans_ad($trans_ad);
  Description: Ensembl transcript_adaptor
  Return     : trans_ad object

=cut

sub trans_ad {
    my $self = shift;
    return $self->{ensembl}->{trans_ad};
}
#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#

=head2 get_gene_list

  Arg (3)    : $chr, $start_coord, $end_coord
  Example    : $gene_names = $Object->get_gene_list($chr, $start_coord, $end_coord);
  Description: get a list of gene_names that lie between a pair of coordinates
  Return     : $gene_names

=cut

sub get_gene_list {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    my $names = '';
    if ($self->{cache}) {
	$names = $self->{cache}->get_gene_list($chr, $start_coord, $end_coord);
    }
    elsif ($self->{ensembl}) {
	$names = $self->{ensembl}->get_gene_list($chr, $start_coord, $end_coord);
    }

    return($names);
}
#---------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

=head2 fetch_transcripts_by_region

  Arg (3)    : $chr, $start_coord, $end_coord
  Example    : $gene_names = $Object->fetch_transcripts_by_region($chr, $start_coord, $end_coord);
  Description: get transcripts that hit a specific genomic region
  Return     : reference to an array of $grass_transcript_objects

=cut

# note that all coordinates are relative to the slice (ie a negative coordinate lies before the slice)
sub fetch_transcripts_by_region {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    my $transcripts = undef;

    if ($self->{cache}) {
	$transcripts = $self->{cache}->fetch_transcripts_by_region($chr, $start_coord, $end_coord);
    }
    elsif ($self->{ensembl}) {
	$transcripts = $self->{ensembl}->fetch_transcripts_by_region($chr, $start_coord, $end_coord);
    }
    return($transcripts);
}
#--------------------------------------------------------------------------------------------#

=head2 thin_out_translist

  Arg (1)    : reference to an array of grass transcript objects
  Example    : ($transcripts, ccds_only) = $Object->thin_out_translist($transcripts);
  Description: pick the best transcripts from a list. also return a ccds_only flag - whether the list has been filtered to only contain ccds entries.
  Return     : reference to an array of $grass_transcript_objects, $ccds_only (1 or 0)

=cut

# this is filtering grass transcript objects (not ensembl transcript objects)
sub thin_out_translist {
    my ($self, $translist) = @_;

    my @new_translist = ();
    my $ccds_only = 0;
    # check to see if any have a consensus coding sequence id (only available for human and mouse)
    # this will reduce the number of rg_gene_anno entries by about 1/3rd
    my $ccds_total = 0;
    foreach my $trans(@$translist){
	if ($trans->ccds_id()) {
	    $ccds_total++;
	    if ($self->{debug}) { print $trans->display_id . " (CCDS)\n"; }
	    push @new_translist, $trans;
	}
    }

    if ($self->{debug}) { print "\n ccds count $ccds_total, trans count " . (scalar(@$translist)) . "\n"; }
    if (scalar(@new_translist)) {
	$ccds_only = 1; # mark this end as using ccds entries only
	return(\@new_translist, $ccds_only);
    }

    # if it can't be restricted, use what we've got
    return($translist, $ccds_only);
}
#--------------------------------------------------------------------------------------------#

1;
