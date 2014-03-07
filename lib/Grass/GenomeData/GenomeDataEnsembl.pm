## Grass::GenomeData::GenomeDataEnsembl

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

Class to access genome data from a remote ensembl DB

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Grass::GenomeData::GenomeDataEnsembl;

use Grass::GenomeData::Transcript;
use Grass::GenomeData::Exon;

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    : 
  Example    : $object = new Grass::GenomeData::GenomeDataEnsembl();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    if ($args{-species})      { $self->species($args{-species}); }
    if ($args{-ensembl_api})  { $self->ensembl_api($args{-ensembl_api}); }

    return $self;
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 species

  Arg (0)    : $species
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

  Arg (0)    : $ensembl_api
  Example    : $ensembl_api = $Object->ensembl_api($ensembl_api);
  Description: get/set the ensembl ensembl_api object to use
  Return     : ensembl_api object

=cut

sub ensembl_api {
    my $self = shift;
    my $ensembl_api = shift if @_;
    if ($ensembl_api) {
	$self->set_ensembl($ensembl_api);
    }
    return $self->{ensembl_api};
}
#-----------------------------------------------------------------------#

=head2 registry

  Arg (0)    : $registry
  Example    : $registry = $Object->registry($registry);
  Description: Ensembl registry
  Return     : registry object

=cut

sub registry {
    my $self = shift;
    return $self->{registry};
}
#-----------------------------------------------------------------------#

=head2 slice_ad

  Arg (0)    : $slice_ad
  Example    : $slice_ad = $Object->slice_ad($slice_ad);
  Description: Ensembl slice_adaptor
  Return     : slice_ad object

=cut

sub slice_ad {
    my $self = shift;
    return $self->{slice_ad};
}
#-----------------------------------------------------------------------#

=head2 trans_ad

  Arg (0)    : $trans_ad
  Example    : $trans_ad = $Object->trans_ad($trans_ad);
  Description: Ensembl transcript_adaptor
  Return     : trans_ad object

=cut

sub trans_ad {
    my $self = shift;
    return $self->{trans_ad};
}
#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
sub set_ensembl {
    my ($self, $ensembl_api) = @_;

    unless ($self->{species}) { print "Can not set Ensembl DB. Species not supplied\n"; exit; }

    # put this here because need to pass in which ensembl version to use
    my $sent = 'use lib qw(' . $ensembl_api . '/ensembl/modules/  );
                use Bio::EnsEMBL::Registry;';
    eval $sent; # delays it until run time so that $ensembl is set first

    my $registry = 'Bio::EnsEMBL::Registry';

    if (($self->{species}) eq 'CEREVISIAE') { $self->{species} = 'Saccharomyces cerevisiae'; }
    $registry->clear();
    $registry->load_registry_from_db( -host => 'ensembldb.ensembl.org',
				      -user => 'anonymous');

    unless ($registry->get_adaptor(($self->{species}),'core','slice')) {
	print "could not get connection to remote ensembl registry\n"; 
	exit;
    }

    $self->{ensembl_api} = $ensembl_api;
    $self->{registry} = $registry;
}
#------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
# get a list of genes between 2 coordinates
sub get_gene_list {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    my @names = ();

    # get the slice adaptor if we dont already have one for the species
    unless ($self->{slice_ad}) { 
	$self->{slice_ad} = $self->{registry}->get_adaptor(($self->{species}),'core','slice'); 
	unless ($self->{slice_ad}) { print "could not get slice for $species core\n"; return(""); }
    }

    my $slice = $self->{slice_ad}->fetch_by_region('chromosome', $chr, $start_coord, $end_coord);
    unless ($slice) { $slice = $self->{slice_ad}->fetch_by_region(undef, chr, $start_coord, $end_coord); } # look at every type of structure, not just chromosomes
    
    my $genes = $slice->get_all_Genes();
    foreach my $gene(@$genes) {
	my $name = '';
	my @links = @{$gene->get_all_DBEntries};
	foreach my $link(@links){
	    if ($link->dbname =~ /^HGNC/) { $name = $link->display_id; 
					    last; 
	    } # picks up 'HGNC' and 'HGNC_curated_gene' names
	}
	unless ($name) { $name = $gene->stable_id; }
	push @names, $name;
    }
    my $names = '';
    if (scalar(@names)) { $names = join ',',@names; }

    return($names);
}
#------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
# note that all coordinates are relative to the slice (ie a negative coordinate lies before the slice)
sub fetch_transcripts_by_region {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    # get the slice adaptor
    unless ($self->{slice_ad}) {
	$self->{slice_ad} = $self->{registry}->get_adaptor($self->{species},'core','slice');
	unless ($self->{slice_ad}) { print "could not get slice for " . $self->{species} . " core\n"; exit; }
    }

    if (($self->{species} eq 'DOG') && ($chr eq 'M')) { $chr = 'mt'; }

    # get the slice
    if ($self->{debug}) { print "\nSLICE " . $chr . ':' . $start_coord . '-' . $end_coord . "\n"; }
    my $slice = $self->{slice_ad}->fetch_by_region('chromosome', $chr, $start_coord, $end_coord);
    unless ($slice) { $slice = $self->{slice_ad}->fetch_by_region(undef, $chr, $start_coord, $end_coord); } # look at every type of structure, not just chromosomes
    unless ($slice) { print "slice fetching errors\n"; return(); }

    # get the transcript adaptor
    unless ($self->{trans_ad}) {
	$self->{trans_ad} = $self->{registry}->get_adaptor(($self->{species}),'core','transcript');
    }
	
    # get all the transcripts that are in the slice
    my $transcripts = $self->{trans_ad}->fetch_all_by_Slice($slice);

    $self->{slice_ad}->dbc->disconnect_if_idle();
    $self->{trans_ad}->dbc->disconnect_if_idle();

    # populate a set of Grass Transcript/Exon objects
    my $grass_transcripts = $self->populate_grass_transcripts($transcripts);

    return($grass_transcripts);
}
#--------------------------------------------------------------------------------------------#
sub populate_grass_transcripts {
    my ($self, $ensembl_transcripts) = @_;

    my @grass_transcripts = ();
    unless ($self->{gene_ad}) {
	$self->{gene_ad} = $self->{registry}->get_adaptor($self->{species},'core','gene');
    }

    foreach my $ensembl_transcript(@$ensembl_transcripts) {

	# see if this has ccds/entrez names
	$ccds_id = '';
	$entrez_id = '';
	my @trans_links = @{$ensembl_transcript->get_all_DBEntries};
	foreach my $link(@trans_links) {
	    if    ($link->dbname eq 'CCDS')       { $ccds_id = $link->display_id; }
	    elsif ($link->dbname eq 'EntrezGene') { $entrez_id = $link->primary_id; }
	    last if ($ccds_id && $entrez_id);
	}

	# get gene info
	my $gene = $self->{gene_ad}->fetch_by_transcript_stable_id($ensembl_transcript->stable_id);
	my $gene_name = $gene->stable_id();
	my @gene_links = @{$gene->get_all_DBEntries};
	foreach my $gene_link(@gene_links){
	    if ($gene_link->dbname =~ /^HGNC/) { $gene_name = $gene_link->display_id; last; } # picks up 'HGNC' and 'HGNC_curated_gene' names
	}
	$self->{gene_ad}->dbc->disconnect_if_idle();

	# force translation of the transcript to get the translation length
	my $aa_length = '';
	my $accession = '';
	eval { $ensembl_transcript->translate;  };
	if ($@) { print "translating ensembl transcript failed... $@\n"; }
	elsif ($ensembl_transcript->translate) { 
	    $aa_length = $ensembl_transcript->translate->length; 
	    $accession = $ensembl_transcript->translation->stable_id;
	}
	if ($self->{debug}) { 
            print "name: " . $ensembl_transcript->display_id . "\n";
            print "id: " . $ensembl_transcript->stable_id . "\n";
	    print "Translation_length ensembl: $aa_length\n"; 
	    print "  start ensembl: "  . $ensembl_transcript->start . "\n";
	    print "  end ensembl  : "  . $ensembl_transcript->end . "\n";
	    print "  cds start ensembl: "  . $ensembl_transcript->coding_region_start . "\n";
	    print "  cds end ensembl  : "  . $ensembl_transcript->coding_region_end . "\n";
	}

	# populate the object
        # use gene_name in place of gene->stable_id because otherwise we can't compare to the cache version (which doesn't supply gene->stable_id)
	my $grass_transcript = new Grass::GenomeData::Transcript(-display_id => $ensembl_transcript->display_id,
								 -stable_id  => $ensembl_transcript->stable_id,
								 -length     => $ensembl_transcript->length,
								 -strand     => $ensembl_transcript->strand,
								 -start      => $ensembl_transcript->start,
								 -end        => $ensembl_transcript->end,
								 -coding_region_start => $ensembl_transcript->coding_region_start,
								 -coding_region_end   => $ensembl_transcript->coding_region_end,
								 -ccds_id   => $ccds_id,
								 -entrez_id => $entrez_id,
								 -gene_name      => $gene_name,
								 -gene_biotype   => $gene->biotype,
								 -gene_stable_id => $gene_name,
								 -translation_length => $aa_length,
								 -accession          => $accession);

	# get the ensembl exons for this transcript and populate equivalent grass exon objects
	my @grass_exons = ();
	my @ensembl_exons = @{$ensembl_transcript->get_all_Exons()};
	foreach my $ensembl_exon(@ensembl_exons) {
	    my $grass_exon = $self->populate_grass_exon($ensembl_exon, $ensembl_transcript);
	    push @grass_exons, $grass_exon;
	}
	# add to the transcript object
	$grass_transcript->exons(\@grass_exons);

	# push to the array of new transcripts
	push @grass_transcripts, $grass_transcript;
    }
    return(\@grass_transcripts);
}
#--------------------------------------------------------------------------------------------#
sub populate_grass_exon {
    my ($self, $ensembl_exon, $ensembl_transcript) = @_;

    my $seq = '';
    if ($ensembl_exon->seq) { $seq = $ensembl_exon->seq->seq; }
    # populate the object
    my $grass_exon = new Grass::GenomeData::Exon(-phase  => $ensembl_exon->phase,
						 -start  => $ensembl_exon->start,
						 -end    => $ensembl_exon->end,
						 -coding_region_start => $ensembl_exon->coding_region_start($ensembl_transcript),
						 -coding_region_end   => $ensembl_exon->coding_region_end($ensembl_transcript),
						 -seq    => $seq );
    return($grass_exon);
}
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

1;
