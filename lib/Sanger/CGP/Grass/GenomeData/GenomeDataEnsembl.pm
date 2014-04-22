## Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl

#
# Author las
#
=head1 NAME

GenomeDataEnsembl

=head1 SYNOPSIS

use Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl;

$GD_object = new Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl(-species     => $species,
						      -ensembl_api => $ensembl_api,
						      -gene_id_required => 0,
						      -use_all_biotypes => 1);

=head1 DESCRIPTION

Class to access genome data from a remote ensembl DB server

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl;

use strict;

use Sanger::CGP::Grass::GenomeData::Transcript;
use Sanger::CGP::Grass::GenomeData::Exon;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    : 
  Example    : $object = new Sanger::CGP::Grass::GenomeData::GenomeDataEnsembl();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    if (defined($args{-gene_id_required})) { $self->gene_id_required($args{-gene_id_required}); }
    if (defined($args{-use_all_biotypes})) { $self->use_all_biotypes($args{-use_all_biotypes}); } # set this before setting genome_cache file
    if ($args{-species})      { $self->species($args{-species}); }
    if ($args{-ensembl_api})  { $self->ensembl_api($args{-ensembl_api}); }

    return $self;
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 gene_id_required

  Arg (1)    : 1/0
  Example    : $gene_id_required = $Object->gene_id_required($gene_id_required);
  Description: define whether gene_id id is required
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
  Description: define whether all transcript biotypes are required
  Return     : 1/0

=cut

sub use_all_biotypes {
    my $self = shift;
    $self->{use_all_biotypes} = shift if @_;
    return $self->{use_all_biotypes};
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
	$self->_set_ensembl($ensembl_api);
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

  Arg (1)    : $slice_ad
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

  Arg (1)    : $trans_ad
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
sub _set_ensembl {
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

=head2 get_gene_list

  Arg (3)    : $chr, $start_coord, $end_coord
  Example    : $gene_names = $Object->get_gene_list($chr, $start_coord, $end_coord);
  Description: get a list of gene_names that lie between a pair of coordinates
  Return     : $gene_names

=cut

# get a list of genes between 2 coordinates
sub get_gene_list {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    my @names = ();

    # get the slice adaptor if we dont already have one for the species
    unless ($self->{slice_ad}) { 
	$self->{slice_ad} = $self->{registry}->get_adaptor(($self->{species}),'core','slice'); 
	unless ($self->{slice_ad}) { print "could not get slice for $self->{species} core\n"; return(""); }
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

=head2 fetch_transcripts_by_region

  Arg (3)    : $chr, $start_coord, $end_coord
  Example    : $gene_names = $Object->fetch_transcripts_by_region($chr, $start_coord, $end_coord);
  Description: get transcripts that hit a specific genomic region
  Return     : reference to an array of $grass_transcript_objects

=cut

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
    my $grass_transcripts = $self->_populate_grass_transcripts($transcripts);

    return($grass_transcripts);
}
#--------------------------------------------------------------------------------------------#
sub _populate_grass_transcripts {
    my ($self, $ensembl_transcripts) = @_;

    my @grass_transcripts = ();
    unless ($self->{gene_ad}) {
	$self->{gene_ad} = $self->{registry}->get_adaptor($self->{species},'core','gene');
    }

    foreach my $ensembl_transcript(@$ensembl_transcripts) {

	# if required, skip if the transcript biotype is only putative coding
	unless ($self->{use_all_biotypes}) {
	    if ($self->{debug}) { print "TRANS BIOTYPE: " . $ensembl_transcript->biotype . "\n"; }
	    next unless (grep { $_ eq $ensembl_transcript->biotype } ('protein_coding', 'lincRNA', 'rRNA', 'snRNA', 'snoRNA', 'miRNA'));
	    next unless (($ensembl_transcript->status()) eq 'KNOWN');
	}

	if ($self->{debug}) { print " DOIT\n"; }

	# see if this has ccds/entrez names
	my $ccds_id = '';
	my $entrez_id = '';
	my @trans_links = @{$ensembl_transcript->get_all_DBEntries};
	foreach my $link(@trans_links) {
	    if    ($link->dbname eq 'CCDS')       { $ccds_id = $link->display_id; }
	    elsif ($link->dbname eq 'EntrezGene') { $entrez_id = $link->primary_id; }
	    last if ($ccds_id && $entrez_id);
	}

	# get gene info
	my ($gene_name, $gene_id, $gene_biotype, $gene_status) = $self->_getGeneInfoForTranscript($ensembl_transcript); 

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
	my $grass_transcript = new Sanger::CGP::Grass::GenomeData::Transcript(-display_id => $ensembl_transcript->display_id,
								 -stable_id  => $ensembl_transcript->stable_id,
								 -length     => $ensembl_transcript->length,
								 -strand     => ($ensembl_transcript->strand . ''),
								 -start      => $ensembl_transcript->start,
								 -end        => $ensembl_transcript->end,
								 -coding_region_start => $ensembl_transcript->coding_region_start,
								 -coding_region_end   => $ensembl_transcript->coding_region_end,
								 -ccds_id   => $ccds_id,
								 -entrez_id => $entrez_id,
								 -gene_name      => $gene_name,
								 -gene_biotype   => $gene_biotype,
								 -gene_stable_id => $gene_id,
								 -translation_length => $aa_length,
								 -accession          => $accession);

	# get the ensembl exons for this transcript and populate equivalent grass exon objects
	my @grass_exons = ();
	my @ensembl_exons = @{$ensembl_transcript->get_all_Exons()};
	foreach my $ensembl_exon(@ensembl_exons) {
	    my $grass_exon = $self->_populate_grass_exon($ensembl_exon, $ensembl_transcript);
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

sub _getGeneInfoForTranscript {
	my ($self,$transcript) = @_;

	my $ga = $self->{registry}->get_adaptor($self->{species},'core','gene');
	my $gene = $ga->fetch_by_transcript_stable_id($transcript->stable_id);
	my $biotype = $gene->biotype();
	my @links = @{$gene->get_all_DBEntries};
	my $gene_name_hgnc = '';
	my $gene_name_vega = '';
	my $gene_name_uniprot = '';
	my $gene_name_ensembl = '';
	my $gene_name_rfam = '';
	foreach my $gene_link(@links){
	    if ($gene_link->dbname =~ /^HGNC/) { $gene_name_hgnc = $gene_link->display_id; } # picks up 'HGNC' and 'HGNC_curated_gene' names
	    if ($gene_link->dbname =~ /Clone_based_vega_gene/) { $gene_name_vega = $gene_link->display_id; } # picks up 'Vega' and 'Vega_curated_gene' names
	    if ($gene_link->dbname =~ /Clone_based_ensembl/) { $gene_name_ensembl = $gene_link->display_id; } # picks up 'Ensembl' names
	    if ($gene_link->dbname =~ /Uniprot_gn/) { $gene_name_uniprot = $gene_link->display_id; } # picks up 'Uniprot' names
	    if ($gene_link->dbname =~ /RFAM/) { $gene_name_rfam = $gene_link->display_id; } # picks up 'Uniprot' names
	}

	my $gene_name = $gene->stable_id();
	if ($gene_name_hgnc)       { $gene_name = $gene_name_hgnc; }
	elsif ($gene_name_vega)    { $gene_name = $gene_name_vega; }
	elsif ($gene_name_uniprot) { $gene_name = $gene_name_uniprot; }
	elsif ($gene_name_ensembl) { $gene_name = $gene_name_ensembl; }
	elsif ($gene_name_rfam)    { $gene_name = $gene_name_rfam; }
	my $gene_id = $gene->stable_id();
	my $gene_status = $gene->status();
	$ga->dbc->disconnect_if_idle();

	if ($self->{gene_id_required}) { return ($gene_name, $gene_id, $biotype, $gene_status); }
	else                           { return ($gene_name, $gene_name, $biotype, $gene_status); }
}
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
sub _populate_grass_exon {
    my ($self, $ensembl_exon, $ensembl_transcript) = @_;

    my $seq = '';
    if ($ensembl_exon->seq) { $seq = $ensembl_exon->seq->seq; }
    # populate the object
    my $grass_exon = new Sanger::CGP::Grass::GenomeData::Exon(-phase  => $ensembl_exon->phase,
						 -end_phase => $ensembl_exon->end_phase,
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
