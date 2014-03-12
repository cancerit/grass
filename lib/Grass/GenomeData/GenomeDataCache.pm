## Grass::GenomeData::GenomeDataCache

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

Class to access genome data from a cached flat file. 

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Grass::GenomeData::GenomeDataCache;

use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;
use Bio::DB::Sam;
use Data::Dumper;

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    : 
  Example    : $object = new Grass::GenomeData::GenomeDataCache();
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
    if ($args{-genome_cache}) { $self->genome_cache($args{-genome_cache}); }

    return $self;
}

#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 gene_id_required

  Arg (0)    : 1/0
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

=head2 genome_cache

  Arg (0)    : $genome_cache
  Example    : $genome_cache = $Object->genome_cache($genome_cache);
  Description: get/set the ensembl genome_cache object to use
  Return     : genome_cache object

=cut

sub genome_cache {
    my $self = shift;
    my $genome_cache = shift if @_;
    if ($genome_cache) {
	$self->check_genome_cache($genome_cache);
    }
    return $self->{genome_cache};
}
#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#
sub check_genome_cache {
    my ($self, $cache_file) = @_;

    my $root_file_name = $cache_file;
    $root_file_name =~ s/.cache.gz\s*//;

    my $cache_tbi_file = $root_file_name . '.cache.gz.tbi';
    my $transcript_file = $root_file_name . '.transcript.fa';
    my $transcript_fai_file = $root_file_name . '.transcript.fa.fai';

    unless (-e "$cache_file") { print "genome cache file $cache_file not found\n"; exit; }
    unless (-e "$cache_tbi_file") { print "genome cache tbi file $cache_tbi_file not found\n"; exit; }
    unless (-e "$transcript_file") { print "genome transcript file $transcript_file not found\n"; exit; }
    unless (-e "$transcript_fai_file") { print "genome transcript fai file $transcript_fai_file not found\n"; exit; }

    $self->{genome_cache} = $cache_file;
}
#------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
sub get_gene_list {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    my $names = '';

    my $cache = $self->{genome_cache};

    # eg tabix Homo_sapiens.GRCh37.74.vagrent.cache.gz 1:10000000-10100000 | cut -f5 | uniq
    my $coord_string = $chr . ':' . $start_coord . '-' . $end_coord;
    my $results = `tabix $cache $coord_string | cut -f5 | uniq`;
    chomp $results;

    if ($results) {
	my @results = split ',', $results;
	$names = join ',', @results;
    }

    return($names);
}
#------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
# note that all coordinates are relative to the slice (ie a negative coordinate lies before the slice)
sub fetch_transcripts_by_region {
    my ($self, $chr, $start_coord, $end_coord) = @_;

    my @vagrent_transcripts = ();
    my @transcripts = ();

    my $cache = $self->{genome_cache};
    if (($cache =~ /Canis/) && ($chr eq 'M')) { $chr = 'mt'; } # dogs have motochondrial chr naming problem

    my $coord_string = $chr . ':' . $start_coord . '-' . $end_coord;
    my $results = `tabix $cache $coord_string | cut -f6,7`;
    my @results = split "\n", $results;
    foreach my $res(@results) {
	my @res = split "\t", $res;
	my $length = $res[0];
	my $trans = $res[1];
	my $VAR1;
	eval $trans;
	my $vagrent_transcript = $VAR1;
	push @vagrent_transcripts, [$length, $vagrent_transcript];
    }

    # populate a set of Grass Transcript/Exon objects
    my $grass_transcripts = $self->populate_grass_transcripts(\@vagrent_transcripts, $chr, $start_coord);

    return($grass_transcripts);
}
#--------------------------------------------------------------------------------------------#
sub populate_grass_transcripts {
    my ($self, $vagrent_transcripts, $chr, $start_coord) = @_;

    my @grass_transcripts = ();

    foreach my $entry(@$vagrent_transcripts) {
	my $length = $entry->[0];
	my $vagrent_transcript = $entry->[1];

	my $biotype = $vagrent_transcript->getGeneType;
	if ($biotype eq 'ProteinCoding') { $biotype = 'protein_coding'; }

	my $coding_region_length = $vagrent_transcript->getCdsLength;
	my $translation_length = 0;
	if ($coding_region_length) { $translation_length = int($coding_region_length/3) - 1; } # cds length divided by three and remove the stop codon  

	my $start = $vagrent_transcript->getGenomicMinPos - $start_coord + 1;
	my $end = $vagrent_transcript->getGenomicMaxPos - $start_coord + 1;

	# cds min and max are coordinates on the cDNA so have to work out what's going on from the exons
	# get the array of exon objects at the same time
	my ($grass_exons, $coding_region_start, $coding_region_end) = $self->process_exons($vagrent_transcript, $chr, $start_coord);

	if ($self->{debug}) { 
	    print "name: " . $vagrent_transcript->getAccession . "\n";
	    print "Translation_length cache: |$coding_region_length| |$translation_length|\n"; 
	    print "  start cache: $start\n";
	    print "  end cache  : $end\n";
	    print "  cds start cache: $coding_region_start\n";
	    print "  cds end cache  : $coding_region_end\n";
	}

	my $gene_id = $vagrent_transcript->getGeneName;
	if ($self->{gene_id_required} && ($vagrent_transcript->getGeneId)) { $gene_id = $vagrent_transcript->getGeneId; }

	# populate the object (convert coordinates so they are relative to the coordinate slice requested) 
	# Need to force strand to string form so its the same as for ensembl and can compare the two
	my $grass_transcript = new Grass::GenomeData::Transcript(-display_id => $vagrent_transcript->getAccession,
								 -stable_id  => $vagrent_transcript->getAccession,
								 -length     => $length,
								 -strand     => ($vagrent_transcript->getStrand . ''),
								 -start      => $start,
								 -end        => $end,
								 -coding_region_start => $coding_region_start,
								 -coding_region_end   => $coding_region_end,
								 -ccds_id   => $vagrent_transcript->getCCDS,
								 -entrez_id => '',
								 -gene_name      => $vagrent_transcript->getGeneName,
								 -gene_biotype   => $biotype,
								 -gene_stable_id => $gene_id,
								 -translation_length => $translation_length,
								 -accession          => $vagrent_transcript->getProteinAccession,
								 -exons => $grass_exons);
	# push to the array of new transcripts
	push @grass_transcripts, $grass_transcript;
    }
    return(\@grass_transcripts);
}
#--------------------------------------------------------------------------------------------#
sub process_exons {
    my ($self, $vagrent_transcript, $chr, $start_coord) = @_;

    # need this, relative to slice
    my $transcript_cds_start = undef;
    my $transcript_cds_end = undef;

    my $transcript_cds_start_on_cDNA = $vagrent_transcript->getCdsMinPos;
    my $transcript_cds_end_on_cDNA =  $vagrent_transcript->getCdsMaxPos;
    my $transcript_cds_phase =  $vagrent_transcript->getCdsPhase; # 0 usually, except when transcript/cds is truncated at the start
    my $strand =  $vagrent_transcript->getStrand;

    # retrieve the transcript sequence from the other flat file (this contains the sequence of the transcript on the transcript's strand)
    my $transcript_seq = $self->get_seq($vagrent_transcript->getAccession);

    # get the vagrent exons for this transcript and populate equivalent grass exon objects
    # also, convert the transcript cds start/end coord to coordinates relative to the slice
    my $exon_phase = -1;
    my $exon_end_phase = -1;
    my $check_coding = 0;
    my @grass_exons = ();
    my @vagrent_exons = ();
    my $vagrent_exons = $vagrent_transcript->{_exons};

    foreach my $vagrent_exon(@$vagrent_exons) {

	# to keep track of where we are in the coding sequence
	my $contains_cds_start = 0;
	my $contains_cds_end = 0;

        # need this, relative to slice
	my $exon_cds_start = undef;
	my $exon_cds_end = undef;

	my $exon_cDNA_start = $vagrent_exon->getRnaMinPos();
	my $exon_cDNA_end = $vagrent_exon->getRnaMaxPos();
	my $exon_genomic_start = $vagrent_exon->getMinPos();
	my $exon_genomic_end = $vagrent_exon->getMaxPos();
	my $exon_start = $exon_genomic_start - $start_coord + 1; # relative to slice
	my $exon_end   = $exon_genomic_end - $start_coord + 1; # relative to slice

	# retrieve the exon sequence from transcript sequence
	my $exon_seq = substr($transcript_seq, ($exon_cDNA_start - 1), ($exon_cDNA_end - $exon_cDNA_start + 1));

	# check where we are in the coding sequence...
	if ($transcript_cds_start_on_cDNA && ($transcript_cds_start_on_cDNA >= $exon_cDNA_start) && ($transcript_cds_start_on_cDNA <= $exon_cDNA_end)) { $contains_cds_start = 1; }
	if ($transcript_cds_end_on_cDNA   && ($transcript_cds_end_on_cDNA   >= $exon_cDNA_start) && ($transcript_cds_end_on_cDNA   <= $exon_cDNA_end)) { $contains_cds_end = 1; }

	if ($contains_cds_start) { $check_coding = 1; } # this is the start coding exon
	if ($contains_cds_end)   { $check_coding = 3; } # this is the end coding exon
	if ($check_coding && 
            !$contains_cds_start && 
            !$contains_cds_end)  { $check_coding = 2; } # this is a mid coding exon

	# treat exons that are start exons BUT the cds covers the full length as mids
	if ($contains_cds_start && 
	    $contains_cds_end && 
	    ($exon_cDNA_start == $transcript_cds_start_on_cDNA) && 
	    ($exon_cDNA_end == $transcript_cds_end_on_cDNA)) { $check_coding = 2; }
	if ($contains_cds_start && !($contains_cds_end)   && ($exon_cDNA_start == $transcript_cds_start_on_cDNA)) { $check_coding = 2; }
	if ($contains_cds_end   && !($contains_cds_start) && ($exon_cDNA_end == $transcript_cds_end_on_cDNA))     { $check_coding = 2; }

	# get the transcript cds start/end coords relative to the slice
	if ($contains_cds_start) {
	    if    ($strand eq '1')  { $transcript_cds_start = $exon_start + ($transcript_cds_start_on_cDNA - $exon_cDNA_start); }
	    elsif ($strand eq '-1') { $transcript_cds_end   = $exon_end   - ($transcript_cds_start_on_cDNA - $exon_cDNA_start); }
	}
	if ($contains_cds_end) {
	    if    ($strand eq '1')  { $transcript_cds_end   = $exon_end   - ($exon_cDNA_end - $transcript_cds_end_on_cDNA); }
	    elsif ($strand eq '-1') { $transcript_cds_start = $exon_start + ($exon_cDNA_end - $transcript_cds_end_on_cDNA); }
	}

	# set the exon coding region start and end
	if ($check_coding) {
	    if ($contains_cds_start) {
		if    ($strand eq '1')  { $exon_cds_start = $transcript_cds_start; }
		elsif ($strand eq '-1') { $exon_cds_end   = $transcript_cds_end; }
	    }
	    if ($contains_cds_end) {
		if    ($strand eq '1')  { $exon_cds_end   = $transcript_cds_end; }
		elsif ($strand eq '-1') { $exon_cds_start = $transcript_cds_start; }
	    }
	    unless (defined($exon_cds_start)) { $exon_cds_start = $exon_start; }
	    unless (defined($exon_cds_end))   { $exon_cds_end   = $exon_end; }
	}

	# the phase of the current exon is the same as the end phase of the previous exon
	$exon_phase = $exon_end_phase;
	# BUT if this is a 5' truncated transcript and the exon start and transcript start are the same, use the transcript cds phase as the start phase
	if ($contains_cds_start && ($exon_cDNA_start == $transcript_cds_start_on_cDNA)) { $exon_phase = $transcript_cds_phase; }

	# work out the end_phase for this exon
	# if this is a coding exon, calc how many coding bases there are to the end of the exon, divide by 3, and get the remainder
	#    this number is used for calculating phase of the NEXT exon, not this one.
	if ($check_coding == 1) {
	    $exon_end_phase = ($exon_cDNA_end - $transcript_cds_start_on_cDNA + 1) % 3;
	    if ($self->{debug}) { print "first_exon: phase $exon_phase end_phase $exon_end_phase\n"; }
	}
	elsif ($check_coding == 2) {
	    $exon_end_phase = ($exon_cDNA_end - ($exon_cDNA_start - $exon_phase) + 1) % 3;
	    if ($self->{debug}) { print "mid_exon: phase $exon_phase end_phase $exon_end_phase\n"; }
	}
	elsif ($check_coding == 3) {
	    $exon_end_phase = -1;
	    if ($self->{debug}) { print "end_exon: phase $exon_phase end_phase $exon_end_phase\n"; }
	}

	if ($self->{debug}) { print "exon_coord |gstart $exon_genomic_start| |gend $exon_genomic_end| |start $exon_start| |end $exon_end| |cdcsstart $exon_cds_start| |cdsend $exon_cds_end| phases |$exon_phase| |$exon_end_phase|\n"; }

	# initialise grass exon object
	my $grass_exon = new Grass::GenomeData::Exon(-phase  => $exon_phase,
						     -end_phase  => $exon_end_phase,
						     -start  => $exon_start,
						     -end    => $exon_end,
						     -coding_region_start => $exon_cds_start,
						     -coding_region_end   => $exon_cds_end,
						     -seq    => $exon_seq );
	push @grass_exons, $grass_exon;
    }

    return(\@grass_exons, $transcript_cds_start, $transcript_cds_end);
}
#--------------------------------------------------------------------------------------------#
sub get_seq {
    my ($self, $transcript) = @_;

    unless ($self->{fai}) { 
	my $fasta_ref = $self->{genome_cache};
	$fasta_ref =~ s/cache.gz/transcript.fa/;
	$self->{fai} = Bio::DB::Sam::Fai->load( $fasta_ref ); 
    }
    my $seq = $self->{fai}->fetch($transcript);

    return($seq);
}
#--------------------------------------------------------------------------------------------#

1;
