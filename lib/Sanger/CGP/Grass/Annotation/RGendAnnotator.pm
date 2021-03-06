package Sanger::CGP::Grass::Annotation::RGendAnnotator;

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
use Sanger::CGP::Grass::Anno;
use Sanger::CGP::Grass::AnnoPoint;
use Carp qw(croak);
use Try::Tiny qw(try catch);

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

use constant REGION_TYPE_CDS => 'CDS';
use constant REGION_TYPE_MAYBE_CDS => 'maybeCDS'; # overlapping CDS
use constant REGION_TYPE_EXON => 'exon';
use constant REGION_TYPE_INTRON => 'intron';
use constant REGION_TYPE_NONCDS => 'nonCDS';
use constant REGION_TYPE_5NONCDS => '5nonCDS'; # overlapping 5UTR
use constant REGION_TYPE_3NONCDS => '3nonCDS'; # overlapping 3UTR
use constant REGION_TYPE_UTRINTRON => 'UTRintron';
use constant REGION_TYPE_UTREXON => 'UTRexon';
use constant REGION_TYPE_UTR => 'UTR';
use constant REGION_TYPE_5UTRINTRON => '5UTRintron';
use constant REGION_TYPE_5UTREXON => '5UTRexon';
use constant REGION_TYPE_5UTR => '5UTR';
use constant REGION_TYPE_3UTRINTRON => '3UTRintron';
use constant REGION_TYPE_3UTREXON => '3UTRexon';
use constant REGION_TYPE_3UTR => '3UTR';
use constant REGION_TYPE_UPSTREAM => 'upstream';
use constant REGION_TYPE_DOWNSTREAM => 'downstream';

=head1 NAME

RGendAnnotator

=head1 SYNOPSIS

use Sanger::CGP::Grass::Annotation::RGendAnnotator;

my $End = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry   => $entry,
						-end     => $end,
						-within  => $within,
						-genome_data => $genome_data_access_object);

my $anns = $End->annotate();

=head1 DESCRIPTION

Class to do the annotation of one end of a rearrangement

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::Annotation::RGendAnnotator();
  Description: make a new RGannotator object
  Return     : RGannotator object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    # set the input parameters

    # get the object by id
    if ($args{-entry}) { $self->entry($args{-entry}); }
    if ($args{-end}) { $self->end($args{-end}); }
    if ($args{-genome_data}) { $self->genome_data($args{-genome_data}); }
    if ($args{-within}) { $self->within($args{-within}); }
    if (defined($args{-entrez_required})) { $self->entrez_required($args{-entrez_required}); }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 entry

  Arg (1)    : $entry
  Example    : $entry = $Object->entry($entry);
  Description: get/set the data entry object to use - contains info on the rearrangement
  Return     : data entry object

=cut

sub entry {
    my $self = shift;
    $self->{entry} = shift if @_;
    return $self->{entry};
}
#-----------------------------------------------------------------------#

=head2 end

  Arg (1)    : $end
  Example    : $end = $Object->end($end);
  Description: get/set the end (1 or 2)
  Return     : end object

=cut

sub end {
    my $self = shift;
    $self->{end} = shift if @_;
    if ($self->{end} eq 'L') { $self->{end} = 1; }
    if ($self->{end} eq 'H') { $self->{end} = 2; }
    return $self->{end};
}
#-----------------------------------------------------------------------#

=head2 genome_data

  Arg (1)    : $genome_data
  Example    : $genome_data = $Object->genome_data($genome_data);
  Description: get/set the genome_data object to use
  Return     : genome_data object

=cut

sub genome_data {
    my $self = shift;
    $self->{genome_data} = shift if @_;
    return $self->{genome_data};
}
#-----------------------------------------------------------------------#

=head2 within

  Arg (1)    : $within
  Example    : $within = $Object->within($within);
  Description: distance either side, within which to look for genes/transcripts
  Return     : within object

=cut

sub within {
    my $self = shift;
    $self->{within} = shift if @_;
    return $self->{within};
}
#-----------------------------------------------------------------------#

=head2 entrez_required

  Arg (1)    : $entrez_required
  Example    : $entrez_required = $Object->entrez_required($entrez_required);
  Description: get the entrez id
  Return     : entrez_required object

=cut

sub entrez_required {
    my $self = shift;
    $self->{entrez_required} = shift if @_;
    return $self->{entrez_required};
}
#-----------------------------------------------------------------------#

=head2 chr

  Arg (1)    : $chr
  Example    : $chr = $Object->chr($chr);
  Description: which chr to get from Ensembl
  Return     : chr object

=cut

sub chr {
    my $self = shift;
    $self->{chr} = shift if @_;
    return $self->{chr};
}
#-----------------------------------------------------------------------#

=head2 pos_start

  Arg (1)    : $pos_start
  Example    : $pos_start = $Object->pos_start($pos_start);
  Description: which pos_start to get from Ensembl
  Return     : pos_start object

=cut

sub pos_start {
    my $self = shift;
    $self->{pos_start} = shift if @_;
    return $self->{pos_start};
}
#-----------------------------------------------------------------------#

=head2 pos_end

  Arg (1)    : $pos_end
  Example    : $pos_end = $Object->pos_end($pos_end);
  Description: which pos_end to get from Ensembl
  Return     : pos_end object

=cut

sub pos_end {
    my $self = shift;
    $self->{pos_end} = shift if @_;
    return $self->{pos_end};
}
#-----------------------------------------------------------------------#

=head2 rg_annos

  Arg (0)    :
  Example    : $rg_annos = $Object->rg_annos($rg_annos);
  Description: array of all partial rg_anno (fusion) objects for this rearrangement end (filled in only for this end)
  Return     : array of partially populated rg_anno objects

=cut

sub rg_annos {
    my $self = shift;
    return $self->{rg_annos};
}
#-----------------------------------------------------------------------#

=head2 ccds_only

  Arg (0)    :
  Example    : $ccds_only = $Object->ccds_only();
  Description: getter - whether this end returns ccds transcripts only (restricts the number of transcripts)
  Return     : 1/0

=cut

sub ccds_only {
    my $self = shift;
    return $self->{ccds_only};
}
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

=head2 annotate

  Arg (0)    :
  Example    : $Object->annotate();
  Description: Gets all the transcripts over this rearrangement end and determines which region the rearrangement disrupts,
               produces an array of pairs of RGannoPoint objects (5 and 3 prime ends for each)
  Return     : Reference to an array of anno objects (each containing a pair of annoPoint objects)

=cut


sub annotate {
    my $self = shift;

    $self->{rg_annos} = undef;

    my $transList = $self->_getTranscripts();

    if(defined($transList) && defined($transList->[0])){

	# if there are more than 4 transcripts, see if any are ccds entries and use those
	# changed to 1 since want ccds if there is one
	if ((scalar(@$transList)) > 1) { ($transList, $self->{ccds_only}) = $self->{genome_data}->thin_out_translist($transList); }

  	foreach my $trans(sort { $b->{translation_length} <=> $a->{translation_length}
                                                       ||
                                         $b->{length} <=> $a->{length} } @$transList){
	    my $name = $trans->gene_name;
	    my $id = $trans->gene_stable_id;
	    my $biotype = $trans->gene_biotype;
	    my $trans_id = $trans->display_id;
	    my $aa_length = $trans->translation_length;
	    if ($self->{debug}) { print "\nTranscript:  $trans_id, $name ($id)  TRANS LENGTH: " . $trans->length() . " ( aa $aa_length)\n"; }

	    # See if breakpoint is in a gene (if so, whether it is exon or intron), 5'UTR, 3'UTR, upstream or downstream of coding region
	    # check which strand the transcript is on
	    # if exon, check the phase
	    my $anno = $self->_buildRGanno($trans);
	    push @{$self->{rg_annos}}, $anno;
  	}
    }
    return($self->{rg_annos});
}

#--------------------------------------------------------------------------------------------#
# note that all coordinates are relative to the slice (ie a negative coordinate lies before the slice)
sub _getTranscripts {
    my $self = shift;

    # get the coordinates (adjusted according to the 'within' parameter)
    $self->_getEndCoords();

    # get the transcripts overlapping this region
    if ($self->{debug}) { print "\nSLICE " . $self->{chr} . ':' . $self->{pos_start} . '-' . $self->{pos_end} . "\n"; }
    my $transcripts = [];
    try {
      $transcripts = $self->{genome_data}->fetch_transcripts_by_region($self->{chr}, ($self->{pos_start}), ($self->{pos_end}));
    }
    catch {
      croak $_ if($_ !~ m/^Tabix::tabix_read: iter is not of type ti_iter_t/);
    };
    return($transcripts);
}

#--------------------------------------------------------------------------------------------#
# get the coordinates for the requested end, adjusted according to the 'within' parameter
sub _getEndCoords {
    my $self = shift;

    my ($chr, $pos_start, $pos_end);

    unless ($self->{within}) { $self->{within} = 0; }

    my $end = $self->{end};
    my $strand1 = $self->{entry}->strand1();
    my $strand2 = $self->{entry}->strand2();

    if    (($end == 1) && ($strand1 eq '+') && ($self->{entry}->pos1_end() >= $self->{entry}->pos1_start())) {
	$chr = ($self->{entry}->chr1());
	$pos_start =  $self->{entry}->pos1_start();
	$pos_end   = ($self->{entry}->pos1_end()) + ($self->{within});
    }
    elsif (($end == 2) && ($strand2 eq '+') && ($self->{entry}->pos2_end() >= $self->{entry}->pos2_start())) {
	$chr = ($self->{entry}->chr2());
	$pos_start =  $self->{entry}->pos2_start();
	$pos_end   = ($self->{entry}->pos2_end()) + ($self->{within});
    }
    elsif (($end == 1) && ($strand1 eq '-') && ($self->{entry}->pos1_end() >= $self->{entry}->pos1_start())) {
	$chr = ($self->{entry}->chr1());
	$pos_start = ($self->{entry}->pos1_start()) - ($self->{within});
	$pos_end   =  $self->{entry}->pos1_end();
    }
    elsif (($end == 2) && ($strand2 eq '-') && ($self->{entry}->pos2_end() >= $self->{entry}->pos2_start())) {
	$chr = ($self->{entry}->chr2());
	$pos_start = ($self->{entry}->pos2_start()) - ($self->{within});
	$pos_end   =  $self->{entry}->pos2_end();
    }


    # this is to deal with brassII coordinates which are sometimes flipped over
    if    (($end == 1) && ($strand1 eq '+') && ($self->{entry}->pos1_start() > $self->{entry}->pos1_end())) {
	$chr = ($self->{entry}->chr1());
	$pos_end   =  $self->{entry}->pos1_start();
	$pos_start = ($self->{entry}->pos1_end()) + ($self->{within});
    }
    elsif (($end == 2) && ($strand2 eq '+') && ($self->{entry}->pos2_start() > $self->{entry}->pos2_end())) {
	$chr = ($self->{entry}->chr2());
	$pos_end   =  $self->{entry}->pos2_start();
	$pos_start = ($self->{entry}->pos2_end()) + ($self->{within});
    }
    elsif (($end == 1) && ($strand1 eq '-') && ($self->{entry}->pos1_start() > $self->{entry}->pos1_end())) {
	$chr = ($self->{entry}->chr1());
	$pos_end   = ($self->{entry}->pos1_start()) - ($self->{within});
	$pos_start =  $self->{entry}->pos1_end();
    }
    elsif (($end == 2) && ($strand2 eq '-') && ($self->{entry}->pos2_start() > $self->{entry}->pos2_end())) {
	$chr = ($self->{entry}->chr2());
	$pos_end   = ($self->{entry}->pos2_start()) - ($self->{within});
	$pos_start =  $self->{entry}->pos2_end();
    }


    $self->{chr} = $chr;
    $self->{pos_start} = $pos_start;
    $self->{pos_end} = $pos_end;
}
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

sub _buildRGanno {
    my $self = shift;
    my $tran = shift;

    return unless ($tran);

    # convert rearrangement coordinates to slice coordinate system (-n -> 0 -> +m where 1 is the base we requested)
    my $pos_s = ($self->{within} + 1); # add 1 to shift it onto the actual base
    my $pos_e = $self->{pos_end} - $self->{pos_start} + ($self->{within} + 1);# add 1 to shift it onto the actual base

    # get the strand, region type and phase(if relevant)
    my $strand = $tran->strand();

    my ($region, $phase, $start_exon, $end_exon, $start_intron, $end_intron, $start_region, $end_region, $LHlength, $start_base, $end_base, $region_count, $upstream_seq, $downstream_seq) = $self->_getRegion($pos_s,$pos_e,$tran);

    my ($start_number, $end_number);
    if ($start_exon) { $start_number = $start_exon; }
    elsif ($start_intron) { $start_number = $start_intron; }
    if ($end_exon)   { $end_number = $end_exon; }
    elsif ($end_intron)   { $end_number = $end_intron; }

    # get various names...
    my $label = $tran->gene_name;
    my $id = $tran->gene_stable_id;
    my $biotype = $tran->gene_biotype;
    my $entrez_id = ($tran->entrez_id || '');
    my $transcript = $tran->display_id();
    my $transcript_id = $tran->stable_id();
    my $ccds = ($tran->ccds_id() || '');
    my $acc = ($tran->accession() || '');
    my $length = int($tran->length);
    my $translation_length = int($tran->translation_length);


    # fix up $phase, needs to be forced text
    if(defined $phase) {
      $phase .= q{};
    }
    else {
      $phase = q{};
    }

    # make half an RGanno object and populate it with the RGannoPoint objects (L5/L3 or H5/H3) for this end
    my $annopoint1 = undef;
    my $annopoint2 = undef;
    $annopoint1 = new Sanger::CGP::Grass::AnnoPoint(-gene  => $label,
				       -gene_id       => $id,
				       -biotype       => $biotype,
				       -transcript    => $transcript,
				       -transcript_id => $transcript_id,
				       -trans_region_count => $region_count,
				       -trans_length  => $length,
				       -translation_length => $translation_length,
				       -entrez_id     => $entrez_id,
				       -strand        => $strand,
				       -phase         => $phase,
				       -up2           => $upstream_seq,
				       -down2         => $downstream_seq,
				       -region        => $start_region,
				       -region_number => $start_number,
				       -start_base    => $start_base,
				       -end_base      => $end_base);

    $annopoint2 = new Sanger::CGP::Grass::AnnoPoint(-gene => $label,
				       -gene_id       => $id,
				       -biotype       => $biotype,
				       -transcript    => $transcript,
				       -transcript_id => $transcript_id,
				       -trans_region_count => $region_count,
				       -trans_length  => $length,
				       -translation_length => $translation_length,
				       -entrez_id     => $entrez_id,
				       -strand        => $strand,
				       -phase         => $phase,
				       -up2           => $upstream_seq,
				       -down2         => $downstream_seq,
				       -region        => $end_region,
				       -region_number => $end_number,
				       -start_base    => $start_base,
				       -end_base      => $end_base);
    my $anno = undef;
    if    ($self->{end} == 1) {
	$anno = new Sanger::CGP::Grass::Anno(-L5      => $annopoint1,
				-L3      => $annopoint2,
				-Ltype   => $region,
				-Llength => $LHlength,
				-id_rg   => ($self->{entry}->name()) );
    }
    elsif ($self->{end} == 2) {
	$anno = new Sanger::CGP::Grass::Anno(-H5      => $annopoint1,
				-H3      => $annopoint2,
				-Htype   => $region,
				-Hlength => $LHlength,
				-id_rg   => ($self->{entry}->name()) );
    }



if ($self->{debug}) {
    print "End " . $self->end() . ", ";
    print "Strand " . $tran->strand() . ", ";
    print "Region " . $region . ", ";
    print "Region(s) " . $start_region . ", ";
    print "Region(e) " . $end_region . ", ";
    if (defined($phase)) { print "Phase  " . $phase . ","; }
    if ($start_exon) { print "ExonStart " . $start_exon . ", "; }
    if ($end_exon)   { print "ExonEnd " . $end_exon . " "; }
    if ($start_intron) { print "IntronStart " . $start_intron . ", "; }
    if ($end_intron)   { print "IntronEnd " . $end_intron . " "; }
    print "\n";
}
    return($anno);
}

#--------------------------------------------------------------------------------------------#

# see if pos_start/pos_end overlap with transcript_start/transcript_end or the transcript coding region
# coordinates are all relative to the slice
sub _getRegion {
    my ($self,$pos_start,$pos_end,$tran) = @_;

    my $region = '';
    my $start_region = '';
    my $end_region = '';
    my $phase = '';
    my $start_exon = '';
    my $end_exon = '';
    my $start_intron = '';
    my $end_intron = '';
    my $start_base = 0; # whether an exact changepoint is the first base of an exon
    my $end_base = 0; # whether an exact changepoint is the last base of an exon
    my $region_count = 0;
    my $upstream_seq = '';
    my $downstream_seq = '';

    my $trans_start = $tran->start();
    my $trans_end = $tran->end();
    my $cds_start = $tran->coding_region_start();
    my $cds_end = $tran->coding_region_end();
    if ($self->{debug}) {
	if ($cds_start && $cds_end) { print $tran->display_id . " RG:$pos_start,$pos_end  TRANS:$trans_start,$trans_end CDS:$cds_start,$cds_end\n"; }
	elsif ($trans_start && $trans_end) { print $tran->display_id . " RG:$pos_start,$pos_end  TRANS:$trans_start,$trans_end CDS:\n"; }
    }
    # get the degree of overlap of pos_start and pos_end with cds, or transcript
    # and whether the start/end range is entirely within/without the cds or the transcript, or just overlapping
    my ($overlap_sc, $overlap_ec, $overlap_st, $overlap_et, $overlap_c, $span_c, $inside_c, $overlap_u, $span_u, $inside_u);

   ($overlap_sc)                                           = $self->_overlap($pos_start,$pos_start,$cds_start  ,$cds_end  , $tran->strand());
    unless ($overlap_sc)                   { ($overlap_st) = $self->_overlap($pos_start,$pos_start,$trans_start,$trans_end, $tran->strand()); }

   ($overlap_ec)                                           = $self->_overlap($pos_end  ,$pos_end  ,$cds_start  ,$cds_end  , $tran->strand());
    unless ($overlap_ec)                   { ($overlap_et) = $self->_overlap($pos_end  ,$pos_end  ,$trans_start,$trans_end, $tran->strand()); }

   ($overlap_c, $span_c, $inside_c)                        = $self->_overlap($pos_start,$pos_end  ,$cds_start  ,$cds_end  , $tran->strand());
    unless ($overlap_c)   { ($overlap_u,$span_u,$inside_u) = $self->_overlap($pos_start,$pos_end  ,$trans_start,$trans_end, $tran->strand()); }

    # is end entirely within coding region?
    if ($overlap_c && $inside_c) {
	$region = REGION_TYPE_CDS;
	$start_region = REGION_TYPE_CDS;
	$end_region = REGION_TYPE_CDS;
	if ($self->{debug}) { print "cds: $pos_start,$pos_end," . ($cds_start || '') . ',' . ($cds_end || '') . ',' . $tran->strand() . "\n"; }
	# is end entirely in intron/exon? (get which exon/intron it is in)
	($start_exon, $end_exon, $start_intron, $end_intron, $phase, $start_base, $end_base, $region_count, $upstream_seq, $downstream_seq) = $self->_getIntronExon($pos_start,$pos_end,$tran);
	if    ($start_exon   && $end_exon   && ($start_exon   == $end_exon))   { $region = REGION_TYPE_EXON; }
	elsif ($start_intron && $end_intron && ($start_intron == $end_intron)) { $region = REGION_TYPE_INTRON; }
	if ($start_exon)      { $start_region = REGION_TYPE_EXON; }
	elsif ($start_intron) { $start_region = REGION_TYPE_INTRON; }
	if ($end_exon)        { $end_region = REGION_TYPE_EXON; }
	elsif ($end_intron)   { $end_region = REGION_TYPE_INTRON; }
    }

    # is end overlapping coding region?
    elsif ($overlap_c) {
	$region = REGION_TYPE_MAYBE_CDS;
	if ($self->{debug}) { print "cds: $pos_start,$pos_end," . ($cds_start || '') . ',' . ($cds_end || '') . ',' . $tran->strand() . "\n"; }
	# check to see which introns/exons it spans
	($start_exon, $end_exon, $start_intron, $end_intron, $phase, $start_base, $end_base, $region_count, $upstream_seq, $downstream_seq) = $self->_getIntronExon($pos_start,$pos_end,$tran);
	if ($start_exon)      { $start_region = REGION_TYPE_EXON; }
	elsif ($start_intron) { $start_region = REGION_TYPE_INTRON; }
	if ($end_exon)        { $end_region = REGION_TYPE_EXON; }
	elsif ($end_intron)   { $end_region = REGION_TYPE_INTRON; }
    }

    # is end entirely within UTR region?
    elsif ($overlap_u && $inside_u) {
	if ($self->{debug}) { print "utr: $pos_start,$pos_end," . ($cds_start || '') . ',' . ($cds_end || '') . ',' . $tran->strand() . "\n"; }
	# check to see which introns/exons it spans
	($start_exon, $end_exon, $start_intron, $end_intron, $phase, $start_base, $end_base, $region_count, $upstream_seq, $downstream_seq) = $self->_getIntronExon($pos_start,$pos_end,$tran);
	# check to see if 5'/3'
	$region = REGION_TYPE_UTR;
	$start_region = REGION_TYPE_UTR;
	$end_region = REGION_TYPE_UTR;
	if  ($cds_start) {
	    if ( (($pos_start < $cds_start) && ($tran->strand() ==  1)) ||
	         (($pos_start > $cds_start) && ($tran->strand() == -1)) ) {
		$region = REGION_TYPE_5UTR;
		$start_region = REGION_TYPE_5UTR;
		$end_region = REGION_TYPE_5UTR;
		if ($start_intron && $end_intron && ($start_intron == $end_intron)) { $region = REGION_TYPE_5UTRINTRON; }
		if ($start_intron) { $start_region = REGION_TYPE_5UTRINTRON; }
		if ($end_intron)   { $end_region = REGION_TYPE_5UTRINTRON; }
		if ($start_exon && $end_exon && ($start_exon == $end_exon)) { $region = REGION_TYPE_5UTREXON; }
		if ($start_exon) { $start_region = REGION_TYPE_5UTREXON; }
		if ($end_exon)   { $end_region = REGION_TYPE_5UTREXON; }
	    }
	    elsif ( (($pos_end > $cds_end) && ($tran->strand() ==  1)) ||
	            (($pos_end < $cds_end) && ($tran->strand() == -1)) ) {
		$region = REGION_TYPE_3UTR;
		$start_region = REGION_TYPE_3UTR;
		$end_region = REGION_TYPE_3UTR;
		if ($start_intron && $end_intron && ($start_intron == $end_intron)) { $region = REGION_TYPE_3UTRINTRON; }
		if ($start_intron) { $start_region = REGION_TYPE_3UTRINTRON; }
		if ($end_intron)   { $end_region = REGION_TYPE_3UTRINTRON; }
		if ($start_exon && $end_exon && ($start_exon == $end_exon)) { $region = REGION_TYPE_3UTREXON; }
		if ($start_exon) { $start_region = REGION_TYPE_3UTREXON; }
		if ($end_exon)   { $end_region = REGION_TYPE_3UTREXON; }
	    }
	}
	else {
	    if ($start_intron && $end_intron && ($start_intron == $end_intron)) { $region = REGION_TYPE_UTRINTRON; }
	    if ($start_intron) { $start_region = REGION_TYPE_UTRINTRON; }
	    if ($end_intron)   { $end_region = REGION_TYPE_UTRINTRON; }
	    if ($start_exon && $end_exon && ($start_intron == $end_intron)) { $region = REGION_TYPE_UTREXON; }
	    if ($start_exon) { $start_region = REGION_TYPE_UTREXON; }
	    if ($end_exon)   { $end_region = REGION_TYPE_UTREXON; }
	}
    }

    # is end overlapping UTR region?
    elsif ($overlap_u) {
	if ($self->{debug}) { print "utr: $pos_start,$pos_end," . ($cds_start || '') . ',' . ($cds_end || '') . ',' . $tran->strand() . "\n"; }
	# check to see if 5'/3'
	$region = REGION_TYPE_NONCDS;
	if  ($cds_start) {
	    $region = REGION_TYPE_3NONCDS;
	    if ( (($pos_start < $cds_start) && ($tran->strand() ==  1)) ||
	         (($pos_start > $cds_start) && ($tran->strand() == -1)) ) { $region = REGION_TYPE_5NONCDS; }
	}
    }

    # else it is flanking
    else {
	if ($self->{debug}) { print "flanking: $pos_start,$pos_end," . ($cds_start || '') . ',' . ($cds_end || '') . ',' . $tran->strand() . "\n"; }
        # check to see if end is in upstream/downstream
	$region = REGION_TYPE_DOWNSTREAM;
	$start_region = REGION_TYPE_DOWNSTREAM;
	$end_region = REGION_TYPE_DOWNSTREAM;
	if    ( (($pos_start < $trans_start) && ($tran->strand() ==  1)) ||
		(($pos_start > $trans_start) && ($tran->strand() == -1)) ) { $region = REGION_TYPE_UPSTREAM;
									     $start_region = REGION_TYPE_UPSTREAM;
									     $end_region = REGION_TYPE_UPSTREAM; }
    }

    unless ($start_region) {
	if    ($overlap_sc)   { $start_region = REGION_TYPE_CDS; }
	elsif ($overlap_st &&
	       $cds_start &&
	       (  (($pos_start < $cds_start) && ($tran->strand() ==  1)) ||
		  (($pos_start > $cds_start) && ($tran->strand() == -1)) ) ) { $start_region = REGION_TYPE_5UTR; }
	elsif ($overlap_st) { $start_region = REGION_TYPE_3UTR; }
	elsif  ( (($pos_start < $trans_start) && ($tran->strand() ==  1)) ||
		 (($pos_start > $trans_start) && ($tran->strand() == -1)) )  { $start_region = REGION_TYPE_UPSTREAM; }
	else { $start_region = REGION_TYPE_DOWNSTREAM; }
    }

    unless ($end_region) {
	if ($overlap_ec)   { $end_region = REGION_TYPE_CDS; }
	elsif ($overlap_et &&
               $cds_start &&
	       ( (($pos_end < $cds_start) && ($tran->strand() ==  1)) ||
		 (($pos_end > $cds_start) && ($tran->strand() == -1)) ) ) { $end_region = REGION_TYPE_5UTR; }
	elsif ($overlap_et) { $end_region = REGION_TYPE_3UTR; }
	elsif  ( (($pos_end < $trans_start) && ($tran->strand() ==  1)) ||
		 (($pos_end > $trans_start) && ($tran->strand() == -1)) ) { $end_region = REGION_TYPE_UPSTREAM; }
	else { $end_region = REGION_TYPE_DOWNSTREAM; }
    }


    unless ($overlap_c) { $overlap_c = $overlap_u; }
    return($region, $phase, $start_exon, $end_exon, $start_intron, $end_intron, $start_region, $end_region, $overlap_c, $start_base, $end_base, $region_count, $upstream_seq, $downstream_seq);
}

#--------------------------------------------------------------------------------------------#
# find out where the breakpoint falls within a coding region (including exon or intron number)
sub _getIntronExon {
    my ($self,$pos_start,$pos_end,$tran) = @_;

    my @exons = @{($tran->exons())};
    if ($tran->strand() == -1) { @exons = reverse(@exons); } # reverse the exon array order when transcript is on the negative strand
                                                             # so we go the same direction along the chromosome
    my $exon_count = 0;
    my $intron_count = 0;
    my $start_exon = 0;
    my $start_intron = 0;
    my $end_exon = 0;
    my $end_intron = 0;
    my $phase = undef; # set this when the start intron is found (it is start phase of the next exon), or what the phase at an exon base is
    my $start_base = 0;
    my $end_base = 0;

    my $cds_start = undef;
    my $cds_end = undef;

    foreach my $exon(@exons) {
	$exon_count++;
	if ($self->{debug} && ($self->{debug} == 2)) { print "query $pos_start, exon $exon_count " . $exon->start() . "-" . $exon->end() .  " P" . $exon->phase() . ", cds_start "  . ($exon->coding_region_start() || '')  . ", cds_end "  . ($exon->coding_region_end() || '') . "\n"; }
	unless ($start_exon || $start_intron) {
	    if    ($pos_start < $exon->start()) { $start_intron = $intron_count;
						  $phase = $self->_getPhase($tran,$intron_count,\@exons); }
	    elsif ($pos_start <= $exon->end())  { $start_exon   = $exon_count;   }
	}
	unless ($end_exon || $end_intron) {
	    if    ($pos_end < $exon->start()) { $end_intron = $intron_count; }
	    elsif ($pos_end <= $exon->end())  { $end_exon   = $exon_count;   }
	}

	# for RNA stuff, with exact coordinates, the coordinate is often the very first or very last base - mark this
	if ($pos_start == $pos_end) {
	    if ($pos_start == $exon->start()) {
		if ($tran->strand() == -1) { $end_base = 1; }
		else                       { $start_base = 1; }
	    }
	    if ($pos_start == $exon->end()) {
		if ($tran->strand() == -1) { $start_base = 1; }
		else                       { $end_base = 1; }
	    }
	}

	if ($phase && $phase eq '-1') { $phase = undef; } # -1 means there isn't a phase

	$intron_count++;
    }

    # work out the phase if breakpoint is mid exon and within a coding region
    if ($self->{debug} &&
        $tran->coding_region_start() &&
	$tran->coding_region_end()) {
    }
    if (defined($tran->coding_region_start()) &&
	defined($tran->coding_region_end()) &&
	!(defined($phase)) &&
	($pos_start == $pos_end) &&
	($pos_start >= $tran->coding_region_start()) &&
        ($pos_start <= $tran->coding_region_end()) ) {
#	print "checking exon phase\n";
	$phase = $self->_getExonPhase($pos_start, $tran);
    }

    # get the coding sequence upstream and downstream of an exact breakpoint if it falls inside a coding region (exonic or intronic)
    my ($upstream_seq, $downstream_seq);
    if (defined($tran->coding_region_start()) &&
	defined($tran->coding_region_end()) &&
	($pos_start == $pos_end) &&
	($pos_start >= $tran->coding_region_start()) &&
        ($pos_start <= $tran->coding_region_end()) ) {
	($upstream_seq, $downstream_seq) = $self->_getSeqs($pos_start, $tran);
    }

    # on the negative strand, exons have to be numbered the other way, so do the conversion
    if ($tran->strand() == -1) {
	if ($start_exon)   { $start_exon   = scalar(@exons) - $start_exon + 1; }
	if ($end_exon)     { $end_exon     = scalar(@exons) - $end_exon + 1; }
	if ($start_intron) { $start_intron = scalar(@exons) - $start_intron; }
	if ($end_intron)   { $end_intron   = scalar(@exons) - $end_intron; }
    }

#   print "$pos_start $pos_end: start exon/intron  $start_exon/$start_intron. end exon/intron $end_exon/$end_intron\n";

    unless ($start_intron == $end_intron) { $phase = undef; } # reset the intron implied phase if introns aren't the same across the range

    # possible site of rearrangement breakpoint spans more than one intron/exon. Can only say it is in CDS
    #elsif ($start_exon && $end_exon) { $exons = "$start_exon-$end_exon"; }
    #elsif ($start_intron && $end_intron) { $introns = "$start_intron-$end_intron"; }

    return($start_exon, $end_exon, $start_intron, $end_intron, $phase, $start_base, $end_base, $exon_count, $upstream_seq, $downstream_seq);
}

#--------------------------------------------------------------------------------------------#

# this is the start phase for the exon AFTER the intron we are considering
#  0 = starts with complete codon, 1 = 2 bases of a codon present, 2 = 1 base of a codon present.
sub _getPhase {
    my ($self, $tran, $number, $exons) = @_;
    my $phase = undef;

    if ($tran->strand() == -1) {
	if ($exons->[$number - 1]) { $phase = $exons->[$number - 1]->phase(); }  # phase of next exon in exon array (NB array reversed)
    }
    else {
	if ($exons->[$number])     { $phase = $exons->[$number]->phase(); } # next exon in the exon array
    }

    return($phase);
}

#--------------------------------------------------------------------------------------------#
# phase of an exact coordinate within an exon
sub _getExonPhase {
    my ($self, $pos, $tran) = @_;
    my $phase = undef;

    # pos is always (1 + within) - coordinates are relative to this

    # sort out phases if they haven't been done yet
    # add up the coding bases to work out the phase
    my $coding_bases = 0;
    foreach my $exon(@{$tran->exons}) {
	if ($self->{debug}) { print "test " . ($exon->coding_region_start() || '') . " (" . ($exon->phase || '') . ") " . ($exon->coding_region_end() || '') . " (" . ($exon->end_phase || '') . ")\n"; }

	# if the breakpoint is in this exon...
	if ( ($pos >= $exon->start()) &&
	     ($pos <= $exon->end()) ) {
	    my $used_exon_bases = 0;
	    if ($tran->strand() == -1) { $used_exon_bases = $exon->coding_region_end() - $pos + 1; }
	    else                       { $used_exon_bases = $pos - $exon->coding_region_start(); }
	    if ($self->{debug}) { print "add $used_exon_bases (last)\n"; }
	    $coding_bases += $used_exon_bases;
	    last;
	}
	# if the breakpoint is Not in this exon, use it all...
	elsif ($exon->coding_region_start()) {
	    my $total_exon_bases = $exon->coding_region_end() - $exon->coding_region_start() + 1;
	    if ($self->{debug}) { print $exon->coding_region_end() . " - " . $exon->coding_region_start() . " + 1\n"; }
	    if ($self->{debug}) { print "add $total_exon_bases\n"; }
	    $coding_bases += $total_exon_bases;
	}
    }
    # work out phase
    if    (($self->{end} == 1) && ($tran->strand() == -1)) { $phase = ($coding_bases)%3; } # want the end phase (phase of the following base) not the phase
    elsif (($self->{end} == 1) && ($tran->strand() ==  1)) { $phase = ($coding_bases + 1)%3; } # want the end phase (phase of the following base) not the phase
    elsif (($self->{end} == 2) && ($tran->strand() == -1)) { $phase = ($coding_bases - 1)%3; } # want the phase of the current base
    elsif (($self->{end} == 2) && ($tran->strand() ==  1)) { $phase = ($coding_bases )%3; } # want the phase of the current base
    if ($self->{debug}) { print "coding bases $coding_bases phase $phase\n"; }

    return($phase);
}

#--------------------------------------------------------------------------------------------#
# phase of an exact coordinate within an exon
sub _getSeqs {
    my ($self, $pos, $tran) = @_;

    # only get here if the breakpoint falls within a coding region
    # pos is always (1 + within) - coordinates are relative to this

    # get the coding sequence around the breakpoint
    # return 2 bases either side

    my $upstream_seq = '';
    my $downstream_seq = '';
    my $upstream_check = 1;
    my $downstream_check = 0;

    foreach my $exon(@{$tran->exons}) {
	if ($self->{debug}) { print "test " . ($exon->coding_region_start() || '') . " (" . ($exon->phase || '') . ") " . ($exon->coding_region_end() || '') . " (" . ($exon->end_phase || '') . ")\n"; }

	# switch to downstream once the breakpoint is found
	if (($tran->strand() == 1) && $exon->coding_region_start() && ($exon->coding_region_start() >= $pos)) {
	    $upstream_check = 0;
	    $downstream_check = 1;
	}
	elsif (($tran->strand() == -1) && $exon->coding_region_start() && ($exon->coding_region_start() < $pos)) {
	    $upstream_check = 0;
	    $downstream_check = 1;
	}

	my ($len, $start, $start_rev);
	if (defined($exon->coding_region_start()) && defined($exon->coding_region_end())) {
	    $len = $exon->coding_region_end() - $exon->coding_region_start() + 1; # length of entire coding region
	    $start = $exon->coding_region_start() - $exon->start; # number of first coding base in this exon counting from the exon start
	    $start_rev = $exon->end - $exon->coding_region_end(); # number of first coding base in this exon counting back from the exon end
	}
	# if the exon contains the breakpoint, part is upstream, part is downstream (pos is always (1 + within) - slice is relative to this)
	if ( ($pos >= $exon->start()) &&
	     ($pos <= $exon->end()) ) {

	    if ($tran->strand() == -1) {
		my $start_u =  $exon->end() - $exon->coding_region_end();
		my $len_u = $exon->coding_region_end() - $pos;
		my $start_d = $exon->end() - $pos;
		my $len_d = $pos -($exon->coding_region_start()) + 1;
		if ($self->{end} == 1) { # for end1, include the base we are on in the upstream region, for end2 its included in the downstream region
		    $len_u++;
		    $len_d--;
		    $start_d++;
		}
		# upstream portion...
		$upstream_seq   .= substr($exon->seq,$start_u,$len_u);
		if ($self->{debug}) { print "U $len_u $start_u " . substr($exon->seq,$start_u,$len_u) . "\n"; }
		# downstream portion...
		$downstream_seq .= substr($exon->seq,$start_d,$len_d);
		if ($self->{debug}) { print "D $len_d $start_d " . substr($exon->seq,$start_d,$len_d) . "\n"; }
	    }
	    else {
	        my $len_u = $pos -($exon->coding_region_start());
		my $len_d = $exon->coding_region_end() - $pos + 1;
		my $start_d = $pos -($exon->start());
		if ($self->{end} == 1) { # for end1, include the base we are on in the upstream region, for end2 its included in the downstream region
		    $len_u++;
		    $len_d--;
		    $start_d++;
		}
		# upstream portion...
		$upstream_seq   .= substr($exon->seq,$start,$len_u);
		if ($self->{debug}) { print "U  " . substr($exon->seq,$start,$len_u) . "\n"; }
		# downstream portion...
		$downstream_seq   .= substr($exon->seq,$start_d,$len_d);
		if ($self->{debug}) { print "D  " . -($exon->start()) . "  " . substr($exon->seq,$start_d,$len_d) . "\n"; }
	    }
	}
	# if the exon is upstream of the breakpoint...
	elsif ($upstream_check && defined($len))   {
	    if ($tran->strand() == -1) {
		$upstream_seq   .= substr($exon->seq,$start_rev,$len);
		if ($self->{debug}) { print "UP $len $start_rev " . substr($exon->seq,$start_rev,$len) . "\n"; }
	    }
	    else {
		$upstream_seq   .= substr($exon->seq,$start,$len);
		if ($self->{debug}) { print "UP $len $start " . substr($exon->seq,$start,$len) . "\n"; }
	    }
	}
	# if the exon is downstream of the breakpoint...
	elsif ($downstream_check && defined($len)) {
	    if ($tran->strand() == -1) {
		$downstream_seq .= substr($exon->seq,$start_rev,$len);
		if ($self->{debug}) { print "DOWN $len $start_rev " . substr($exon->seq,$start_rev,$len) . "\n"; }
	    }
	    else {
		$downstream_seq .= substr($exon->seq,$start,$len);
		if ($self->{debug}) { print "DOWN $len $start " . substr($exon->seq,$start,$len) . "\n"; }
	    }
	}
	else { if ($self->{debug}) { print "Non-coding\n"; } }
    }
    my ($upstream_seq_pair) = ($upstream_seq =~ /(\D\D)$/);
    my ($downstream_seq_pair) = ($downstream_seq =~ /^(\D\D)/);

    return($upstream_seq_pair,$downstream_seq_pair);
}

#--------------------------------------------------------------------------------------------#

sub _overlap {
    my ($self, $pos1_s, $pos1_e, $pos2_s, $pos2_e, $strand) = @_;
    my $overlap = 0;
    my $span = 0;
    my $inside = 0;

    return unless (defined($pos1_s) && defined($pos1_e) && defined($pos2_s) &&  defined($pos2_e));

    my ($max1,$min1, $max2,$min2);
    if ($pos1_e >= $pos1_s) { $max1 = $pos1_e; $min1 = $pos1_s; }
    else                    { $max1 = $pos1_s; $min1 = $pos1_e; }
    if ($pos2_e >= $pos2_s) { $max2 = $pos2_e; $min2 = $pos2_s; }
    else                    { $max2 = $pos2_s; $min2 = $pos2_e; }

    if (($max1 >= $min2) && ($max2 >= $min1)) { $overlap = 1; } # check to see if they overlap
    if (($min1 <= $min2) && ($max1 >= $max2)) { $span = 1; }    # check to see if 1 completely spans 2
    if (($min2 <= $min1) && ($max2 >= $max1)) { $inside = 1; }  # check to see if 1 is completely inside 2

    # get the amount of overlap
    if ($overlap) {
	if (($self->end()) == 1) {
	    my $length = $max1 - $min2;
	    if ($length > 0) { $overlap = $length; }
	}
	else {
	    my $length = $max2 - $min1;
	    if ($length > 0) { $overlap = $length; }
	}
	if ($self->{debug}) { print " RG: $pos1_s, $pos1_e   CDS/trans: $pos2_s, $pos2_e   strand: $strand    overlap: $overlap\n"; }
    }

    return($overlap,$span,$inside);
}
#--------------------------------------------------------------------------------------------#


1;
