## Grass::Annotation::RGannotator

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

Class to control getting gene annotations for a defined breakpoint and predicting possible fusion outcomes

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Grass::Annotation::RGannotator;

use Grass::Annotation::RGendAnnotator;
use Grass::Annotation::RGcombinationAnnotator;

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    : 
  Example    : $object = new Grass::Annotation::RGannotator();
  Description: make a new RGannotator object
  Return     : RGannotator object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    if ($args{-dataset})  { $self->dataset($args{-dataset}); }
    if ($args{-entry})    { $self->entry($args{-entry}); }
    if ($args{-microhom}) { $self->microhom($args{-microhom}); }
    if ($args{-within})   { $self->within($args{-within}); }
    if ($args{-species})  { $self->species($args{-species}); }
    if ($args{-registry}) { $self->registry($args{-registry}); }
    if ($args{-entrez_required}) { $self->{entrez_required} = 1; }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 dataset

  Arg (0)    : $dataset
  Example    : $dataset = $Object->dataset($dataset);
  Description: array of dataEntry objects
  Return     : array of dataEntry objects

=cut

sub dataset {
    my $self = shift;
    $self->{dataset} = shift if @_;
    return $self->{dataset};
}
#-----------------------------------------------------------------------#

=head2 microhom

  Arg (0)    : $microhom
  Example    : $microhom = $pair->microhom($microhom);
  Description: microhom between breakpoints
  Return     : $microhom

=cut

sub microhom {
    my $self = shift;
    $self->{microhom} = shift if @_;
    return $self->{microhom};
}
#-----------------------------------------------------------------------#

=head2 entry

  Arg (0)    : $entry
  Example    : $entry = $pair->entry($entry);
  Description: entry object containing coordinates (does not exist in database - acts as holder for coordinates)
  Return     : $entry

=cut

sub entry {
    my $self = shift;
    my $entry = shift if @_;
    if ($entry) {
	push @{$self->{dataset}}, $entry;
	$self->{entry} = $entry;
    }
    return $self->{entry};
}
#-----------------------------------------------------------------------#

=head2 registry

  Arg (0)    : $registry
  Example    : $registry = $Object->registry($registry);
  Description: get/set the ensembl registry object to use
  Return     : registry object

=cut

sub registry {
    my $self = shift;
    $self->{registry} = shift if @_;
    return $self->{registry};
}
#-----------------------------------------------------------------------#

=head2 within

  Arg (0)    : $within
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

=head2 rg_annos

  Arg (0)    : $rg_annos
  Example    : $rg_annos = $Object->rg_annos($rg_annos);
  Description: array of all rg_anno (fusion) objects for this annotation set
  Return     : array of rg_anno objects

=cut

sub rg_annos {
    my $self = shift;
    $self->{rg_annos} = shift if @_;
    return $self->{rg_annos};
}
#-----------------------------------------------------------------------#

=head2 getAnnotation

  Arg (0)    : 
  Example    : $Object->getAnnotation();
  Description: Takes the set of entries and looks for all possible fusion genes.
               Call rg_annos to get the results.
  Return     : 

=cut


sub getAnnotation {
    my $self = shift;

    $self->{rg_annos} = undef;

    unless ($self->{within}) { $self->{within} = 0; }

    # go through dataset and get the annotations for each rearrangement end
    my $count = 0;
    foreach my $entry(@{$self->{dataset}}) {

	if ($self->{debug}) { print "entry " . ($entry->name || '') . " " . $entry->chr1 . " " . $entry->pos1_start . " " .$entry->strand1 . " " .$entry->chr2 . " " .$entry->pos2_start . " " .$entry->strand2 . " " . ($entry->shard || '') . "\n"; }
	my ($anns1,$ccds1) = $self->get_end_anns($entry, 1);
	my ($anns2,$ccds2) = $self->get_end_anns($entry, 2);
	
	# see if either end is limited to ccds entries only
	my $ccds = 0;
	if ($ccds1 || $ccds2) { $ccds = 1; }
	
	
	# check through the combinations to see if there are any, same gene, same transcript combinations available
	my $same_gene_trans = {};
	foreach my $ann1 (@$anns1) { 
	    foreach my $ann2 (@$anns2) {
		# if the gene is the same, and the transcript is the same
		# (should they both be in the intron?)
		if (
		    ($ann1->L5->gene_id eq $ann1->L3->gene_id) && 
		    ($ann2->H5->gene_id eq $ann2->H3->gene_id) && 
		    ($ann1->L5->gene_id eq $ann2->H5->gene_id) && 
		    ($ann1->L5->transcript_id eq $ann2->H5->transcript_id)) {
		    my $gene_id = $ann1->L5->gene_id;
		    $same_gene_trans->{$gene_id} = 1;
		}
	    }
	}
	
	# for each rearrangement, send each combination of end annotations to the RGcombinationAnnotator 
	my $annos = [];
	foreach my $ann1 (@$anns1) { 
	    foreach my $ann2 (@$anns2) {
		
		# if the gene is the same, only try to combine same transcripts
		my $gene_id = $ann1->L5->gene_id;
		next if ( ($same_gene_trans->{$gene_id}) && 
			  ($ann1->L5->gene_id eq $ann2->H5->gene_id) && 
			 !($ann1->L5->transcript_id eq $ann2->H5->transcript_id));
		
		my $rganno = $self->combine($ann1,$ann2, $entry);
		next unless ($rganno);
		push @$annos, $rganno;
	    }
	}
	
	# if there are more than 16 combinations (4 X 4), thin them out to keep only the top 2 scores 
#	if ((scalar(@$annos) > 16) && !($ccds)) { $annos = $self->thin_out_annos($annos); }
	
	# keep only 1 of each flag value (longest total genomic length covered by transcript)
	if ($annos && scalar(@$annos)) { 
	    $annos = $self->thin_out_annos2($annos); 
	}
	# if there are annotations for only end 1
	elsif ($anns1 && scalar(@$anns1) && (!$anns2 || !(scalar(@$anns2))))  {
	    $annos = $self->thin_end($anns1);
	}
	# if there are annotations for only end 2
	elsif ($anns2 && scalar(@$anns2) && (!$anns1 || !(scalar(@$anns1))))  {
	    $annos = $self->thin_end($anns2);
	}

	foreach (@$annos) { push @{$self->{rg_annos}}, $_; }

	$count++;
    }
}
#-----------------------------------------------------------------------#
sub get_end_anns {
    my $self = shift;
    my $entry = shift;
    my $end = shift;

    my $rgend = new Grass::Annotation::RGendAnnotator(-entry    => $entry,
						      -end      => $end,
						      -species  => $self->{species},
						      -registry => $self->{registry},
						      -entrez_required => $self->{entrez_required},
						      -within   => $self->{within});
    my $anns = $rgend->annotate();

    return($anns, ($rgend->ccds_only()))
}
#---------------------------------------------------------------------------------------------------------------#
sub combine {
    my $self = shift;
    my $ann1 = shift;
    my $ann2 = shift;
    my $entry = shift;

    # returns an RGanno (fusion) object (keep data for all pairs, even flag = 0)
    my $combi = new Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1,
							      -anno2   => $ann2,
							      -strand1 => $entry->strand1,
							      -strand2 => $entry->strand2,
							      -shard   => $entry->shard);
    $combi->combine(); 
    my $rganno = $combi->anno();
    return($rganno);
}
#---------------------------------------------------------------------------------------------------------------#
# only keep the top 2 highest value flags 
sub thin_out_annos {
    my $self = shift;
    my $annos = shift;
    my @new_annos = ();

    my $flag_count = 0;
    my $last_flag = 0;
    foreach my $anno(sort {$b->id_fusion_flag <=> $a->id_fusion_flag} @$annos) {
	$current_flag = $anno->id_fusion_flag();
	unless ($current_flag eq $last_flag) {
	    $flag_count++;
	    $last_flag = $current_flag;
	}
	if ($flag_count > 2) { last; }
	else {
	    if ($self->{debug}) { print $anno->id_fusion_flag() . "\n"; }
	    push @new_annos, $anno;
	}
    }
    return (\@new_annos);
}
#---------------------------------------------------------------------------------------------------------------#
# only keep the longest of each flag value 
sub thin_out_annos2 {
    my $self = shift;
    my $annos = shift;
    my @new_annos = ();

    my $last_flag = 0;
    my $longest_length = 0;
    my $longest_trans_length = 0;
    my $best_anno4flag = '';

    # go through each flag type in turn, highest to lowest
    if ($self->{debug}) { print "thinning " . scalar(@$annos) . " annos\n"; }
    foreach my $anno(sort {$b->id_fusion_flag <=> $a->id_fusion_flag} @$annos) {
	if ($self->{debug}) { print "|" . $anno->L5->transcript . "| |" . $anno->H5->transcript . "| |" . $anno->Llength  . "| |" . $anno->Hlength . "|\n"; }
	$current_flag = $anno->id_fusion_flag();
	$current_length = ($anno->Llength() || 0) + ($anno->Hlength() || 0);
	$current_trans_length = ($anno->L5->trans_length() || 0) + ($anno->H5->trans_length() || 0);
	if ($self->{debug}) { print "    current $current_flag $current_length\n"; }

	# save the data for the previous flag value
	if ($last_flag && !($current_flag eq $last_flag)) {

	    # save the best anno for the previous flag to an array
	    if ($self->{debug}) { print "saving $last_flag $longest_length\n"; }
	    if ($best_anno4flag) { push @new_annos, $best_anno4flag; }

	    # reset the variables for the next flag to be considered
	    $longest_length = 0;
	    $best_anno4flag = '';
	}

	if ($current_length > $longest_length) {
	    $longest_length = $current_length;
	    $best_anno4flag = $anno;
	    if ($self->{debug}) { print "updating length $longest_length to $current_length\n"; }
	}
	elsif (($current_length == $longest_length) && ($current_trans_length > $longest_trans_length)) {
	    $longest_trans_length = $current_trans_length;
	    $best_anno4flag = $anno;
	    if ($self->{debug}) { print "updating total trans_length $longest_trans_length to $current_trans_length\n"; }
	}

	unless ($best_anno4flag) { 
	    $best_anno4flag = $anno; 
	    $longest_length = $current_length; 
	}

	$last_flag = $current_flag;
    }

    # save the data for the last flag value
    if ($self->{debug}) { print "saving $last_flag $longest_length\n"; }
    if ($best_anno4flag) { push @new_annos, $best_anno4flag; }

    if ($self->{debug}) { print "thinned " . scalar(@new_annos) . " new annos\n"; }

    if ($self->{debug}) {
	foreach my $anno(@new_annos) {
	    unless ($anno && $anno->L5) { next; }
	    print "|" . $anno->id_rg . "|\n";
	}
    }

#    if (scalar(@new_annos) > 2) {exit;}
    if (@new_annos) { return (\@new_annos); }
    else { print "thinning failed\n"; return ($annos); }
}
#---------------------------------------------------------------------------------------------------------------#
# if only 1 end has annotations, keep only the longest 
sub thin_end {
    my $self = shift;
    my $annos = shift;
    my @new_annos = ();

    my $longest = 0;
    my $longest_trans = 0;
    my $best_anno;
    foreach my $anno(@$annos) {
	$current_length = ($anno->Llength() || 0) + ($anno->Hlength() || 0);
	$current_trans_length = 0;
	if ($anno->L5) { $current_trans_length = ($anno->L5->trans_length() || 0); }
	if ($anno->H5) { $current_trans_length = ($anno->H5->trans_length() || 0); }
	if ($current_length && ($current_length > $longest)) {
	    $best_anno = $anno;
	}
	elsif (!$longest && $current_trans_length && (($current_trans_length > $longest_trans))) {
	    $best_anno = $anno;
	}
    }
    push @new_annos, $best_anno;
    if (@new_annos) { return (\@new_annos); }
    else { print "thin end failed\n"; return ($annos); }
}
#---------------------------------------------------------------------------------------------------------------#
1;
