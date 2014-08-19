package Sanger::CGP::Grass::Annotation::RGannotator;

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
use Sanger::CGP::Grass::Annotation::RGendAnnotator;
use Sanger::CGP::Grass::Annotation::RGcombinationAnnotator;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

=head1 NAME

RGannotator

=head1 SYNOPSIS

use Sanger::CGP::Grass::Annotation::RGannotator;

my $rgann = new Sanger::CGP::Grass::Annotation::RGannotator(-dataset         => $dataset,
						-within          => $within,
						-genome_data    => $genome_data_object,
						-list_between    => 1,
						-show_biotype    => 1 );

$rgann->getAnnotation();

my $output_string = $rgann->format_for_printing();

=head1 DESCRIPTION

Class to control getting gene annotations for a pair of breakpoints and predict possible fusion outcomes.

Takes in a set of DataEntry objects.

Returns a set of completed annotations, formated for printing.


=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::Annotation::RGannotator();
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
    if ($args{-within})   { $self->within($args{-within}); }
    if ($args{-genome_data}) { $self->genome_data($args{-genome_data}); }
    if (defined($args{-entrez_required})) { $self->entrez_required($args{-entrez_required}); }
    if (defined($args{-gene_id_required})) { $self->gene_id_required($args{-gene_id_required}); }
    if (defined($args{-show_biotype})) { $self->show_biotype($args{-show_biotype}); }
    if (defined($args{-list_between})) { $self->list_between($args{-list_between}); }
    if (defined($args{-multi_annos})) { $self->multi_annos($args{-multi_annos}); }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 dataset

  Arg (1)    : $dataset
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

=head2 entry

  Arg (1)    : $dataEntry
  Example    : $entry = $object->entry($dataEntry);
  Description: Data entry object containing coordinates (acts as holder for coordinates)
  Return     : $dataEntry

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

  Arg (1)    : $registry
  Example    : $registry = $Object->registry($registry);
  Description: get/set the Ensembl registry object to use
  Return     : registry object

=cut

sub registry {
    my $self = shift;
    $self->{registry} = shift if @_;
    return $self->{registry};
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

  Arg (1)    : 1/0
  Example    : $entrez_required = $Object->entrez_required($entrez_required);
  Description: define whether entrez id is required
  Return     : 1/0

=cut

sub entrez_required {
    my $self = shift;
    $self->{entrez_required} = shift if @_;
    return $self->{entrez_required};
}
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

=head2 show_biotype

  Arg (1)    : 1/0
  Example    : $show_biotype = $Object->show_biotype($show_biotype);
  Description: whether to show the biotype or not (eg protein_coding) for each gene
  Return     : 1/0

=cut

sub show_biotype {
    my $self = shift;
    $self->{show_biotype} = shift if @_;
    return $self->{show_biotype};
}
#-----------------------------------------------------------------------#

=head2 list_between

  Arg (1)    : 1/0
  Example    : $list_between = $Object->list_between($list_between);
  Description: define whether to include a list of all the genes that lie between the 2 breakpoints
  Return     : 1/0

=cut

sub list_between {
    my $self = shift;
    $self->{list_between} = shift if @_;
    return $self->{list_between};
}
#-----------------------------------------------------------------------#

=head2 multi_annos

  Arg (1)    : 1/0
  Example    : $multi_annos = $Object->multi_annos($multi_annos);
  Description: get more than one annotation if there are other potentially interesting ones
  Return     : 1/0

=cut

sub multi_annos {
    my $self = shift;
    $self->{multi_annos} = shift if @_;
    return $self->{multi_annos};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 rg_annos

  Arg (1)    : $rg_annos
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
	my ($anns1,$ccds1) = $self->_get_end_anns($entry, 1);
	my ($anns2,$ccds2) = $self->_get_end_anns($entry, 2);

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

		my $rganno = $self->_combine($ann1,$ann2, $entry);
		next unless ($rganno);
		push @$annos, $rganno;
	    }
	}

	# if there are more than 16 combinations (4 X 4), thin them out to keep only the top 2 scores
#	if ((scalar(@$annos) > 16) && !($ccds)) { $annos = $self->thin_out_annos($annos); }

	# keep only 1 of each flag value (longest total genomic length covered by transcript)
	if ($annos && scalar(@$annos)) {
	    $annos = $self->_thin_out_annos2($annos);
	}
	# if there are annotations for only end 1
	elsif ($anns1 && scalar(@$anns1) && (!$anns2 || !(scalar(@$anns2))))  {
	    $annos = $self->_thin_end($anns1);
	}
	# if there are annotations for only end 2
	elsif ($anns2 && scalar(@$anns2) && (!$anns1 || !(scalar(@$anns1))))  {
	    $annos = $self->_thin_end($anns2);
	}

        # just keep the top annotation unless more are requested
	if ($self->{multi_annos}) {
	    foreach (@$annos) { push @{$self->{rg_annos}}, $_; }
	}
	else {
	    foreach (@$annos) { push @{$self->{rg_annos}}, $_; last; }
	}
	$count++;
    }
    if ($self->{list_between}) { $self->get_list_between(); }
}
#-----------------------------------------------------------------------#
sub _get_end_anns {
    my $self = shift;
    my $entry = shift;
    my $end = shift;

    my $rgend = new Sanger::CGP::Grass::Annotation::RGendAnnotator(-entry    => $entry,
						      -end      => $end,
						      -registry => $self->{registry},
						      -genome_data => $self->{genome_data},
						      -entrez_required => $self->{entrez_required},
						      -gene_id_required => $self->{gene_id_required},
						      -within   => $self->{within});
    my $anns = $rgend->annotate();

    return($anns, ($rgend->ccds_only()));
}
#---------------------------------------------------------------------------------------------------------------#
sub _combine {
    my $self = shift;
    my $ann1 = shift;
    my $ann2 = shift;
    my $entry = shift;

    # returns an RGanno (fusion) object (keep data for all pairs, even flag = 0)
    my $combi = new Sanger::CGP::Grass::Annotation::RGcombinationAnnotator(-anno1   => $ann1,
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
sub _thin_out_annos {
    my $self = shift;
    my $annos = shift;
    my @new_annos = ();

    my $flag_count = 0;
    my $last_flag = 0;
    foreach my $anno(sort {$b->id_fusion_flag <=> $a->id_fusion_flag} @$annos) {
	my $current_flag = $anno->id_fusion_flag();
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
sub _thin_out_annos2 {
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
	my $current_flag = $anno->id_fusion_flag();

#	print "" . $anno->L5->transcript . " $current_flag\n";
#	print "   Llength: " . ($anno->Llength() || 0) . "\n";
#	print "   Hlength: " . ($anno->Hlength() || 0) . "\n";
#	print  "   L5 trans length: " . ($anno->L5->trans_length() || 0) . "\n";
#	print  "   H5 trans length: " . ($anno->H5->trans_length() || 0) . "\n";

	my $current_length = ($anno->Llength() || 0) + ($anno->Hlength() || 0); # putative fusion translation product length
	my $current_trans_length = ($anno->L5->trans_length() || 0) + ($anno->H5->trans_length() || 0); # transcript length
	if ($self->{debug}) { print "    current $current_flag $current_length trans_length $current_trans_length\n"; }

	# save the data for the previous flag value
	if ($last_flag && !($current_flag eq $last_flag)) {

	    # save the best anno for the previous flag to an array
	    if ($self->{debug}) { print "saving $last_flag $longest_length\n"; }
	    if ($best_anno4flag) { push @new_annos, $best_anno4flag; }

	    # reset the variables for the next flag to be considered
	    $longest_length = 0;
	    $longest_trans_length = 0;
	    $best_anno4flag = '';
	}

	# check putative fusion translation product length
	if ($current_length > $longest_length) {
	    if ($self->{debug}) { print "updating length $longest_length to $current_length\n"; }
#	    print "        this overlap is longer\n";
	    $longest_length = $current_length;
	    $longest_trans_length = $current_trans_length;
	    $best_anno4flag = $anno;
	}
	# check transcript length
	elsif (($current_length == $longest_length) && ($current_trans_length > $longest_trans_length)) {
	    if ($self->{debug}) { print "updating total trans_length $longest_trans_length to $current_trans_length\n"; }
#	    print "        this overlap is same but transcript is longer (current $current_trans_length longest $longest_trans_length)\n";
	    $longest_trans_length = $current_trans_length;
	    $best_anno4flag = $anno;
	}

	unless ($best_anno4flag) {
	    $best_anno4flag = $anno;
	    $longest_length = $current_length;
	    $longest_trans_length = $current_trans_length;
#	    print "setting it for the first time\n";
	}

	$last_flag = $current_flag;
    }
#    print "picked " . $best_anno4flag->L5->transcript . "\n";

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
# if only 1 end has annotations, keep only the longest upstream/downstream, or translation product, or transcript
sub _thin_end {
    my $self = shift;
    my $annos = shift;
    my @new_annos = ();

    my $longest = 0;
    my $longest_trans = 0;
    my $longest_upstream = 0;
    my $best_anno;
    foreach my $anno(@$annos) {
 	my $current_length = 0;
	my $current_trans_length = 0;

	# check the putative upstream/downstream translation product length
	my $current_upstream_length = ($anno->Llength() || 0) + ($anno->Hlength() || 0);
	# check the translation product length
	if ($anno->L5) { $current_length = ($anno->L5->translation_length() || 0); }
	if ($anno->H5) { $current_length = ($anno->H5->translation_length() || 0); }
	# check the transcript length
	if ($anno->L5) { $current_trans_length = ($anno->L5->trans_length() || 0); }
	if ($anno->H5) { $current_trans_length = ($anno->H5->trans_length() || 0); }

	# check upstream/downstream translation product length and if there is one, and it is the longest, pick this anno
	if ($current_upstream_length && ($current_upstream_length > $longest_upstream)) {
	    $best_anno = $anno;
	    $longest = $current_length;
	    $longest_trans = $current_trans_length;
	    $longest_upstream = $current_upstream_length;
	}
	# check translation product length and if there is one, and it is the longest, pick this anno
	elsif (!$current_upstream_length && !$longest_upstream && $current_length && ($current_length > $longest)) {
	    $best_anno = $anno;
	    $longest = $current_length;
	    $longest_trans = $current_trans_length;
	    $longest_upstream = $current_upstream_length;
	}
	# check translation product length and if there is None, check the transcript length instead and pick based on that
	elsif (!$current_upstream_length && !$longest_upstream && !$current_length && !$longest && $current_trans_length  && ($current_trans_length > $longest_trans)) {
	    $best_anno = $anno;
	    $longest_trans = $current_trans_length;
	    $longest_upstream = $current_upstream_length;
	}
	# if upstream products are the same length, use the longest translation product
	elsif ($current_upstream_length && $longest_upstream && ($current_upstream_length == $longest_upstream) && ($current_length > $longest)) {
	    $best_anno = $anno;
	    $longest = $current_length;
	    $longest_trans = $current_trans_length;
	}
	# if upstream products and translation products are the same length, use the longest transcript
	elsif ($current_upstream_length && $longest_upstream && ($current_upstream_length == $longest_upstream) &&
               $current_length && $longest && ($current_length == $longest) && ($current_trans_length > $longest_trans)) {
	    $best_anno = $anno;
	    $longest_trans = $current_trans_length;
	}
	# set if currently unset
	unless ($best_anno) {
	    $best_anno = $anno;
	    $longest = $current_length;
	    $longest_trans = $current_trans_length;
	    $longest_upstream = $current_upstream_length;
	}

	if ($self->{debug}) { print "$current_upstream_length $current_length $current_trans_length. use $longest_upstream/$longest/$longest_trans\n"; }
    }
    push @new_annos, $best_anno;
    if (@new_annos) { return (\@new_annos); }
    else { print "thin end failed\n"; return ($annos); }
}
#------------------------------------------------------------------------------------------#

=head2 get_list_between

  Arg (0)    :
  Example    : $Object->get_list_between();
  Description: list of genes that lie between the coordinates in this dataset (only available where distance between the coordinates =< 1000000 bases)
  Return     : reference to an array of gene_names

=cut

sub get_list_between {
    my ($self) = @_;
    my $dataset = $self->{dataset};
    my $within = $self->{within};

    my $gene_names = '';

    # only look for genes if this is a single chromosome deletion/insertion of < 1mb
    if (($dataset->[0]->chr1 eq $dataset->[0]->chr2) && ($dataset->[0]->strand1 eq $dataset->[0]->strand2)) {
	my $distance = 0;
	my $chr = $dataset->[0]->chr1;
	my $start_coord = 0;
	my $end_coord = 0;
	if ($dataset->[0]->pos1_end < $dataset->[0]->pos2_start) {
	    $distance = $dataset->[0]->pos2_start - $dataset->[0]->pos1_end;
	    $start_coord = $dataset->[0]->pos1_end;
	    $end_coord = $dataset->[0]->pos2_start;
	}
	elsif ($dataset->[0]->pos2_end < $dataset->[0]->pos1_start) {
	    $distance = $dataset->[0]->pos1_start - $dataset->[0]->pos2_end;
	    $start_coord = $dataset->[0]->pos2_end;
	    $end_coord = $dataset->[0]->pos1_start;
	}
	return('') if ($distance > 1000000); # only get list of genes if less than 1MB between coordinates

	$gene_names = $self->{genome_data}->get_gene_list($chr, $start_coord, $end_coord);

    }
    $self->{between_genes} = $gene_names;

    return($gene_names);
}
#---------------------------------------------------------------------------------------------------------------#

=head2 format_for_printing

  Arg (0)    :
  Example    : $Object->format_for_printing();
  Description: Format the data processed by this class for printing
  Return     : formatted_string

=cut


sub format_for_printing {
    my ($self) = @_;

    my $show_biotype = $self->{show_biotype};
    my $list_between = $self->{list_between};
    my $between_genes = $self->{between_genes};
    my $entrez_required = $self->{entrez_required};
    my $string = '';
    my $found = 0;

    foreach my $result (@{$self->{rg_annos}}) {
	$found = 1;
	if ($result->L5) {
	    my $first_last_L = '';
	    if ($result->L5->start_base) { $first_last_L = 'first_base'; }
	    if ($result->L5->end_base)   { $first_last_L = 'last_base';  }
	    my $phaseL = $result->L5->phase;
	    unless (defined($phaseL)) { $phaseL = '_'; }
	    $string .= $result->L5->gene . "\t"
		. $result->L5->gene_id . "\t"
		. $result->L5->transcript_id . "\t"
		. $result->L5->strand . "\t"
		. $phaseL . "\t"
		. ($result->Ltype || $result->L5->region || '_') . "\t"
		. ($result->L5->region_number || '_') . "\t"
		. ($result->L5->trans_region_count || '_') . "\t"
		. ($first_last_L || '_') . "\t";
	    if ($show_biotype) {  $string .= '' . ($result->L5->biotype || '_')  . "\t"; }
	}
	else {
	    $string .= "_\t_\t_\t_\t_\t_\t_\t_\t_\t";
	    if ($show_biotype) {  $string .= "_\t"; }
	}

	if ($result->H5) {
	    my $first_last_H = '';
	    if ($result->H5->start_base) { $first_last_H = 'first_base'; }
	    if ($result->H5->end_base)   { $first_last_H = 'last_base'; }
	    my $phaseH = $result->H5->phase;
	    unless (defined($phaseH)) { $phaseH = '_'; }

	    $string .= $result->H5->gene . "\t"
		. $result->H5->gene_id . "\t"
		. $result->H5->transcript_id . "\t"
		. $result->H5->strand . "\t"
		. $phaseH . "\t"
		. ($result->Htype || $result->H5->region || '_') . "\t"
		. ($result->H5->region_number || '_') . "\t"
		. ($result->H5->trans_region_count || '_') . "\t"
		. ($first_last_H || '_') . "\t";
	    if ($show_biotype) {  $string .= '' . ($result->H5->biotype || '_') . "\t"; }
	}
	else {
	    $string .= "_\t_\t_\t_\t_\t_\t_\t_\t_\t";
	    if ($show_biotype) {  $string .= "_\t"; }
	}

	$string .= $result->id_fusion_flag || 0;
	if ($list_between) { $string .= "\t$between_genes";}
	$string .= "\n";
    }

    unless ($found) {
	if ($show_biotype) {  $string .= "_\t_\t"; }
	if ($list_between) { $string .= "_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t0\t$between_genes\n"; }
	else               { $string .= "_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t_\t0\n"; }
    }

    return($string);
}
#------------------------------------------------------------------------------------------------#
1;
