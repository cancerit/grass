package Sanger::CGP::Grass::Annotation::RGcombinationAnnotator;

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
use Carp;
use Sanger::CGP::Grass::Annotation::RGendAnnotator;
use Sanger::CGP::Grass::Anno;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

# list of fusion flags
use constant FF_IFF_DIFF_E_SO_RFOK => '910'; # In frame fusion gene, 2 different genes, broken in coding exons. Same orientation. RF preserved
use constant FF_IFF_DIFF_I_SO_RFOK => '900'; # In frame fusion gene, 2 different genes, broken in coding introns. Same ori. Reading frame preserved
use constant FF_IFF_SAME_E_SO_RFOK => '890'; # In frame fusion gene, same gene, broken in coding exons. Same orientation. RF preserved
use constant FF_IFF_SAME_I_SO_RFOK => '880'; # In frame fusion gene, same gene, broken in 2 different coding introns. Same ori. RF preserved
use constant FF_55F_DIFF_SO => '860'; # 5UTR to 5UTR fusion, different genes. Same orientation
use constant FF_55F_SAME_SO => '840'; # 5UTR to 5UTR fusion, same gene. Same orientation
use constant FF_5CF_DIFF_SO_ARF => '820'; # 5UTR to coding fusion, different genes. 1 break 5UTR, 1 break coding region. Same ori. Ambiguous RF
use constant FF_5CF_SAME_SO_ARF => '800'; # 5UTR to coding fusion, same gene. 1 break in 5UTR, 1 break in coding region. Same ori. Ambiguous RF
use constant FF_F_DIFF_E_SO_ARF => '780'; # Fusion, 2 different genes, broken in coding exons. Same orientation. Ambiguous reading frame
use constant FF_F_SAME_E_SO_ARF => '760'; # Fusion, same gene, broken in coding exons. Same orientation. Ambiguous RF
use constant FF_F_DIFF_C_SO_ARF => '740'; # Fusion, 2 different genes, broken in ambiguous coding region. Same orientation. Ambiguous RF
use constant FF_F_SAME_C_SO_ARF => '720'; # Fusion, same gene, broken in ambiguous coding region. Same orientation. Ambiguous reading frame
use constant FF_TRUNCATED => '715'; # Truncated protein product. Stop codon formed at breakpoint junction.
use constant FF_OFF_DIFF_E_SO_DRF => '710'; # Out of frame fusion gene, 2 different genes, broken in coding exons. Same orientation. RF differet
use constant FF_OFF_DIFF_I_SO_DRF => '700'; # Out of frame fusion gene, 2 different genes, broken in coding introns. Same ori. RF different
use constant FF_OFF_SAME_E_SO_DRF => '690'; # Out of frame fusion gene, same gene, broken in coding exons. Same orientation. RF differet
use constant FF_OFF_SAME_I_SO_DRF => '680'; # Out of frame fusion gene, same gene, broken in 2 different coding introns. Same ori. RF different
use constant FF_53F_DIFF_SO => '660'; # 5UTR to 3UTR fusion, different genes. Same orientation
use constant FF_53F_SAME_SO => '640'; # 5UTR to 3UTR fusion, same gene. Same orientation
use constant FF_3CF_DIFF_SO_ARF => '620'; # 3UTR to coding fusion, 2 different genes. 1 break 3UTR, 1 break coding region. Same ori. Ambiguous RF
use constant FF_3CF_SAME_SO_ARF => '600'; # 3UTR to coding fusion, same gene. 1 break in 3UTR, 1 break in coding region. Same ori. Ambiguous RF
use constant FF_5CF_DIFF_OO_ARF => '580'; # 5UTR to coding fusion, 2 different genes. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
use constant FF_5CF_SAME_OO_ARF => '570'; # 5UTR to coding fusion, same gene. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
use constant FF_3CF_DIFF_OO_ARF => '560'; # 3UTR to coding fusion, 2 different genes. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
use constant FF_3CF_SAME_OO_ARF => '550'; # 3UTR to coding fusion, same gene. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
use constant FF_F_DIFF_I_OO_ARF => '540'; # fusion, 2 different genes, broken in coding introns. Opposite orientation. Ambiguous reading frame
use constant FF_F_SAME_I_OO_ARF => '530'; # fusion, same gene, broken in coding introns. Opposite orientation. Ambiguous reading frame
use constant FF_F_DIFF_E_OO_ARF => '520'; # Fusion, 2 different genes, broken in coding exons. Opposite ori. Ambiguous RF
use constant FF_F_SAME_E_OO_ARF => '510'; # Fusion, same gene, broken in coding exons. Opposite ori. Ambiguous RF
use constant FF_F_DIFF_C_OO_ARF => '500'; # Fusion, 2 different genes, broken in ambiguous coding region. Opposite ori. Ambiguous RF
use constant FF_F_SAME_C_OO_ARF => '490'; # Fusion, same gene, broken in ambiguous coding region. Opposite ori. Ambiguous RF
use constant FF_F_SAME_SI_SO => '460'; # In frame fusion gene, same gene, broken in same coding intron. Same orientation. RF unchanged
use constant FF_F_SAME_SI_OO => '450'; # Fusion, same gene, broken in same coding intron. Opposite ori. Ambiguous RF
use constant FF_35F_DIFF_SO => '440'; # 3UTR to 5UTR fusion, different genes. Same orientation
use constant FF_35F_SAME_SO => '420'; # 3UTR to 5UTR fusion, same gene. Same orientation
use constant FF_33F_DIFF_SO => '400'; # 3UTR to 3UTR fusion, different genes. Same orientation
use constant FF_33F_SAME_SO => '380'; # 3UTR to 3UTR fusion, same gene. Same orientation
use constant FF_55F_DIFF_OO => '360'; # 5UTR to 5UTR fusion, different genes. Opposite orientation
use constant FF_55F_SAME_OO => '350'; # 5UTR to 5UTR fusion, same gene. Opposite orientation
use constant FF_53F_DIFF_OO => '320'; # 5UTR to 3UTR fusion, different genes. Different orientation
use constant FF_53F_SAME_OO => '310'; # 5UTR to 3UTR fusion, same gene. Different orientation
use constant FF_35F_DIFF_OO => '280'; # 3UTR to 5UTR fusion, different genes. Different orientation
use constant FF_35F_SAME_OO => '270'; # 3UTR to 5UTR fusion, same gene. Different orientation
use constant FF_33F_DIFF_OO => '240'; # 3UTR to 3UTR fusion, different genes. Different orientation
use constant FF_33F_SAME_OO => '230'; # 3UTR to 3UTR fusion, same gene. Different orientation
use constant FF_I2OTHER      => '200'; # Intron to something else fusion
use constant FF_OTHER2I => '180'; # something else to intron fusion
use constant FF_OTHER2OTHER  => '100'; # something else to something else fusion

=head1 NAME

RGcombinationAnnotator

=head1 SYNOPSIS

use Sanger::CGP::Grass::Annotation::RGcombinationAnnotator;

my $combi = new Sanger::CGP::Grass::Annotation::RGcombinationAnnotator(-anno1   => $Anno_object1,
							  -anno2   => $Anno_object2,
							  -strand1 => $strand1,
							  -strand2 => $strand2,
							  -shard   => $shard);

$combi->combine();

my $rganno_object = $combi->anno();

=head1 DESCRIPTION

Class to combine gene annotations and predict possible rearrangement annotation (fusion) outcomes.

Requires a partially completed Anno object, generated by RGendAnnotator, for each breakpoint end.

Returns a completed Anno object.


=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#--------------------------------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::Annotation::RGcombinationAnnotator();
  Description: make a new RGcombinationAnnotator object
  Return     : RGcombinationAnnotator object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};

    bless $self,$class;

    $self->{debug} = 0;

    # set the input parameters
    if ($args{-anno1}) { $self->anno1($args{-anno1}); }
    if ($args{-anno2}) { $self->anno2($args{-anno2}); }
    if ($args{-shard}) { $self->shard($args{-shard}); }
    if ($args{-strand1}) { $self->strand1($args{-strand1}); }
    if ($args{-strand2}) { $self->strand2($args{-strand2}); }

    return $self;
}

#-----------------------------------------------------------------------#

=head2 anno1

  Arg (1)    : $anno1
  Example    : $anno1 = $Object->anno1($anno1);
  Description: get/set the anno1 object to use - the lower end of the rearrangement
  Return     : anno1 object

=cut

sub anno1 {
    my $self = shift;
    $self->{anno1} = shift if @_;
    return $self->{anno1};
}
#-----------------------------------------------------------------------#

=head2 anno2

  Arg (1)    : $anno2
  Example    : $anno2 = $Object->anno2($anno2);
  Description: get/set the anno2 object to use - the higher end of the rearrangement
  Return     : anno2 object

=cut

sub anno2 {
    my $self = shift;
    $self->{anno2} = shift if @_;
    return $self->{anno2};
}
#-----------------------------------------------------------------------#

=head2 strand1

  Arg (1)    : $strand1
  Example    : $strand1 = $object->strand1($strand1);
  Description: strand1 or rearrangement (1 or -1)
  Return     : $strand1

=cut

sub strand1 {
    my $self = shift;
    $self->{strand1} = shift if @_;
    return $self->{strand1};
}
#-----------------------------------------------------------------------#

=head2 strand2

  Arg (1)    : $strand2
  Example    : $strand2 = $object->strand2($strand2);
  Description: strand2 or rearrangement (1 or -1)
  Return     : $strand2

=cut

sub strand2 {
    my $self = shift;
    $self->{strand2} = shift if @_;
    return $self->{strand2};
}
#-----------------------------------------------------------------------#

=head2 shard

  Arg (1)    : $shard
  Example    : $shard = $object->shard($shard);
  Description: shard between breakpoints
  Return     : $shard

=cut

sub shard {
    my $self = shift;
    $self->{shard} = shift if @_;
    return $self->{shard};
}
#-----------------------------------------------------------------------#

=head2 stop

  Arg (0)    :
  Example    : $stop = $object->stop();
  Description: returns 1 if the join creates an in frame stop codon
  Return     : 1/0

=cut

sub stop {
    my $self = shift;
    $self->{stop} = shift if @_;
    return ($self->{stop} || 0);
}
#-----------------------------------------------------------------------#

=head2 anno

  Arg (1)    : $anno
  Example    : $anno = $Object->anno($anno);
  Description: completed Anno object containing both endAnno objects and the resulting fusion outcome
  Return     : anno object

=cut

sub anno {
    my $self = shift;
    $self->{anno} = shift if @_;
    return $self->{anno};
}
#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#

=head2 combine

  Arg (0)    :
  Example    : $Object->combine();
  Description: combines the annotation for 2 rearrangement ends (anno1 and anno2)
  Return     :

=cut

sub combine {
    my $self = shift;

    # each RGanno object initially contains one end of an annotation.
    # 1 is the lower end, 2 is the higher end.
    # get the appropriate fusion flag for this pair of ends
    # return the combined RGanno object

    # combine the annopoint and H/L type data into one object
    if ($self->{anno1}->id_rg && ($self->{anno1}->id_rg =~ /^d+$/)) { unless ($self->{anno1}->id_rg == $self->{anno2}->id_rg) {  die "ALERT: id_rgs dont match\n"; } }

    # do I need to swap them over if both genes are on the reverse strand and the product is currently on the plus strand?
    # need to mess around with the end_phase and phase if so

    # need to convert to using a data object that is independant of the database at some point



    # make a new annotation object with the appropriate bits from both input objects
    $self->{anno} = new Sanger::CGP::Grass::Anno(-L5    => $self->{anno1}->L5,
				    -L3    => $self->{anno1}->L3,
				    -Ltype => $self->{anno1}->Ltype,
				    -Llength => $self->{anno1}->Llength,
				    -H5    => $self->{anno2}->H5,
				    -H3    => $self->{anno2}->H3,
				    -Htype => $self->{anno2}->Htype,
				    -Hlength => $self->{anno2}->Hlength,
				    -id_rg => $self->{anno1}->id_rg );

    # work out the flag and set it, also check for in frame stops
    my $flag = $self->_getFlag();
    $self->{anno}->id_fusion_flag($flag);
}

#-------------------------------------------------------------------------------------------------------------------#
sub _getFlag {
    my $self = shift;
    my $flag = -1; # if it is left as -1, I've missed a case


    # check for a stop at the breakpoint (only for intron->intron, or exon->exon exact breakpoints)
    $flag = $self->_check4stops($flag);

    # if the gene is the same ( gene strand will always be the same but breakpoint strands may not be), and there's no stop codon formed
    if (!$self->{stop} && (($self->{anno}->{L5}->gene_id) eq ($self->{anno}->{H5}->gene_id))) {
	if ($self->{debug}) { print "NO STOP. GENE SAME, " . $self->{anno}->{L5}->gene_id . " " . $self->{anno}->{H5}->gene_id . "\n"; }
	# if the breaks both fall entirely within introns...
	    #   (the 'intron' type is only set if L5/3 are the same number and H5/3 are the same, so only have to look are H5/L5)
	if (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) &&
	    ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON)) {
	    $flag = $self->_get_same_gene_intronic_flag($flag);
	}
	# if the breaks both fall entirely within introns...
	elsif ((($self->{anno}->Ltype) eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) &&
	       (($self->{anno}->Htype) eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)) {
	    $flag = $self->_get_same_gene_exonic_flag($flag);
	}
	# coding to coding... (picks up any other combination of coding fusion)
	elsif ( (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)   ||
		 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS))     &&
	        (($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)   ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS)) ) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_F_SAME_C_SO_ARF; }
	    else                                  { $flag = FF_F_SAME_C_OO_ARF; }
	}
	# 5UTR to coding fusions...
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS)) ) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_5CF_SAME_SO_ARF; }
	    else                                  { $flag = FF_5CF_SAME_OO_ARF; }
	}
	# 3UTR to coding fusions...
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS)) ) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_3CF_SAME_SO_ARF; }
	    else                                  { $flag = FF_3CF_SAME_OO_ARF; }
	}
	# UTR fusions...
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) )) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_55F_SAME_SO; }
	    else                                  { $flag = FF_55F_SAME_OO; }
	}
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) )) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_53F_SAME_SO; }
	    else                                  { $flag = FF_53F_SAME_OO; }
	}
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) )) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_35F_SAME_SO; }
	    else                                  { $flag = FF_35F_SAME_OO; }
	}
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) )) {
 	    if ($self->strand1 eq $self->strand2) { $flag = FF_33F_SAME_SO; }
	    else                                  { $flag = FF_33F_SAME_OO; }
	}
    }


    # if the gene is different and no stop codon is formed
    elsif (!$self->{stop}) {
	if ($self->{debug}) { print "NO STOP. GENE DIFF, " . $self->{anno}->{L5}->gene_id . " " . $self->{anno}->{H5}->gene_id . "\n"; }

	# if the breaks both fall entirely within introns
	#   (the 'intron' type is only set if L5/3 are the same number and H5/3 are the same, so only have to look are H5/L5)
	if (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) &&
	    ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON)) {
	    $flag = $self->_get_diff_gene_intronic_flag($flag);
	}
	# exon to exon...
	elsif (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) &&
	       ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)) {
	    $flag = $self->_get_diff_gene_exonic_flag($flag);
	}
	# coding to coding...
	elsif ( (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)   ||
		 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS))     &&
	        (($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)   ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS)) ) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_F_DIFF_C_SO_ARF; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_F_DIFF_C_SO_ARF; }
	    else                                                            { $flag = FF_F_DIFF_C_OO_ARF; }
	}
	# 5UTR to coding fusions...
	elsif (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS)) ) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_5CF_DIFF_SO_ARF; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_5CF_DIFF_SO_ARF; }
	    else                                                            { $flag = FF_5CF_DIFF_OO_ARF; }
	}
	# 3UTR to coding fusions...
	elsif  (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) ) &&
	       ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) ||
		 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_CDS)) ) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_3CF_DIFF_SO_ARF; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_3CF_DIFF_SO_ARF; }
	    else                                                            { $flag = FF_3CF_DIFF_OO_ARF; }
	}
	# UTR fusions...
	if (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) ) &&
	    ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) )) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_55F_DIFF_SO; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_55F_DIFF_SO; }
	    else                                                            { $flag = FF_55F_DIFF_OO; }
	}
	if (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) ) &&
	    ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) )) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_53F_DIFF_SO; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_53F_DIFF_SO; }
	    else                                                            { $flag = FF_53F_DIFF_OO; }
	}
	if (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) ) &&
	    ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_5UTREXON) )) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_35F_DIFF_SO; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_35F_DIFF_SO; }
	    else                                                            { $flag = FF_35F_DIFF_OO; }
	}
	if (( ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) ) &&
	    ( ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTR)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTRINTRON)  ||
                 ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_3UTREXON) )) {
 	    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
		($self->strand1 eq $self->strand2) )                        { $flag = FF_33F_DIFF_SO; }
 	    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
		   ($self->strand1 ne $self->strand2) )                     { $flag = FF_33F_DIFF_SO; }
	    else                                                            { $flag = FF_33F_DIFF_OO; }
	}
    }

    # flag random remaining ones
#    print "test $flag\n";
    if ($flag == -1) {
	if ($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) {
	    $flag = FF_I2OTHER;
	}
	elsif ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) {
	    $flag = FF_OTHER2I;
	}
	else {
	    $flag = FF_OTHER2OTHER;
	}
    }

    if ($self->{debug}) { print "FLAG $flag\n"; }
    return ($flag);
}
#-------------------------------------------------------------------------------------------------------------------#
sub _check4stops {
    my $self = shift;
    my $flag = shift;

    # check for a stop at the breakpoint (only for intron->intron or exact exon breakpoints)
#    print "|". ($self->{anno}->Htype || '_') . "| |" .($self->{anno}->Ltype || '_') ."| |". Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON . "|\n";
    if (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON) &&
	($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_INTRON)) {
	# check for formation of in frame stop codon, change flag to 'truncated product'
	if ($self->{debug}) { print "INTRON\n"; }
	my $seq = $self->_get_join_seq('intron');
	if ($seq) {
	    $self->{stop} = $self->_contains_stop($seq);
	    if ($self->{stop}) { $flag = FF_TRUNCATED; }
	}
    }
    elsif (($self->{anno}->Ltype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON) &&
	   ($self->{anno}->Htype eq Sanger::CGP::Grass::Annotation::RGendAnnotator::REGION_TYPE_EXON)) {
         # check for formation of in frame stop codon, change flag to 'truncated product'
	if ($self->{debug}) { print "EXON\n"; }
	my $seq = $self->_get_join_seq('exon');
	if ($seq) {
	    $self->{stop} = $self->_contains_stop($seq);
	    if ($self->{stop}) { $flag = FF_TRUNCATED; }
	}
    }
    return($flag);
}
#-------------------------------------------------------------------------------------------------------------------#
# get the in frame section of sequence spanning an exon-exon join, where the breakpoints are in exons
# includes any shard
sub _get_join_seq {
    my $self = shift;
    my $intron_exon = shift;
    my $seq = '';

    my $shard = $self->shard;
    # if this is an intron intron join, dont include any shards
    if ($intron_exon eq 'intron') { $shard = ''; }
    # if the gene is on the opposite strand to read 1, need to rev comp the given shard
    elsif ( (($self->{anno}->{L5}->strand eq '1')  && ($self->{strand1} eq '-')) ||
	    (($self->{anno}->{L5}->strand eq '-1') && ($self->{strand1} eq '+')) ) { $shard = _revcomp($self->shard); }

    # need to check what strand the reads are on too...
    return unless ($self->{anno}->{L5} && $self->{anno}->{L5}->strand && $self->{anno}->{H5} && $self->{anno}->{H5}->strand && $self->{strand1} && $self->{strand2});

    # if a fusion on the + strand is formed...
    if ( (($self->{anno}->{L5}->strand eq '1') && ($self->{anno}->{H5}->strand eq '1') && ($self->{strand1} eq '+') && ($self->{strand2} eq '+')) ||
	 (($self->{anno}->{L5}->strand eq '1') && ($self->{anno}->{H5}->strand eq '-1') && ($self->{strand1} eq '+') && ($self->{strand2} eq '-')) ||
	 (($self->{anno}->{L5}->strand eq '-1') && ($self->{anno}->{H5}->strand eq '1') && ($self->{strand1} eq '-') && ($self->{strand2} eq '+')) ||
	 (($self->{anno}->{L5}->strand eq '-1') && ($self->{anno}->{H5}->strand eq '-1') && ($self->{strand1} eq '-') && ($self->{strand2} eq '-'))) {
	if ($self->{debug}) { print "HERE1 " . ($self->{anno}->{L5}->up2 || '')  . " |" . ($shard || '') . "| " . ($self->{anno}->{H5}->down2 || '') . "\n"; }
	if    (defined($self->{anno}->{L5}->phase) && ($self->{anno}->{L5}->phase eq '0')) {
	    $seq = ($shard || '') . ($self->{anno}->{H5}->down2 || '');
	}
	elsif (defined($self->{anno}->{L5}->phase) && ($self->{anno}->{L5}->phase eq '1')) {
	    my $up = '';
	    if ($self->{anno}->{L5}->up2) { ($up) = ($self->{anno}->{L5}->up2 =~ /(\D)$/); }
	    $seq = $up . ($shard || '') . ($self->{anno}->{H5}->down2 || '');
	}
	elsif (defined($self->{anno}->{L5}->phase) && ($self->{anno}->{L5}->phase eq '2')) {
	    if (($self->{anno}->{L5}->up2 || '') &&
                (length($self->{anno}->{L5}->up2)) == 2) {
		$seq = ($self->{anno}->{L5}->up2 || '') . ($shard || '') . ($self->{anno}->{H5}->down2 || '');
	    }
	}
	if ($self->{debug}) { print "seq++ is $seq\n"; }
	return($seq);
    }
    # if a fusion on minus strand is formed...
    elsif ( (($self->{anno}->{L5}->strand eq '-1') && ($self->{anno}->{H5}->strand eq '-1') && ($self->{strand1} eq '+') && ($self->{strand2} eq '+')) ||
	    (($self->{anno}->{L5}->strand eq '-1') && ($self->{anno}->{H5}->strand eq '1') && ($self->{strand1} eq '+') && ($self->{strand2} eq '-')) ||
	    (($self->{anno}->{L5}->strand eq '1') && ($self->{anno}->{H5}->strand eq '-1') && ($self->{strand1} eq '-') && ($self->{strand2} eq '+')) ||
	    (($self->{anno}->{L5}->strand eq '1') && ($self->{anno}->{H5}->strand eq '1') && ($self->{strand1} eq '-') && ($self->{strand2} eq '-'))
	    ) {
	if ($self->{debug}) { print "HERE2 " . ($self->{anno}->{L5}->up2 || '')  . " |" . ($shard || '') . "| " . ($self->{anno}->{H5}->down2 || '') . "\n"; }
	if    (defined($self->{anno}->{L5}->phase) && ($self->{anno}->{L5}->phase eq '0')) {
	    $seq = ($shard || '') . ($self->{anno}->{H5}->down2 || '');
	}
	elsif (defined($self->{anno}->{L5}->phase) && ($self->{anno}->{L5}->phase eq '1')) {
	    my $up = '';
	    if ($self->{anno}->{L5}->up2) { ($up) = ($self->{anno}->{L5}->up2 =~ /(\D)$/); }
	    $seq = $up . ($shard || '') . ($self->{anno}->{H5}->down2 || '');
	}
	elsif (defined($self->{anno}->{L5}->phase) && ($self->{anno}->{L5}->phase eq '2')) {
	    if ($self->{anno}->{L5}->up2 &&
                (length($self->{anno}->{L5}->up2)) == 2) {
		$seq = ($self->{anno}->{L5}->up2 || '') . ($shard || '') . ($self->{anno}->{H5}->down2 || '');
	    }
	}
	if ($self->{debug}) { print "seq-- is $seq\n"; }
	return($seq);
    }

    #$self->{anno}->{L5}->up2
    #$self->{anno}->{L5}->down2
    #$self->{anno}->{L5}->phase
    #$self->{anno}->{L5}->strand
    #$self->{shard}
}
#-------------------------------------------------------------------------------------------------------------------#
# this function takes a bit of sequence, scans through the triplets, returns 1 if any are a stop codon or 0 if not
sub _contains_stop {
    my $self = shift;
    my $seq = shift;
    my $stop = 0;

    return(0) unless $seq;

    my @bases = split '', uc($seq);

    while (scalar(@bases)) {
	my $b1 = shift(@bases);
	(my $b2 = shift(@bases)) if (scalar(@bases));
	(my $b3 = shift(@bases)) if (scalar(@bases));
	if ($self->{debug}) { print "code is " . ($b1 || '_') . ($b2 || '_') . ($b3 || '_') . "\n"; }
	next unless (($b1 eq 'T') && $b2 && $b3);
	if ( (($b2 eq 'A') && ($b3 eq 'G')) ||
             (($b2 eq 'A') && ($b3 eq 'A')) ||
             (($b2 eq 'G') && ($b3 eq 'A')) ) {
	    $stop = 1;
	}
    }
    return($stop);
}
#-------------------------------------------------------------------------------------------------------------------#
sub _revcomp {
    my $shard = shift;
    my $revcomp = '';
    if ($shard) { $shard = lc($shard); }
    else { return($revcomp); }

    my @array = split '', $shard;
    foreach (@array) {
	if    ($_ eq 'a') { $_ = 't'; }
	elsif ($_ eq 't') { $_ = 'a'; }
	elsif ($_ eq 'g') { $_ = 'c'; }
	elsif ($_ eq 'c') { $_ = 'g'; }
    }

    $revcomp = join '', reverse(@array);
    return($revcomp);
}
#-------------------------------------------------------------------------------------------------------------------#
# if both breakpoint fall in introns... (shards wont make any difference since intronic)
sub _get_same_gene_intronic_flag {
    my $self = shift;
    my $flag = shift;

    # introns the same...
    if ($self->{anno}->{L5}->region_number &&
	$self->{anno}->{H5}->region_number &&
	$self->{anno}->{L5}->region_number eq $self->{anno}->{H5}->region_number) {
	if ($self->strand1 eq $self->strand2) { $flag = FF_F_SAME_SI_SO; }
	else                                  { $flag = FF_F_SAME_SI_OO; }
    }
    # introns different, phase the same
    elsif ($self->{anno}->{L5}->region_number &&
	   $self->{anno}->{H5}->region_number &&
	   defined($self->{anno}->{L5}->phase) &&
	   defined($self->{anno}->{H5}->phase) &&
	   ($self->{anno}->{L5}->phase eq $self->{anno}->{H5}->phase)) {
	if ($self->strand1 eq $self->strand2) { $flag = FF_IFF_SAME_I_SO_RFOK; }
	else                                  { $flag = FF_F_SAME_I_OO_ARF; }
    }
    # introns different, phase different
    else {
	if ($self->strand1 eq $self->strand2) { $flag = FF_OFF_SAME_I_SO_DRF; }
	else                                  { $flag = FF_F_SAME_I_OO_ARF; }
    }
    return($flag);
}
#-------------------------------------------------------------------------------------------------------------------#
# if both breakpoint fall in exons... (shards will make a difference)
sub _get_same_gene_exonic_flag {
    my $self = shift;
    my $flag = shift;

    # phase the same
    if ($self->{anno}->{L5}->region_number &&
	$self->{anno}->{H5}->region_number &&
	defined($self->{anno}->{L5}->phase) &&
	defined($self->{anno}->{H5}->phase) &&
	($self->{anno}->{L5}->phase eq $self->{anno}->{H5}->phase)) {

	if ($self->strand1 eq $self->strand2) {
	    $flag = FF_IFF_SAME_E_SO_RFOK;

	    # if shard, check for phase shifting
	    # if shard, and shard is not divisible by 3, and there are no stops created, phase will be pushed out of frame
	    if ($self->{shard} && ((length($self->{shard}))%3)) { $flag = FF_OFF_SAME_E_SO_DRF; }
	}
	else { $flag = FF_F_SAME_E_OO_ARF; }
    }
    # phase different
    elsif ($self->{anno}->{L5}->region_number &&
	   $self->{anno}->{H5}->region_number &&
	   defined($self->{anno}->{L5}->phase) &&
	   defined($self->{anno}->{H5}->phase) &&
	   ($self->{anno}->{L5}->phase ne $self->{anno}->{H5}->phase)) {

	if ($self->strand1 eq $self->strand2) {
	    $flag = FF_OFF_SAME_E_SO_DRF;

	    # if shard, check for phase shifting
	    my $in_phase = 0;
	    if ($self->{shard}) { $in_phase = _check_phase_with_shard($self->{anno}->{L5}->phase, $self->{anno}->{H5}->phase, length($self->{shard})); }
	    if ($in_phase)      { $flag = FF_IFF_SAME_E_SO_RFOK; }
	}
	else { $flag = FF_F_SAME_E_OO_ARF; }
    }
    else {
	if ($self->strand1 eq $self->strand2) { $flag = FF_F_SAME_E_SO_ARF; }
	else                                  { $flag = FF_F_SAME_E_OO_ARF; }
    }
    return($flag);
}
#-------------------------------------------------------------------------------------------------------------------#
# if both breakpoint fall in introns... (shards wont make any difference since intronic)
sub _get_diff_gene_intronic_flag {
    my $self = shift;
    my $flag = shift;

    # gene strand the same, breakpoint strand same, phase the same
    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
	   ($self->{strand1} eq $self->{strand2}) &&
	defined($self->{anno}->{L5}->phase) &&
	defined($self->{anno}->{H5}->phase) &&
	($self->{anno}->{L5}->phase  eq $self->{anno}->{H5}->phase)) {
	$flag = FF_IFF_DIFF_I_SO_RFOK;
    }
    # gene strand the same, breakpoint strand same (phase different)
    elsif (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
	   ($self->{strand1} eq $self->{strand2})) {
	$flag = FF_OFF_DIFF_I_SO_DRF;
    }
    # gene strand different, breakpoint strand different (flips one of the genes), phase the same
    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
	   ($self->{strand1} ne $self->{strand2}) &&
	   defined($self->{anno}->{L5}->phase) &&
	   defined($self->{anno}->{H5}->phase) &&
	   ($self->{anno}->{L5}->phase  eq $self->{anno}->{H5}->phase)) {
	$flag = FF_IFF_DIFF_I_SO_RFOK;
    }
    # gene strand different, breakpoint strand different (flips one of the genes) (phase different)
    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
	   ($self->{strand1} ne $self->{strand2})) {
	$flag = FF_OFF_DIFF_I_SO_DRF;
    }
    # strand different, phase unknown
    else {
	$flag = FF_F_DIFF_I_OO_ARF;
    }
    return($flag);
}
#-------------------------------------------------------------------------------------------------------------------#
# if both breakpoint fall in exons... (shards will make a difference)
sub _get_diff_gene_exonic_flag {
    my $self = shift;
    my $flag = shift;

    if (($self->{anno}->{L5}->strand eq $self->{anno}->{H5}->strand) &&
	($self->{strand1} eq $self->{strand2})) {
	$flag = FF_F_DIFF_E_SO_ARF;
	# phase the same
	if ($self->{anno}->{L5}->region_number &&
	    $self->{anno}->{H5}->region_number &&
	    defined($self->{anno}->{L5}->phase) &&
	    defined($self->{anno}->{H5}->phase) &&
	    ($self->{anno}->{L5}->phase eq $self->{anno}->{H5}->phase)) {
	    $flag = FF_IFF_DIFF_E_SO_RFOK;

	    # if shard, check for phase shifting
	    if ($self->{shard} && ((length($self->{shard}))%3)) { $flag = FF_OFF_DIFF_E_SO_DRF; }
	}
	# phase different
	elsif ($self->{anno}->{L5}->region_number &&
	       $self->{anno}->{H5}->region_number &&
	       defined($self->{anno}->{L5}->phase) &&
	       defined($self->{anno}->{H5}->phase) &&
	       ($self->{anno}->{L5}->phase ne $self->{anno}->{H5}->phase)) {
	    $flag = FF_OFF_DIFF_E_SO_DRF;

	    # if shard, check for phase shifting
	    my $in_phase = 0;
	    if ($self->{shard}) { $in_phase = _check_phase_with_shard($self->{anno}->{L5}->phase, $self->{anno}->{H5}->phase, length($self->{shard})); }
	    if ($in_phase)      { $flag = FF_IFF_DIFF_E_SO_RFOK; }
	}
    }
    elsif (($self->{anno}->{L5}->strand ne $self->{anno}->{H5}->strand) &&
	   ($self->{strand1} ne $self->{strand2}) ) {
	$flag = FF_F_DIFF_E_SO_ARF;
	# phase the same
	if ($self->{anno}->{L5}->region_number &&
	    $self->{anno}->{H5}->region_number &&
	    defined($self->{anno}->{L5}->phase) &&
	    defined($self->{anno}->{H5}->phase) &&
	    ($self->{anno}->{L5}->phase eq $self->{anno}->{H5}->phase)) {
	    $flag = FF_IFF_DIFF_E_SO_RFOK;

	    # if shard, check for phase shifting
	    if ($self->{shard} && ((length($self->{shard}))%3)) { $flag = FF_OFF_DIFF_E_SO_DRF; }
	}
	# phase different
	elsif ($self->{anno}->{L5}->region_number &&
	       $self->{anno}->{H5}->region_number &&
	       defined($self->{anno}->{L5}->phase) &&
	       defined($self->{anno}->{H5}->phase) &&
	       ($self->{anno}->{L5}->phase ne $self->{anno}->{H5}->phase)) {
	    $flag = FF_OFF_DIFF_E_SO_DRF;

	    # if shard, check for phase shifting
	    my $in_phase = 0;
	    if ($self->{shard}) { $in_phase = _check_phase_with_shard($self->{anno}->{L5}->phase, $self->{anno}->{H5}->phase, length($self->{shard})); }
	    if ($in_phase)      { $flag = FF_IFF_DIFF_E_SO_RFOK; }
	}
    }
    else {
	$flag = FF_F_DIFF_E_OO_ARF;
    }
    return($flag);

}
#-------------------------------------------------------------------------------------------------------------------#
sub _check_phase_with_shard {
    my ($phase1, $phase2, $shard_length) = @_;
    my $in_frame = 0;

    # shard codon remainder
    my $shard_codon_remainder = $shard_length%3;

    if ($shard_codon_remainder) {
	my $base_count = $phase1 + (3 - $phase2) + $shard_codon_remainder;
	if ($base_count%3) {                } # base count still not divisible by 3, so codons not complete
	else               { $in_frame = 1; } # shard completes codons and pushes into frame
    }
    else {} # no phase shift

    return($in_frame);
}
#-------------------------------------------------------------------------------------------------------------------#


1;
