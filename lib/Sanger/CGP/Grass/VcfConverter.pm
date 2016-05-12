package Sanger::CGP::Grass::VcfConverter;

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


use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

use strict;
use Sanger::CGP::Vcf::VcfUtil;

use constant SEP => "\t";
use constant NL => "\n";

1;

=head1 NAME

VcfConverter

=head1 SYNOPSIS

use Sanger::CGP::Grass::VcfConverter;

my $Entry = new Sanger::CGP::Grass::VcfConverter(-name => $name, -contigs => $contigs );

$VcfConverter->convert($wt_vcf_Sample_object, $mt_vcf_Sample_object, $ref_to_array_of_vcf_VcfProcessLog_objects, $fa_reference, $code_source_name_and_version);

=head1 DESCRIPTION

A class for converting grass bedpe output to vcf format
This should have already been processed by FlankingBases class.
Takes in an input file of coordinates in bedpe format (zero-based start coordinate, one-based end coordinate), as output by FlankingBases.
(With brassI entries, the 2nd strand will have already been flipped to brassII-like constructed orientations)

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::VcfConverter();
  Description: make a new VcfConverter object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless $self,$class;

    $self->{add_header} = [];

    if ($args{-infile})  { $self->infile($args{-infile}); }
    if ($args{-contigs}) { $self->{_contigs} = $args{-contigs}; }
    if ($args{-flip_strand}) { $self->flip_strand($args{-flip_strand}); }
    if ($args{-add_header}) { $self->add_header($args{-add_header}); }

    $self->{_format} = 'RC';

    return $self;
}
#-----------------------------------------------------------------------#

=head2 infile

  Arg (1)    : infile name
  Example    : $infile = $object->infile($infile);
  Description: name of the filtered brassI marked groups bedpe infile
  Return     : infile

=cut

sub infile {
    my $self = shift;
    $self->{infile} = shift if @_;
    return $self->{infile};
}
#-----------------------------------------------------------------------#

=head2 add_header

  Arg (1)    : array ref of strings
  Example    : ['brassVer=5.0.1']
  Description: entries to be added to VCF header
  Return     : array ref of strings

=cut

sub add_header {
    my $self = shift;
    $self->{add_header} = shift if @_;
    return $self->{add_header};
}
#-----------------------------------------------------------------------#

=head2 flip_strand

  Arg (1)    : 1/0
  Example    : $flip_strand = $object->flip_strand($flip_strand);
  Description: Whether to flip the second strand or not
  Return     : 1/0

=cut

sub flip_strand {
    my $self = shift;
    $self->{flip_strand} = shift if @_;
    return $self->{flip_strand};
}
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 convert

  Arg (0)    :
  Example    : $outfile = $object->convert();
  Description: convert the filtered bedpe infile to vcf format
  Return     : outfile

=cut

sub convert {
  my ($self,$wt_sample, $mt_sample, $process_logs, $reference_name, $input_source) = @_;

  unless ($self->{infile}) { print "No bedpe infile supplied. Exiting\n"; exit(1); }

  my $filetype = '';
  if ($self->{infile} =~ /^(\S+?)\.bedpe$/) { $self->{outfile} = $1 . '.vcf'; $filetype = 'bedpe'; }
  elsif ($self->{infile} =~ /^(\S+?)\.tab$/){ $self->{outfile} = $1 . '.vcf'; $filetype = 'tab'; }
  else                                   { $self->{outfile} = $self->{infile} . '.vcf'; }

  my $header = $self->gen_header($wt_sample, $mt_sample, $process_logs, $reference_name, $input_source);

  my $in = $self->{infile};
  my $out = $self->{outfile};
  open my $fh_in, "<$in" or die $!;
  open my $fh_out, ">$out" or die $!;

  print $fh_out $header;

  while (my $line = <$fh_in>) {
    next unless ($line);
    next if ($line =~ /^#/);
    my ($rec1, $rec2) = $self->gen_record($line, $filetype,$wt_sample, $mt_sample);
    print $fh_out $rec1.NL;
    print $fh_out $rec2.NL;
  }

  close $fh_in;
  close $fh_out;
}
#-----------------------------------------------------------------------#

sub gen_header{
  my($self,$wt_sample, $mt_sample, $process_logs, $reference_name, $input_source) = @_;

  my $contigs = $self->{_contigs};

  my $info = [
    {key => 'INFO', ID => 'SVTYPE',    Number => 1, Type => 'String',   Description => 'Type of structural variant. (All sequence is on the plus strand and in the forward direction).'},
    {key => 'INFO', ID => 'MATEID',    Number => 1, Type => 'String',   Description => 'ID of mate breakend'},
    {key => 'INFO', ID => 'HOMSEQ',    Number => 1, Type => 'String',   Description => 'Sequence of base pair identical micro-homology at event breakpoints. Plus strand sequence displayed.'},
    {key => 'INFO', ID => 'IMPRECISE', Number => 0, Type => 'Flag',     Description => 'Imprecise structural variation'},
    {key => 'INFO', ID => 'CIPOS',     Number => 2, Type => 'Integer',  Description => 'Confidence interval around POS for imprecise variants'},
    {key => 'INFO', ID => 'CIEND',     Number => 2, Type => 'Integer',  Description => 'Confidence interval around END for imprecise variants'},
    {key => 'INFO', ID => 'HOMLEN',    Number => 1, Type => 'Integer',  Description => 'Length of base pair identical micro-homology at event breakpoints'},
    {key => 'INFO', ID => 'REPS',      Number => '.', Type => 'String', Description => 'Repeat features that may contribute to rearrangement'},
    {key => 'INFO', ID => 'NPSNO',     Number => 1, Type => 'Integer',  Description => 'Number of normal panel samples with sequencing reads for this rearrangement'},
    {key => 'INFO', ID => 'NPRNO',     Number => 1, Type => 'Integer',  Description => 'Number of sequencing reads from normal panel samples for this rearrangement'},
    {key => 'INFO', ID => 'RBLAT',     Number => 1, Type => 'Integer',  Description => 'Blat score for first bkpt range versus second bkpt range'},
    {key => 'INFO', ID => 'BKDIST',    Number => 1, Type => 'Integer',  Description => 'Distance between breakpoints (-1 if difference chromosomes)'},
    {key => 'INFO', ID => 'BALS',      Number => '.', Type => 'String', Description => 'IDs of complementary rearrangements involved in balanced translocations'},
    {key => 'INFO', ID => 'INVS',      Number => '.', Type => 'String', Description => 'IDs of complementary rearrangements involved in inversion events'},
    {key => 'INFO', ID => 'FFV',       Number => 1, Type => 'Integer',  Description => 'Fusion Flag value. Best one reported.'},
    {key => 'INFO', ID => 'CNCH',      Number => 1, Type => 'String',   Description => 'Passes iterative, weighted CN boundary analysis.'},
    {key => 'INFO', ID => 'OCC',       Number => 1, Type => 'Integer',  Description => 'How many time the breakpoint appears in this dataset'},
    {key => 'INFO', ID => 'SID',       Number => 1, Type => 'String',   Description => 'Unique identifier from gene annotation source or unknown'},
    {key => 'INFO', ID => 'GENE',      Number => 1, Type => 'String',   Description => 'HUGO gene symbol or Unknown'},
    {key => 'INFO', ID => 'TID',       Number => 1, Type => 'String',   Description => 'Transcript id for transcript associated with breakpoint'},
    {key => 'INFO', ID => 'AS',        Number => 1, Type => 'String',   Description => 'Breakpoint annotation on this strand (+ or -)'},
    {key => 'INFO', ID => 'EPH',       Number => 1, Type => 'Integer',  Description => 'End phase - phase of the first coding base following the breakpoint'},
    {key => 'INFO', ID => 'PH',        Number => 1, Type => 'Integer',  Description => 'Phase at the breakpoint'},
    {key => 'INFO', ID => 'FL',        Number => 1, Type => 'String',   Description => 'If the breakpoint is at the first or last base of an exon (common in RNAseq data)'},
    {key => 'INFO', ID => 'RGN',       Number => 1, Type => 'String',   Description => 'Region where nucleotide variant (breakpoint) occurs in relation to a gene'},
    {key => 'INFO', ID => 'RGNNO',     Number => 1, Type => 'Integer',  Description => 'Number of intron/exon where nucleotide variant (breakpoint) occurs in relation to a gene'},
    {key => 'INFO', ID => 'RGNC',      Number => 1, Type => 'Integer',  Description => 'Count of total number of introns/exons in stated transcript'},
    {key => 'INFO', ID => 'TSRDS',     Number => '.', Type => 'String',  Description => 'Reads from the tumour sample ('.$mt_sample->name.') that span to this rearrangement'},
    {key => 'INFO', ID => 'NSRDS',     Number => '.', Type => 'String',  Description => 'Reads from the normal sample ('.$wt_sample->name.') that span to this rearrangement'},
    {key => 'INFO', ID => 'TRDS',      Number => '.', Type => 'String',  Description => 'Reads from the tumour sample ('.$mt_sample->name.') that contribute to this rearrangement'},
    {key => 'INFO', ID => 'NRDS',      Number => '.', Type => 'String',  Description => 'Reads from the normal sample ('.$wt_sample->name.') that contribute to this rearrangement'},
      ];

  # details info layout for the tumour and the normal column
  my $format = [
    {key => 'FORMAT', ID => 'RC', Number => 1, Type => 'Integer', Description => 'Count of countributing reads'},
      ];

  my @other;
  for my $head_line(@{$self->add_header}) {
    die "ERROR, Additional header line doesn't match expected format of '^[^=]+=.+$': $head_line\n" unless($head_line =~ m/^[^=]+=.+$/);
    my ($key, $value) = split /=/, $head_line;
    push @other, {key => $key, value => $value};
  }
  return Sanger::CGP::Vcf::VcfUtil::gen_tn_vcf_header( $wt_sample, $mt_sample, $contigs, $process_logs, $reference_name, $input_source, $info, $format, \@other);
}

#-----------------------------------------------------------------------#

sub gen_record {
  my($self, $line, $bedpe_or_tab,$wt_sample, $mt_sample) = @_;
  chomp $line;

  my $non_templated = '';
  my ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $repeats, $np_sample_count, $tumour_count, $normal_count, $np_count, $distance, $sample, $sample_type, $names, $count, $bal_trans, $inv, $occL, $occH, $copynumber_flag, $range_blat, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag, $up, $down, $microhom, $string, $samp, $readnames, $nreads, $treads);

  # each line splits into 2 records
  my $rec1 = '';
  my $rec2 = '';

  # check if its a brassI or brassII linev (for brassII start and end pos may need swapping over)
  if ($bedpe_or_tab eq 'bedpe') {
# brass I:
# chr1  start1  end1    chr2    start2  end2    id/name brass_score     strand1 strand2 repeats np_sample_count tumour_count     normal_count    np_count        bkpt_distance   sample  sample_type     names   count   bal_trans inv     occL    occH    copynumber_flag range_blat      gene    gene_id transcript_id   strand  end_phase       region  region_number   total_region_count      first/last      gene    gene_id transcript_id   strand  phase region  region_number   total_region_count      first/last      fusion_flag     Upstream_base   Downstream_base(rev_strand)
    my @entries = split "\t", $line;
    foreach (@entries) {
	    if ($_ eq '_') { $_ = ''; }
	    if ($_ eq '.') { $_ = ''; }
    }

    # get the up and down coordinate off the end of the line (FlankingBases added these)
    $down = pop @entries;
    $up = pop @entries;

    # brassI annotated (look for the tumour/normal count fields) - second strand may have already been flipped to brassII-like constructed style
    if (($entries[11] =~ /^\d+$/) &&
        ($entries[12] =~ /^\d+$/) &&
        (scalar(@entries) > 30)) {
	    ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $repeats, $np_sample_count, $tumour_count, $normal_count, $np_count, $distance, $sample, $sample_type, $names, $count, $bal_trans, $inv, $occL, $occH, $copynumber_flag, $range_blat, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag) = @entries;
	    if ($self->{flip_strand}) {
        if ($strand2 eq '+') { $strand2 = '-'; }
        elsif ($strand2 eq '-') { $strand2 = '+'; }
	    }
    }
    # brassI unannotated (look for the tumour/normal count fields) - second strand may have already been flipped to brassII-like constructed style
    elsif (($entries[11] =~ /^\d+$/) &&
           ($entries[12] =~ /^\d+$/)) {
	    ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $repeats, $np_sample_count, $tumour_count, $normal_count, $np_count, $distance, $sample, $sample_type, $names, $count, $bal_trans, $inv, $occL, $occH, $copynumber_flag, $range_blat) = @entries;
	    if ($self->{flip_strand}) {
        if ($strand2 eq '+') { $strand2 = '-'; }
        elsif ($strand2 eq '-') { $strand2 = '+'; }
	    }
    }

    # brassII annotated(look for the non-templated/microhm  fields) with readnames
    elsif ((!$entries[12] || ($entries[12] =~ /^[ATGCN]+$/i)) &&
           (!$entries[13] || ($entries[13] =~ /^[ATGCN]+$/i)) &&
           ($entries[14]) && ($entries[14] =~ /\d+:\d+:\d+/) &&
           (scalar(@entries) > 15))  {
	    ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $samp, $string, $non_templated, $microhom, $readnames, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag) = @entries;
    }
    # brassII unannotated (look for the non-templated/microhm  fields) with readnames
    elsif ((!$entries[12] || ($entries[12] =~ /^[ATGCN]+$/i)) &&
           (!$entries[13] || ($entries[13] =~ /^[ATGCN]+$/i)) &&
           ($entries[14]) && ($entries[14] =~ /\d+:\d+:\d+/))  {
	    ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $samp, $string, $non_templated, $microhom, $readnames) = @entries;
    }
    # brassII annotated(look for the non-templated/microhm  fields)
    elsif ((!$entries[12] || ($entries[12] =~ /^[ATGCN]+$/i)) &&
           (!$entries[13] || ($entries[13] =~ /^[ATGCN]+$/i)) &&
           (scalar(@entries) > 14))  {
	    ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $samp, $string, $non_templated, $microhom, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag) = @entries;
    }
    # brassII unannotated (look for the non-templated/microhm  fields)
    elsif ((!$entries[12] || ($entries[12] =~ /^[ATGCN]+$/i)) &&
           (!$entries[13] || ($entries[13] =~ /^[ATGCN]+$/i)) &&
           !($entries[14]))  {
	    ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $samp, $string, $non_templated, $microhom) = @entries;
    }

    # get the contributing samples and readnames if present
    if ($samp && $readnames) { ($nreads, $treads, $normal_count, $tumour_count) = get_readname_strings($samp, $readnames, $wt_sample, $mt_sample); }

    # add 1 to the start of each range since converting from zero based bedpe format
    $start1++;
    $start2++;
  }
  elsif ($bedpe_or_tab eq 'tab') {
# 1       +       32144920        32144921        14      -       73658034        73658033        .       T       99      5013985 PD4107a,        Chr.1  32144920(21)--T--73658034(33)  Chr.14-  (score 99) annotation stuff??? up? down?
    my @entries = split "\t", $line;
    foreach (@entries) {
	    if ($_ eq '_') { $_ = ''; }
	    if ($_ eq '.') { $_ = ''; }
    }

    # get the up and down coordinate off the end of the line (FlankingBases added these)
    $down = pop @entries;
    $up = pop @entries;

    # all tab files are brassII
    ($chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $non_templated, $microhom, $score, $name, $sample, $string, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag) = @entries;
  }

  # set the name for each breakpoint end
  my $name1 = $name . '_1';
  my $name2 = $name . '_2';

  # work out the imprecise range for each end
  my $imprecise1 = abs($end1 - $start1);
  my $imprecise2 = abs($end2 - $start2);

  # construst the forward and reverse alt fields - add 1 to all start fields since converting from zero based bedpe format
  # the flanking bases are currently fetched (in FlankingBases) on the plus and displayed on the plus strand so dont need revcomp'ing
  #
  # BRASSII issues:
  #   For each range, the lowest coordinate will be in the second field not the first if the strand is negative
  #   Microhom and non-templated sequences are only available for BrassII entries, not BrassI
  #   In the bedpe file, microhom and non-templated sequences are shown on the same strand as the firest coordinate range (ie as it is in the assembled reads)
  #     1000 genomes states that all variants are shown on the plus strand (http://www.1000genomes.org/faq/what-strand-are-variants-your-vcf-file)
  #     so convert the homseq/non-templated to the plus strand
  # the second entry for a rearrangement reads in the opposite direction so strands are flipped
  # BrassII only - microhom sequences in the bedpe file are shown on the same strand as the first coordinate range so convert to plus strand (stick with forward direction)
  #    the second instance of the microhomology sequence in the rearrangement may need reverse complementing, depending on strands

  # work out the alt field and what ref base (start or end) for record 1 and record2
  # If there is a microhomology, the second entry for the rearrangement may be the leftmost base of the microhom, not the flanking base
  my $microhom1 = '';
  my $microhom2 = '';
  my $pos1 = '';
  my $pos2 = '';
  my $alt1 = '';
  my $alt2 = '';
	#     *---------    -----------*
	#               |__|
  if (($strand1 eq '+') && ($strand2 eq '+')) {
    $microhom1 = $microhom;
    $microhom2 = $microhom;
    if ($microhom2) { $down = (substr $microhom2, 0, 1); }
    if ($end2 >= $start2) { $alt1 = $up.$non_templated.'['.$chr2.':'.$start2.'[';   $pos2 = $start2; }
    else                  { $alt1 = $up.$non_templated.'['.$chr2.':'.$end2.'[';     $pos2 = $end2; }
    if ($end1 >= $start1) { $alt2 = ']'.$chr1.':'.$start1.']'.$non_templated.$down; $pos1 = $start1; }
    else                  { $alt2 = ']'.$chr1.':'.$end1.']'.$non_templated.$down;   $pos1 = $end1; }
  }
	#     *---------    *-----------
	#               |______________|
  elsif (($strand1 eq '+') && ($strand2 eq '-')) {
    $microhom1 = $microhom;
    $microhom2 = _revcomp($microhom);
    if ($end2 >= $start2) { $alt1 = $up.$non_templated.']'.$chr2.':'.$end2.']'; $pos2 = $start2; }
    else                  { $alt1 = $up.$non_templated.']'.$chr2.':'.$start2.']'; $pos2 = $end2; }
    if ($end1 >= $start1) { $alt2 = $down.$non_templated.']'.$chr1.':'.$end1.']'; $pos1 = $start1; }
    else                  { $alt2 = $down.$non_templated.']'.$chr1.':'.$start1.']'; $pos1 = $end1; }
  }
	#     ---------*    *-----------
	#    |__________________________|
  elsif (($strand1 eq '-') && ($strand2 eq '-')) {
    $microhom1 = _revcomp($microhom);
    $microhom2 = _revcomp($microhom);
    if ($microhom1) { $up = (substr $microhom1, 0, 1); }
    if ($end2 >= $start2) { $alt1 = ']'.$chr2.':'.$start2.']'._revcomp($non_templated).$up; $pos2 = $start2; }
    else                  { $alt1 = ']'.$chr2.':'.$end2.']'._revcomp($non_templated).$up; $pos2 = $end2; }
    if ($end1 >= $start1) { $alt2 = $down._revcomp($non_templated).'['.$chr1.':'.$start1.'['; $pos1 = $start1; }
    else                  { $alt2 = $down._revcomp($non_templated).'['.$chr1.':'.$end1.'['; $pos1 = $end1; }
  }
	#     ---------*    -----------*
	#    |_____________|
  elsif (($strand1 eq '-') && ($strand2 eq '+')) {
    $microhom1 = _revcomp($microhom);
    $microhom2 = $microhom;
    if ($microhom1) { $up = (substr $microhom1, 0, 1); }
    if ($microhom2) { $down = (substr $microhom2, 0, 1); }
    if ($end2 >= $start2) { $alt1 = '['.$chr2.':'.$end2.'['._revcomp($non_templated).$up; $pos2 = $start2; }
    else                  { $alt1 = '['.$chr2.':'.$start2.'['._revcomp($non_templated).$up; $pos2 = $end2; }
    if ($end1 >= $start1) { $alt2 = '['.$chr1.':'.$end1.'['._revcomp($non_templated).$down; $pos1 = $start1; }
    else                  { $alt2 = '['.$chr1.':'.$start1.'['._revcomp($non_templated).$down; $pos1 = $end1; }
  }

#  print "$name $strand1 $strand2 $microhom $microhom1 $microhom2 $imprecise1 $imprecise2 ($end1 - $start1) ($end2 - $start2)\n";
#  exit;

  # CORE FIELDS
  # (CHR  POS     ID      REF     ALT     QUAL    FILTER)
  $rec1 .= $chr1.SEP.$pos1.SEP.$name1.SEP.$up.SEP.$alt1.SEP.$score.SEP.'.'.SEP;
  $rec2 .= $chr2.SEP.$pos2.SEP.$name2.SEP.$down.SEP.$alt2.SEP.$score.SEP.'.'.SEP;

#  if    (($strand1 eq '+') && ($end1 >= $start1)) { $rec1 .= $chr1.SEP.$start1.SEP.$name1.SEP.$up.SEP.$alt1.SEP.$score.SEP.'.'.SEP; }
#  elsif (($strand1 eq '+') && ($start1 > $end1))  { $rec1 .= $chr1.SEP.$end1.SEP.$name1.SEP.$up.SEP.$alt1.SEP.$score.SEP.'.'.SEP; }
#  elsif (($strand1 eq '-') && ($end1 >= $start1)) { $rec1 .= $chr1.SEP.$end1.SEP.$name1.SEP.$up.SEP.$alt1.SEP.$score.SEP.'.'.SEP; }
#  elsif (($strand1 eq '-') && ($start1 > $end1))  { $rec1 .= $chr1.SEP.$start1.SEP.$name1.SEP.$up.SEP.$alt1.SEP.$score.SEP.'.'.SEP; }

#  if    (($strand2 eq '+') && ($end2 >= $start2)) { $rec2 .= $chr2.SEP.$start2.SEP.$name2.SEP.$down.SEP.$alt2.SEP.$score.SEP.'.'.SEP; }
#  elsif (($strand2 eq '+') && ($start2 > $end2))  { $rec2 .= $chr2.SEP.$end2.SEP.$name2.SEP.$down.SEP.$alt2.SEP.$score.SEP.'.'.SEP; }
#  elsif (($strand2 eq '-') && ($end2 >= $start2)) { $rec2 .= $chr2.SEP.$start2.SEP.$name2.SEP.$down.SEP.$alt2.SEP.$score.SEP.'.'.SEP; }
#  elsif (($strand2 eq '-') && ($start2 > $end2))  { $rec2 .= $chr2.SEP.$end2.SEP.$name2.SEP.$down.SEP.$alt2.SEP.$score.SEP.'.'.SEP; }

  # INFO FIELDS

  # treat all as complex events
  my $svtype = 'BND';

  # split the copynumber entry into breakpoints
  $copynumber_flag ||= 0;

  # put both names in for each id for bal_trans and inv
  my @bal_trans_formatted = ();
  my $bal_trans_formatted = '';
  if ($bal_trans) {
    my @entries = split ',', $bal_trans;
    foreach my $root_name(@entries) { $bal_trans_formatted .= $root_name.'_1|'.$root_name.'_2'; }
  }
  $bal_trans_formatted = join ',', @bal_trans_formatted if (@bal_trans_formatted);

  my @inv_formatted = ();
  my $inv_formatted = '';
  if ($inv) {
    my @entries = split ',', $inv;
    foreach my $root_name(@entries) { $inv_formatted .= $root_name.'_1|'.$root_name.'_2'; }
  }
  $inv_formatted = join ',', @inv_formatted if (@inv_formatted);

  # assemble records

  # RECORD1
  $rec1 .= 'SVTYPE='.$svtype.';';
  $rec1 .= 'MATEID='.$name2.';';
  unless (($start1 == $end1) && ($start2 == $end2)) {
    $rec1 .= 'IMPRECISE;';
    if ($strand1 eq $strand2) { $rec1 .= 'CIPOS=0,'.$imprecise1.';'.'CIEND=0,'.$imprecise2.';'; }
    else                      { $rec1 .= 'CIPOS=0,'.$imprecise1.';'.'CIEND='.$imprecise2.',0;'; }
  }
  $rec1 .= 'REPS='.$repeats.';' if ($repeats);
  $rec1 .= 'NPSNO='.$np_sample_count.';' if ($np_sample_count);
  $rec1 .= 'NPRNO='.$np_count.';' if ($np_count);
  $rec1 .= 'RBLAT='.$range_blat.';' if ($range_blat);
  $rec1 .= 'BKDIST='.$distance.';' if ($distance);
  $rec1 .= 'BALS='.$bal_trans_formatted.';' if ($bal_trans_formatted);
  $rec1 .= 'INVS='.$inv_formatted.';' if ($inv_formatted);
  $rec1 .= 'FFV='.$fusion_flag.';' if (defined($fusion_flag) && $gene1 && $gene2);
  $rec1 .= 'CNCH='.$copynumber_flag.';' if ($copynumber_flag);
  $rec1 .= 'OCC='.$occL.';' if ($occL);
  if ($microhom1) {
    $rec1 .= 'HOMSEQ='.$microhom1.';';
    $rec1 .= 'HOMLEN='.length($microhom1).';';
  }
  $rec1 .= 'NRDS='.$nreads.';' if ($nreads);
  $rec1 .= 'TRDS='.$treads.';' if ($treads);
  # gene annotation
  if ($gene1) {
    $rec1 .= 'GENE='.$gene1.';';
    $rec1 .= 'SID='.$gene_id1.';';
    $rec1 .= 'TID='.$transcript_id1.';';
    $rec1 .= 'AS='.$astrand1.';';
    $rec1 .= 'EPH='.$end_phase1.';' if ($end_phase1);
    $rec1 .= 'RGN='.$region1.';' if ($region1);
    $rec1 .= 'RGNNO='.$region_number1.';' if ($region_number1);
    $rec1 .= 'RGNC='.$total_region_count1.';' if ($total_region_count1);
    $rec1 .= 'FL='.$firstlast1.';' if ($firstlast1);
  }

  # RECORD2
  # read back the other way so in effect, strands change (sequence is still on the plus strand)

  $rec2 .= 'SVTYPE='.$svtype.';';
  $rec2 .= 'MATEID='.$name1.';';
  unless (($start1 == $end1) && ($start2 == $end2)) {
    $rec2 .= 'IMPRECISE;';
    if ($strand1 eq $strand2) { $rec2 .= 'CIPOS=0,'.$imprecise2.';'.'CIEND=0,'.$imprecise1.';'; }
    else                      { $rec2 .= 'CIPOS=0,'.$imprecise2.';'.'CIEND='.$imprecise1.',0;'; }
  }
  $rec2 .= 'REPS='.$repeats.';' if ($repeats);
  $rec2 .= 'NPSNO='.$np_sample_count.';' if ($np_sample_count);
  $rec2 .= 'NPRNO='.$np_count.';' if ($np_count);
  $rec2 .= 'RBLAT='.$range_blat.';' if ($range_blat);
  $rec2 .= 'BKDIST='.$distance.';' if ($distance);
  $rec2 .= 'BALS='.$bal_trans_formatted.';' if ($bal_trans_formatted);
  $rec2 .= 'INVS='.$inv_formatted.';' if ($inv_formatted);
  $rec2 .= 'FFV='.$fusion_flag.';' if (defined($fusion_flag) && $gene1 && $gene2);
  $rec2 .= 'CNCH='.$copynumber_flag.';' if ($copynumber_flag);
  $rec2 .= 'OCC='.$occH.';' if ($occH);
  if ($microhom2) {
    $rec2 .= 'HOMSEQ='.$microhom2.';';
    $rec2 .= 'HOMLEN='.length($microhom2).';';
  }
  $rec2 .= 'NRDS='.$nreads.';' if ($nreads);
  $rec2 .= 'TRDS='.$treads.';' if ($treads);

  # gene annotation
  if ($gene2) {
    $rec2 .= 'GENE='.$gene2.';';
    $rec2 .= 'SID='.$gene_id2.';';
    $rec2 .= 'TID='.$transcript_id2.';';
    $rec2 .= 'AS='.$astrand2.';';
    $rec2 .= 'PH='.$phase2.';' if ($phase2);
    $rec2 .= 'RGN='.$region2.';' if ($region2);
    $rec2 .= 'RGNNO='.$region_number2.';' if ($region_number2);
    $rec2 .= 'RGNC='.$total_region_count2.';' if ($total_region_count2);
    $rec2 .= 'FL='.$firstlast2.';' if ($firstlast2);
  }

  if($names) {
    $rec1 .= 'TSRDS='.$names;
    $rec2 .= 'TSRDS='.$names;
  }

  # FORMAT FIELDS
  # in format put read counts for: normal, tumour

	# FORMAT
	$normal_count = q{0} unless(defined $normal_count);
	$tumour_count = q{0} unless(defined $tumour_count);
	$rec1 .= SEP.$self->{_format}.SEP.$normal_count.SEP.$tumour_count;
	$rec2 .= SEP.$self->{_format}.SEP.$normal_count.SEP.$tumour_count;

	return ($rec1,$rec2);
}
#-----------------------------------------------------------------------#
sub get_readname_strings {
  my ($samp, $readnames, $wt_sample, $mt_sample) = @_;

  my ($nreads, $treads);
  my $normal_count = 0;
  my $tumour_count = 0;
  return ($nreads, $treads) unless ($samp && $readnames);

  my @samples = split ',', $samp;
  my @tnreadnames = split /\|/, $readnames;

  foreach my $samp_name(@samples) {
    my $readname_group = (shift @tnreadnames) || '';
    my @readname_group = split ',', $readname_group;
    my $rcount = scalar(@readname_group);
    if ($samp_name eq $wt_sample->name) { $nreads = $readname_group; $normal_count = $rcount; }
    elsif ($samp_name eq $mt_sample->name) { $treads = $readname_group; $tumour_count = $rcount; }
#    print "NAME $samp_name READS $readname_group\n";
  }

  return($nreads, $treads, $normal_count, $tumour_count);
}
#-----------------------------------------------------------------------#
sub _revcomp {
  my $seq = shift;
  if ($seq) { $seq = uc($seq); }
  else { return(''); }

  my @seq = split '',$seq;
  foreach (@seq) {
    if    ($_ eq 'A') { $_ = 'T'; }
    elsif ($_ eq 'T') { $_ = 'A'; }
    elsif ($_ eq 'G') { $_ = 'C'; }
    elsif ($_ eq 'C') { $_ = 'G'; }
  }
  my $done_seq = join '', reverse(@seq);
  return($done_seq);
}
#-----------------------------------------------------------------------#
sub _rev {
  my $seq = shift;
  if ($seq) { $seq = uc($seq); }
  else { return(''); }

  my @seq = split '',$seq;
  my $done_seq = join '', reverse(@seq);
  return($done_seq);
}
#-----------------------------------------------------------------------#
1;
