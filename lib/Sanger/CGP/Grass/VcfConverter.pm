### Sanger::CGP::Grass::VcfConverter

#
# Author las
#
=head1 NAME

VcfConverter

=head1 SYNOPSIS

use Sanger::CGP::Grass::VcfConverter;

my $Entry = new Sanger::CGP::Grass::VcfConverter(-name       => $name,


 );

=head1 DESCRIPTION

    A class for converting grass bedpe output to vcf format

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Sanger::CGP::Grass::VcfConverter;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

use strict;
use Sanger::CGP::Vcf::VcfUtil;

use constant SEP => "\t";
use constant NL => "\n";

1;

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

    if ($args{-infile})  { $self->infile($args{-infile}); }
    if ($args{-contigs}) { $self->{_contigs} = $args{-contigs}; }

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
    if ($self->{infile} =~ /^(\S+?)\.bedpe$/) { $self->{outfile} = $1 . '.vcf'; $filetype = 'Ibedpe'; }
    elsif ($self->{infile} =~ /^(\S+?)\.tab$/){ $self->{outfile} = $1 . '.vcf'; $filetype = 'IItab'; }
    else                                   { $self->{outfile} = $self->{infile} . '.vcf'; }

    my $header = $self->gen_header($wt_sample, $mt_sample, $process_logs, $reference_name, $input_source);

    my $in = $self->{infile};
    my $out = $self->{outfile};
    open my $fh_in, "<$in" or die $!;
    open my $fh_out, ">$out" or die $!;

    print $fh_out $header.NL;

    while (my $line = <$fh_in>) {
	next unless ($line);
	next if ($line =~ /^#/);
	my ($rec1, $rec2) = $self->gen_record($line, $filetype);
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
	{key => 'INFO', ID => 'SVTYPE',    Number => 1, Type => 'String',   Description => 'Type of structural variant'},
	{key => 'INFO', ID => 'MATEID',    Number => 1, Type => 'String',   Description => 'ID of mate breakend'},
	{key => 'INFO', ID => 'HOMSEQ',    Number => 1, Type => 'String',   Description => 'Sequence of base pair identical micro-homology at event breakpoints'},
	{key => 'INFO', ID => 'IMPRECISE', Number => 0, Type => 'Flag',     Description => 'Imprecise structural variation'},
	{key => 'INFO', ID => 'CIPOS',     Number => 2, Type => 'Integer',  Description => 'Confidence interval around POS for imprecise variants'},
	{key => 'INFO', ID => 'CIEND',     Number => 2, Type => 'Integer',  Description => 'Confidence interval around END for imprecise variants'},
	{key => 'INFO', ID => 'HOMLEN',    Number => 1, Type => 'Integer',  Description => 'Length of base pair identical micro-homology at event breakpoints'},
	{key => 'INFO', ID => 'REPS',      Number => '.', Type => 'String',   Description => 'Repeat features that may contribute to rearrangement'},
	{key => 'INFO', ID => 'NPSNO',     Number => 1, Type => 'Integer',  Description => 'Number of normal panel samples with sequencing reads for this rearrangement'},
	{key => 'INFO', ID => 'NPRNO',     Number => 1, Type => 'Integer',  Description => 'Number of sequencing reads from normal panel samples for this rearrangement'},
	{key => 'INFO', ID => 'RBLAT',     Number => 1, Type => 'Integer',  Description => 'Blat score for first bkpt range versus second bkpt range'},
	{key => 'INFO', ID => 'BKDIST',    Number => 1, Type => 'Integer',  Description => 'Distance between breakpoints (-1 if difference chromosomes)'},
	{key => 'INFO', ID => 'BALS',      Number => '.', Type => 'String', Description => 'IDs of complementary rearrangements involved in balanced translocations'},
	{key => 'INFO', ID => 'INVS',      Number => '.', Type => 'String', Description => 'IDs of complementary rearrangements involved in inversion events'},
	{key => 'INFO', ID => 'FFV',       Number => 1, Type => 'Integer',  Description => 'Fusion Flag value. Best one reported.'},
	{key => 'INFO', ID => 'CNCH',      Number => 1, Type => 'String',   Description => 'Copynumber changepoints lie within CN_WITHIN of breakpoint (A=ASCAT,B=BATTENBERG,N=NGS_PICNIC). Only reported if both rearrangement breakpoints have adjacent changepoints'},
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
	];
	
    # details info layout for the tumor and the normal column
    my $format = [
#	{key => 'FORMAT', ID => 'GT', Number => 1, Type => 'String', Description => 'Genotype'},
	{key => 'FORMAT', ID => 'RC', Number => 1, Type => 'Integer', Description => 'Count of countributing reads'},
	];
    
    return Sanger::CGP::Vcf::VcfUtil::gen_tn_vcf_header( $wt_sample, $mt_sample, $contigs, $process_logs, $reference_name, $input_source, $info, $format, []);
}

#-----------------------------------------------------------------------#

sub gen_record {
    my($self, $line, $brassIorII) = @_;

    chomp $line;

    my $non_templated = '';
    my ($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $repeats, $np_sample_count, $tumor_count, $normal_count, $np_count, $distance, $sample, $sample_type, $names, $count, $bal_trans, $inv, $occL, $occH, $copynumber_flag, $range_blat, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag, $up, $down, $microhom);

    # each line splits into 2 records
    my $rec1 = '';
    my $rec2 = '';

    # check if its a brassI or brassII linev (for brassII start and end pos may need swapping over)
    if ($brassIorII eq 'Ibedpe') {
# brass I:
# chr1  start1  end1    chr2    start2  end2    id/name brass_score     strand1 strand2 repeats np_sample_count tumor_count     normal_count    np_count        bkpt_distance   sample  sample_type     names   count   bal_trans inv     occL    occH    copynumber_flag range_blat      gene    gene_id transcript_id   strand  end_phase       region  region_number   total_region_count      first/last      gene    gene_id transcript_id   strand  phase region  region_number   total_region_count      first/last      fusion_flag     Upstream_base   Downstream_base(rev_strand)
	my @entries = split "\t", $line;
	foreach (@entries) {
	    if ($_ eq '_') { $_ = ''; }
	    if ($_ eq '.') { $_ = ''; }
	}
	($chr1, $start1, $end1, $chr2, $start2, $end2, $name, $score, $strand1, $strand2, $repeats, $np_sample_count, $tumor_count, $normal_count, $np_count, $distance, $sample, $sample_type, $names, $count, $bal_trans, $inv, $occL, $occH, $copynumber_flag, $range_blat, $gene1, $gene_id1, $transcript_id1, $astrand1, $end_phase1, $region1, $region_number1, $total_region_count1, $firstlast1, $gene2, $gene_id2, $transcript_id2, $astrand2, $phase2, $region2, $region_number2, $total_region_count2, $firstlast2, $fusion_flag, $up, $down) = @entries;
    }

    # set the name for each breakpoint end
    my $name1 = $name . '_1';
    my $name2 = $name . '_2';

    # construst the forward and reverse alt fields
    my $alt1 = '';
    my $alt2 = '';
    if (($strand1 eq '+') && ($strand2 eq '+')) {
	$alt1 = $up.$non_templated.'['.$chr2.':'.$start2.'[';
	$alt2 = $down.$non_templated.'['.$chr1.':'.$start1.'[';
    }
    elsif (($strand1 eq '+') && ($strand2 eq '-')) {
	$alt1 = $up.$non_templated.']'.$chr2.':'.$start2.']';
	$alt2 = '['.$chr1.':'.$start1.'['.$non_templated.$down;
    }
    elsif (($strand1 eq '-') && ($strand2 eq '-')) {
	$alt1 = ']'.$chr2.':'.$start2.']'.$non_templated.$up;
	$alt2 = ']'.$chr1.':'.$start1.']'.$non_templated.$down;
    }
    elsif (($strand1 eq '-') && ($strand2 eq '+')) {
	$alt1 = '['.$chr2.':'.$start2.'['.$non_templated.$up;
	$alt2 = $down.$non_templated.']'.$chr1.':'.$start1.']';
    }

    # work out the imprecise range for each end
    my $imprecise1 = $end1 - $start1;
    my $imprecise2 = $end2 - $start2;
    my $imprecise_length = $imprecise1 + $imprecise2;

    # CORE FIELDS
    # (CHR  POS     ID      REF     ALT     QUAL    FILTER)
    $rec1 .= $chr1.SEP.$start1.SEP.$name1.SEP.$up.SEP.$alt1.SEP.$score.SEP.'.'.SEP;
    $rec2 .= $chr2.SEP.$start2.SEP.$name2.SEP.$down.SEP.$alt2.SEP.$score.SEP.'.'.SEP;

    # INFO FIELDS

    # treat all as complex events
    my $svtype = 'BND';

    # split the copynumber entry into breakpoints
    my ($cnch1, $cnch2) = split '_', $copynumber_flag;

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
    $rec1 .= 'SVTYPE='.$svtype.';';
    $rec1 .= 'MATEID='.$name2.';';
    unless (($start1 == $end1) && ($start2 == $end2)) {
	$rec1 .= 'IMPRECISE;';
	$rec1 .= 'CIPOS=0,'.$imprecise1.';';
	$rec1 .= 'CIEND=0,'.$imprecise2.';';
    }
    $rec1 .= 'REPS='.$repeats.';' if ($repeats);
    $rec1 .= 'NPSNO='.$np_sample_count.';';
    $rec1 .= 'NPRNO='.$np_count.';';
    $rec1 .= 'RBLAT='.$range_blat.';';
    $rec1 .= 'BKDIST='.$distance.';';
    $rec1 .= 'BALS='.$bal_trans_formatted.';' if ($bal_trans_formatted);
    $rec1 .= 'INVS='.$inv_formatted.';' if ($inv_formatted);
    $rec1 .= 'FFV='.$fusion_flag.';' if (defined($fusion_flag) && $gene1 && $gene2);
    $rec1 .= 'CNCH='.$cnch1.';' if ($cnch1);
    $rec1 .= 'OCC='.$occL.';';
    $rec1 .= 'HOMSEQ='.$microhom.';' if ($microhom);
    $rec1 .= 'HOMLEN='.length($microhom).';' if ($microhom);
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

    $rec2 .= 'SVTYPE='.$svtype.';';
    $rec2 .= 'MATEID='.$name1.';';
    unless (($start1 == $end1) && ($start2 == $end2)) { 
	$rec2 .= 'IMPRECISE;';
	$rec2 .= 'CIPOS=0,'.$imprecise2.';';
	$rec2 .= 'CIEND=0,'.$imprecise1.';';
    }
    $rec2 .= 'REPS='.$repeats.';' if ($repeats);
    $rec2 .= 'NPSNO='.$np_sample_count.';';
    $rec2 .= 'NPRNO='.$np_count.';';
    $rec2 .= 'RBLAT='.$range_blat.';';
    $rec2 .= 'BKDIST='.$distance.';';
    $rec2 .= 'BALS='.$bal_trans_formatted.';' if ($bal_trans_formatted);
    $rec2 .= 'INVS='.$inv_formatted.';' if ($inv_formatted);
    $rec2 .= 'FFV='.$fusion_flag.';' if (defined($fusion_flag) && $gene1 && $gene2);
    $rec2 .= 'CNCH='.$cnch2.';' if ($cnch2);
    $rec2 .= 'OCC='.$occH.';';
    $rec2 .= 'HOMSEQ='.$microhom.';' if ($microhom);
    $rec2 .= 'HOMLEN='.length($microhom).';' if ($microhom);
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

    # FORMAT FIELDS
    # in format put read counts for: normal, tumour
	
	# FORMAT
	$rec1 .= SEP.$self->{_format}.SEP.$normal_count.SEP.$tumor_count;
	$rec2 .= SEP.$self->{_format}.SEP.$normal_count.SEP.$tumor_count;
	
	return ($rec1,$rec2);
}
#-----------------------------------------------------------------------#
1;