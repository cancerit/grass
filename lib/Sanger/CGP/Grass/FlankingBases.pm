package Sanger::CGP::Grass::FlankingBases;

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
use File::Copy qw(move);
use Bio::DB::HTS;

use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

=head1 NAME

FlankingBases

=head1 SYNOPSIS

use Sanger::CGP::Grass::FlankingBases;

my $FlankingBases = new Sanger::CGP::Grass::FlankingBases(-infile => $testfile, -ref => $ref );

# process file
$FlankingBases->process();

=head1 DESCRIPTION

Class that gets the bases upstream and downstream of the rearrangement (needed for vcf file generation).
Takes in an input file of coordinates in bedpe format (zero-based start coordinate, one-based end coordinate).
With brassI entries, the 2nd strand may have already been flipped to brassII-like constructed orientation. If not, use flip_strand option.

Takes in a reference in fasta format, and associated fai index file.

2 output fields are appended to end of line, upstream then downstream.

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    :
  Example    : $object = new Sanger::CGP::Grass::FlankingBases();
  Description: make a new object
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};


    bless $self,$class;

    # defaults
    $self->{debug} = 0;

    if ($args{-infile})  { $self->infile($args{-infile}); }
    if ($args{-ref})     { $self->ref($args{-ref}); }
    if ($args{-flip_strand}) { $self->flip_strand($args{-flip_strand}); }

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

=head2 ref

  Arg (1)    : ref
  Example    : $ref = $object->ref(ref);
  Description: reference genome in fasta format, required to retrieve L and H range sequence.
               fai index must also be present.
  Return     : ref

=cut

sub ref {
    my $self = shift;
    $self->{ref} = shift if @_;
    return $self->{ref};
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
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

=head2 process

  Arg (0)    :
  Example    : $object->process();
  Description: process the infile and put the output in the field number specified
  Return     :

=cut

sub process {
    my ($self) = @_;

    # check the file is there and valid
    my $ok = $self->_check_file($self->{infile});
    unless ($ok) { print "FlankingBases: Check failed\n"; return; }

    $ok = $self->_check_file($self->{ref}, 'ref');
    unless ($ok) { print "FlankingBases: Check ref failed\n"; return; }

    $ok = $self->_read_data();
    unless ($ok) { print "FlankingBases: Read data failed\n"; return; }

    $self->_print_file();
}
#-----------------------------------------------------------------------#
sub _check_file {
    my ($self, $file, $ref) = @_;

    unless ($file && (-e "$file")) {
	print "file $file not found\n";
	return(0);
    }

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line =~ /\S/);

	if ($ref) { $self->_check_ref_line($line); }
	else      { $self->_check_line($line); }
	last;
    }
    close $fh;

    return(1);
}
#-----------------------------------------------------------------------#
sub _check_line {
    my ($self, $line) = @_;

    chomp $line;

    my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2);

    # where there are digits in the end 2 fields, these are tumour/normal read counts and it means this is a brassI entry
    if ($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([\+-]?1?)\t([\+-]?1?)\t\S*\t\d+\t\d+/) {
	($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
	$self->{brassI} = 1;
	my $is_brass1 = 1;
	if ($self->{flip_strand}) {
	    if ($strand2 eq '+') { $strand2 = '-'; }
	    elsif ($strand2 eq '-') { $strand2 = '+'; }
	}
	return($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2,$is_brass1);
    }
    # where there are not digits in the end 2 fields, these are tumour/normal read counts and it means this is a brassII entry
    elsif ($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([\+-]?1?)\t([\+-]?1?)/) {
	($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
	$self->{brassII} = 1;
	my $is_brass1 = 0;
	return($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2,$is_brass1);
    }
    else {
	print "entry: $line\nNot in bedpe format (chr1<TAB>start1<TAB>end1<TAB>chr2<TAB>start2<TAB>end2<TAB>name<TAB>score<TAB>strand1<TAB>strand2)\n";
	return(0);
    }
}
#-----------------------------------------------------------------------#
sub _check_ref_line {
    my ($self, $line) = @_;

    chomp $line;
    return($line);
}
#-----------------------------------------------------------------------#
sub _read_data {
    my $self = shift;

    my $file = $self->{infile};

    my $data = {};

    # load an indexed fasta file
    my $fai = Bio::DB::HTS::Fai->load( $self->{ref} );

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {
	next if ($line =~ /^\s*#/);
	next unless ($line);

	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2,$is_brass1) = $self->_check_line($line);
	return(0) unless ($chr1);
	#print "HERE: $chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2\n";

  # input is bedpe format so add 1 to each start range
  $start1++;
  $start2++;

	my ($Lbase_pos, $Hbase_pos);

	# brassI coordinates - assuming coordinates are the base either side of potential rearrangement range
        #                Should I take the ref(flanking) base to the upstream side of each range instead?
	#     *---------
	#                   -----------
	if    ($strand1 eq '+') {
	    if ($start1 <= $end1) { $Lbase_pos = $start1; }
	    else                  { $Lbase_pos = $end1; }
	    if ($self->{debug}) { print "plus strand,  L base $Lbase_pos\n"; }
	}
	#     ---------*
	#                   -----------
	elsif ($strand1 eq '-') {
	    if ($start1 <= $end1) { $Lbase_pos = $end1; }
	    else                  { $Lbase_pos = $start1; }
	    if ($self->{debug}) { print "minus strand,  L base $Lbase_pos\n"; }
	}

	#     ---------
	#                   -----------*
	if    ($strand2 eq '+') {
	    if ($start2 <= $end2) { $Hbase_pos = $end2; }
	    else                  { $Hbase_pos = $start2; }
#	    if ($start2 <= $end2) { $Hbase_pos = $start2; }
#	    else                  { $Hbase_pos = $end2; }
	    if ($self->{debug}) { print "plus strand,  H base $Hbase_pos\n"; }
	}
	#     ---------
	#                   *-----------
	elsif ($strand2 eq '-') {
	    if ($start2 <= $end2) { $Hbase_pos = $start2; }
	    else                  { $Hbase_pos = $end2; }
#	    if ($start2 <= $end2) { $Hbase_pos = $end2; }
#	    else                  { $Hbase_pos = $start2; }
	    if ($self->{debug}) { print "minus strand,  H base $Hbase_pos\n"; }
	}

	my $Lstring = "$chr1:" . $Lbase_pos . "-" . $Lbase_pos;
	my $Hstring = "$chr2:" . $Hbase_pos . "-" . $Hbase_pos;
	if ($self->{debug}) { print "L $Lstring, H $Hstring\n"; }

	if ($self->{data}->{$name}->{1}) { die "duplicate entry $name in infile. Exiting.\n"; }
	if ($self->{data}->{$name}->{2}) { die "duplicate entry $name in infile. Exiting.\n"; }

	$self->{data}->{$name}->{1} = $fai->fetch("$Lstring");
	$self->{data}->{$name}->{2} = $fai->fetch("$Hstring");

	# reverse complement if on the negative strand?
        #NO, keep it all on the plus strand since VCF displays all variants on the plus strand (http://www.1000genomes.org/faq/what-strand-are-variants-your-vcf-file)
    }
    close($fh);

    return(1);
}
#-----------------------------------------------------------------------#

sub _print_file {
    my ($self) = @_;

    my $infile = $self->{infile};
    my $temp_file = $self->{infile} . '.temp';

    open my $fh, "<$infile" or die $!;
    open my $fh_temp, ">$temp_file" or die $!;

    while (my $line = <$fh>) {
	if ($line =~ /^\s*#/) { print $fh_temp $line; next; }
	next unless ($line);
	chomp $line;
	my @line = split "\t", $line;

	my $name = $line[6];

	 push @line, $self->{data}->{$name}->{1};
	 push @line, $self->{data}->{$name}->{2};

	my $done_line = join "\t", @line;
	print $fh_temp "$done_line\n";
    }
    close $fh;
    close $fh_temp;

    # check the size of the outfile is the same or greater
    my $infile_size = -s $infile;
    my $outfile_size = -s $temp_file;

    # move the new file to the old file name if the file is the expected size
    if ($outfile_size >= $infile_size) {
	move $temp_file, $infile;
    }
    else { print "WARN: FlankingBase addition failed!\n"; }
}

#-----------------------------------------------------------------------#

1;
