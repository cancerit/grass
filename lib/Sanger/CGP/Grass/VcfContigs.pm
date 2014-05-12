### Sanger::CGP::Grass::VcfContigs

#
# Author las
#
=head1 NAME

VcfConverter

=head1 SYNOPSIS

use Sanger::CGP::Grass::VcfContigs;

my $Entry = new Sanger::CGP::Grass::VcfContigs(-fai => $fai_file );

=head1 DESCRIPTION

    A class for converting contigs in an fai file to Vcf Contig objects

=head1 CONTACT

  Contact Lucy Stebbings, las

=head1 APPENDIX


=cut

package Sanger::CGP::Grass::VcfContigs;

use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Grass;
our $VERSION = Sanger::CGP::Grass->VERSION;

use strict;

1;

#-----------------------------------------------------------------------#

=head2 new

  Arg (0)    : 
  Example    : $object = new Sanger::CGP::Grass::VcfContigs();
  Description: make a set of Vcf::Contig objects from an fai file
  Return     : object

=cut

# new object
sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless $self,$class;

    if ($args{-fai})      { $self->fai($args{-fai}); }
    if ($args{-assembly}) { $self->assembly($args{-assembly}); }
    if ($args{-species})  { $self->species($args{-species}); }

    return $self;
}
#-----------------------------------------------------------------------#

=head2 fai

  Arg (1)    : fai name
  Example    : $fai = $object->fai($fai);
  Description: name of the fasta fai file to get contig info from 
  Return     : fai

=cut

sub fai {
    my $self = shift;
    $self->{fai} = shift if @_;
    return $self->{fai};
}
#-----------------------------------------------------------------------#

=head2 assembly

  Arg (1)    : assembly name
  Example    : $assembly = $object->assembly($assembly);
  Description: name of the sequence assembly
  Return     : assembly

=cut

sub assembly {
    my $self = shift;
    $self->{assembly} = shift if @_;
    return $self->{assembly};
}
#-----------------------------------------------------------------------#

=head2 species

  Arg (1)    : species name
  Example    : $species = $object->species($species);
  Description: name of the species (HUMAN,MOUSE etc)
  Return     : species

=cut

sub species {
    my $self = shift;
    $self->{species} = shift if @_;
    return $self->{species};
}
#-----------------------------------------------------------------------#

=head2 contigs

  Arg (1)    : contigs name
  Example    : $contigs = $object->contigs();
  Description: returns a set of contigs 
  Return     : reference to a sorted array of VCF::contig objects

=cut

sub contigs {
    my $self = shift;
    return $self->{contigs};
}
#-----------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

=head2 generate

  Arg (0)    : 
  Example    : $contigs = $object->generate();
  Description: generates the set of contigs from a fai file
  Return     : reference to a sorted array of VCF::contig objects

=cut

sub generate {
    my ($self) = @_;

    my $file = $self->{fai};
    my $contigs;

    open my $fh, "<$file" or die $!;
    while (my $line = <$fh>) {

	my @line = split ' ', $line;
	my $name = $line[0];
	my $contig = new Sanger::CGP::Vcf::Contig( -name => $name,
						   -length => $line[1],
						   -assembly => $self->{assembly},
						   -species => $self->{species});
	if(exists $contigs->{$name}){
	    croak ("Trying to merge contigs with conflicting data:\n".Dumper($contigs->{$name})."\n".Dumper($contig)) unless $contigs->{$name}->compare($contig);
	} else {
	    $contigs->{$name} = $contig;
	}
    }
    close $fh;

    foreach (sort _chr_sort keys %$contigs) { push @{$self->{contigs}}, $contigs->{$_}; }

    return $self->{contigs};
}
#--------------------------------------------------------------------------#
sub _chr_sort {

  # Extract the digits
  my ($number_a) = $a =~ /(\d+)/;
  my ($number_b) = $b =~ /(\d+)/;

  # Compare and return
  if ($number_a && $number_b) { return $number_a <=> $number_b }
  else                        { return $a cmp $b; }
}
#--------------------------------------------------------------------------#
#-----------------------------------------------------------------------#
1;
