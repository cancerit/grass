#!/usr/bin/perl

# adds gene and mutational consequence annotation to bedpe formatted rearrangement data

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

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';

use Grass::GenomeData;
use Grass::DataEntry;
use Grass::Annotation::RGannotator;
 
use Getopt::Long;

$| = 1;

my $within = 0; # how close to a range should the gene be before we flag it?
my $coord = 0;
my $file = '';
my $file2 = '';
my $outfile = '';
my $list_between = 0;
my $show_biotype = 0;
my $use_all_biotypes = 0;
my $gene_id_required = 0;

my $genome_cache_file = '';

my $ensembl_api = ''; # path to ensembl_api install eg /software/pubseq/PerlModules/Ensembl/www_58_1
my $species = ''; # eg HUMAN

my $help = 0;

GetOptions( 'within:s'      => \$within,
	    'coord:s'       => \$coord, # formats 1:-:123-456,2:+:345-678 or 1:-:123,2:+:345 or 1:-:123,2:+:345,atgatatat (shard at end)
	    'file:s'        => \$file,
	    'r_file:s'      => \$file2,
	    'outfile:s'     => \$outfile,
	    'genome_cache:s'=> \$genome_cache_file,
	    'ensembl_api:s' => \$ensembl_api,
	    'species:s'     => \$species,
	    'list_between'  => \$list_between,
	    'show_biotype'  => \$show_biotype,
	    'use_all_biotypes' => \$use_all_biotypes,
	    'gene_id_required' => \$gene_id_required,
	    'help'          => \$help,
         );

# check inputs
if ($help) { usage(); }

# set up access to genome data  from EnsemblDB (if species/ensembl_api supplied) or a cached flat file version (if genome_cache supplied)
my $genome_data = new Grass::GenomeData(-species      => $species,
					-ensembl_api  => $ensembl_api,
					-genome_cache => $genome_cache_file,
                                        -use_all_biotypes => $use_all_biotypes,
					-gene_id_required => $gene_id_required);

my $field = 0; # field of file that contains the coordinate string

if ($coord) { do_coord($within, $species, $coord, $list_between, $show_biotype, $genome_data); }
elsif ($file)   { $field = 0; do_file($within, $species, $file,    $outfile, $field, '', $list_between, $show_biotype, $genome_data); }
elsif ($file2)  { $field = 2; do_file($within, $species, $file2,   $outfile, $field, 'refract', $list_between, $show_biotype, $genome_data); } # 2 is the field array index that contains the refract coordinate string

#------------------------------------------------------------------------------------------------#
sub do_coord {
    my ($within, $species, $coord, $list_between, $show_biotype, $genome_data) = @_;

    my ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard) = parse_coords($coord);
    # set up the dataEntry object 
    my $entry = new Grass::DataEntry(-name       => $coord,
				     -chr1       => $chr1,
				     -strand1    => $strand1,
				     -pos1_start => $pos1_start,
				     -pos1_end   => $pos1_end,
				     -chr2       => $chr2,
				     -strand2    => $strand2,
				     -pos2_start => $pos2_start,
				     -pos2_end   => $pos2_end,
				     -shard      => $shard);
    # get fusion prediction
    my $dataset = [];
    push @$dataset, $entry;

    my $rgann = new Grass::Annotation::RGannotator(-within   => $within,
						   -dataset  => $dataset,
						   -species  => $species,
						   -list_between => $list_between,
						   -show_biotype => $show_biotype,
						   -genome_data => $genome_data);
    $rgann->getAnnotation();
    my $output_string = $rgann->format_for_printing();

    # print the results to screen
    if ($show_biotype) {
	print "gene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tfusion_flag\n";
    }
    else {
	print "gene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tfusion_flag\n";
    }
    print $output_string;

}
#------------------------------------------------------------------------------------------------#
sub do_file {
    my ($within, $species, $infile, $outfile, $field, $is_refract, $list_between, $show_biotype, $genome_data) = @_;
    unless ($outfile) { 
	$outfile = $infile;
	if ($outfile =~ /\./) { $outfile =~ s/\.(\D+)$/_ann.$1/; }
	else                  { $outfile = $infile . '_ann'; }
	if ($infile eq $outfile) { $outfile = $infile . '_ann'; }
    }

    # open the file
    open my $fh_in, "<$infile" or die $!;
    open my $fh_out, ">$outfile" or die $!;

    while (my $line = <$fh_in>) {
	chomp $line;

	# do any headers
	if ($line =~ /k?e?y?\s*EXPECTED/i) { 
	    if ($show_biotype) {
		print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tfusion_flag\n";
	    }
	    else {
		print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tfusion_flag\n";
	    }
	    next;
	}
	elsif ($line =~ /study\tsample/i) { 
	    if ($show_biotype) {
		print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tfusion_flag\n";
	    }
	    else {
		print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tfusion_flag\n";
	    }
	    next;
	}
	next if ($line =~ /^#/);

	# process main data line
	my @line = split "\t", $line;

        # decide what sort of file we have and get the coordinate details
	my $name;
	my $score;
	my ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard, $microhom);
	my ($chr1b, $strand1b, $pos1_startb, $pos1_endb, $chr2b, $strand2b, $pos2_startb, $pos2_endb, $shardb);
	my ($chr1c, $strand1c, $pos1_startc, $pos1_endc, $chr2c, $strand2c, $pos2_startc, $pos2_endc, $shardc);

	if (($line[0] =~ /^\S+$/) && ($line[1] =~ /^[+-]$/) && ($line[2] =~ /^\d+$/) && ($line[3] =~ /^\d+$/) && ($line[4] =~ /^\S+$/) && ($line[5] =~ /^[+-]$/) && ($line[6] =~ /^\d+$/) && ($line[7] =~ /^\d+$/) && ($line[8] =~ /^[\.ATGCNatgcn]+$/) && ($line[9] =~ /^[\.ATGCNatgcn]+$/)) { # brassII marked.rg format
	    ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $microhom, $shard) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9]);
	    $shard =~ s/\.//;
	    if ($shard) { $name = $chr1 . ':' . $strand1 . ':' . $pos1_start . '-' . $pos1_end . ',' . $chr2 . ':' . $strand2 . ':' . $pos2_start . '-' . $pos2_end . ',' . $shard; }
	    else        { $name = $chr1 . ':' . $strand1 . ':' . $pos1_start . '-' . $pos1_end . ',' . $chr2 . ':' . $strand2 . ':' . $pos2_start . '-' . $pos2_end; }
	}
	elsif (($line[0] =~ /^\S+$/) && ($line[1] =~ /^[+-]$/) && ($line[2] =~ /^\d+$/)&& ($line[3] =~ /^\d+$/)) { # brassI marked.rg format
	    ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7]);
	    $name = $chr1 . ':' . $strand1 . ':' . $pos1_start . '-' . $pos1_end . ',' . $chr2 . ':' . $strand2 . ':' . $pos2_start . '-' . $pos2_end;
	}
	elsif (($line[0] =~ /^\S+$/) && ($line[1] =~ /^\d+$/)&& ($line[2] =~ /^\d+$/) && ($line[3] =~ /^\S+$/) && ($line[4] =~ /^\d+$/)&& ($line[5] =~ /^\d+$/) && ($line[6] =~ /^\S+$/) && ($line[7] =~ /^\S+$/) && ($line[8] =~ /^[+-]$/) && ($line[9] =~ /^[+-]$/) ) { # brassI filter bedpe format
	    ($chr1, $pos1_start, $pos1_end, $chr2, $pos2_start, $pos2_end, $score, $name, $strand1, $strand2) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9]);
	    $name = $chr1 . ':' . $strand1 . ':' . $pos1_start . '-' . $pos1_end . ',' . $chr2 . ':' . $strand2 . ':' . $pos2_start . '-' . $pos2_end;
	}
	elsif ($is_refract) { # refract format - coordinate usually in second field, possibly more than one pair of coordinates, ?? in place of shards or unknown coordinate
	    $name = $line[$field];
	    my ($coord1,$coord2,$coord3) = split_refract_string($line[$field]);
	    if ($coord1) { 
		($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end,$shard) = parse_coords($coord1); 
	    }
	    if ($coord2) { 
		($chr1b, $strand1b, $pos1_startb, $pos1_endb, $chr2b, $strand2b, $pos2_startb, $pos2_endb,$shardb) = parse_coords($coord2);
	    }
	    if ($coord3) { 
		($chr1c, $strand1c, $pos1_startc, $pos1_endc, $chr2c, $strand2c, $pos2_startc, $pos2_endc,$shardc) = parse_coords($coord3);
	    }
	    unless ($coord1) {	    
		print $fh_out "$line\n";
		next;
	    }
	}
	elsif ($field)  { # standard format coordinate pair in specified field
	    ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end,$shard) = parse_coords($line[$field]); 
	    $name = $line[$field];
	}
	else  { # standard format coordinate pair in first field
	    ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end,$shard) = parse_coords($line[0]); 
	    $name = $line[0];
	}

	unless ($chr1 && $strand1 && $pos1_start && $pos1_end && $chr2 && $strand2 && $pos2_start && $pos2_end) {
	    print $fh_out "$line\n";
	    next;
	}

        # get the results string for each coordinate pair
	my ($out_string,$out_stringb,$out_stringc);
	$out_string = process_file_coords($line, $name, $chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard, $within, $species, $list_between, $show_biotype, $genome_data);

	if ($is_refract) { # only get these extra coodinates with refract output.
	    if ($pos1_startb) { 
		$out_stringb = process_file_coords($line, $name, $chr1b, $strand1b, $pos1_startb, $pos1_endb, $chr2b, $strand2b, $pos2_startb, $pos2_endb, $shardb, $within, $species, 0, $show_biotype, $genome_data); 
	    }
	    if ($pos1_startc) { 
		$out_stringc = process_file_coords($line, $name, $chr1c, $strand1c, $pos1_startc, $pos1_endc, $chr2c, $strand2c, $pos2_startc, $pos2_endc, $shardc, $within, $species, 0, $show_biotype, $genome_data); 
	    }

	    # merge the data from the compsite refract coordiate string. just use the first line of annotation
	    my @out_string = ();
	    my @out_stringb = ();
	    my @out_stringc = ();
	    if ($out_string =~ /^([^\n]+)\n?/)                    {  $out_string = $1;   @out_string  = split "\t", $out_string; }
	    if ($out_stringb && ($out_stringb =~ /^([^\n]+)\n?/)) {  $out_stringb = $1;  @out_stringb = split "\t", $out_stringb; }
	    if ($out_stringc && ($out_stringc =~ /^([^\n]+)\n?/)) {  $out_stringc = $1;  @out_stringc = split "\t", $out_stringc; }

	    my $count = 0;
	    foreach my $element(@out_string) {
		if ($out_stringb[$count] && ($out_stringb[$count] ne $element)) { $out_string[$count] .= '/' . $out_stringb[$count]; }
		if ($out_stringc[$count] && ($out_stringc[$count] ne $element)) { $out_string[$count] .= '/' . $out_stringc[$count]; }
		$count++;
	    }
	    $out_string = join "\t", @out_string;
	    $out_string .= "\n";
	}

	print $fh_out $out_string;
    }
}
#------------------------------------------------------------------------------------------------#
sub parse_coords {
    my $coord_string = shift;

    my ($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard);
    my ($coord1, $coord2);

    if ($coord_string =~ /;/) {
	# parse coordinate - formats 1-:123-456;2+:345-678 or 1-:123;2+:345
	($coord1, $coord2, $shard) = split ';', $coord_string;
	return(0) unless ($coord1 && $coord2);
	
	my ($chr_strand1,$pos_string1) = split ':', $coord1;
	my ($chr_strand2,$pos_string2) = split ':', $coord2;
	($chr1, $strand1) = ($chr_strand1 =~ /(\S+?)(\S)$/);
	($chr2, $strand2) = ($chr_strand2 =~ /(\S+?)(\S)$/);
	return(0) unless ($chr1 && $strand1 && $pos_string1 && $chr2 && $strand2 && $pos_string2);

	$pos_string1 =~ s/[,\s]//g; # take out commas and spaces in the coordinates
	$pos1_start = $pos_string1;
	$pos1_end = $pos_string1;
	if ($pos_string1 =~ /(\d+)-(\d+)/) { $pos1_start = $1; $pos1_end = $2; }
	
	$pos_string2 =~ s/[,\s]//g;
	$pos2_start = $pos_string2;
	$pos2_end = $pos_string2;
	if ($pos_string2 =~ /(\d+)-(\d+)/) { $pos2_start = $1; $pos2_end = $2; }

	unless ($chr1 && $strand1 && $pos1_start && $pos1_end && $chr2 && $strand2 && $pos2_start && $pos2_end) {
	    print "$coord_string not parsed correctly. skip\n";
	    return(0); 
	}
    }
    elsif ($coord_string =~ /,/) {
	# parse coordinate - formats 1:-:123-456,2:+:345-678 or 1:-:123,2:+:345
	($coord1, $coord2, $shard) = split ',', $coord_string;
	
	my ($pos_string1,$pos_string2);
	($chr1,$strand1,$pos_string1) = split ':', $coord1;
	($chr2,$strand2,$pos_string2) = split ':', $coord2;
	
	$pos_string1 =~ s/[,\s]//g; # take out commas and spaces in the coordinates
	$pos1_start = $pos_string1;
	$pos1_end = $pos_string1;
	if ($pos_string1 =~ /(\d+)-(\d+)/) { $pos1_start = $1; $pos1_end = $2; }
	
	$pos_string2 =~ s/[,\s]//g;
	$pos2_start = $pos_string2;
	$pos2_end = $pos_string2;
	if ($pos_string2 =~ /(\d+)-(\d+)/) { $pos2_start = $1; $pos2_end = $2; }
    }

    if ($shard) {
	$shard =~ s/ //g;
	unless ($shard =~ /^[aAtTgGcCnN]+$/) { print "WARN: shard sequence $shard is not valid. Can not use it\n"; $shard = undef; }
    }

    return($chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard);
}
#------------------------------------------------------------------------------------------------#
sub  split_refract_string {
    my $string = shift;

    my @coords  = ();

    my @elements = split ';', $string;


    my $last_coord = '';
    my $last_shard = '';
    foreach my $element(@elements) {
	my $current_coord1 = '';
	my $current_coord2 = '';
	my ($chr, $strand, $p1, $p2, $shard) = parse_element($element);
	if ($chr && $strand && $p1) { $current_coord1 = $chr . $strand . ':' . $p1; }
	if ($chr && $strand && $p2) { $current_coord2 = $chr . $strand . ':' . $p2; }

	if ($current_coord1 && $current_coord2 && $last_coord && $last_shard) { 
#	    print "1: $last_coord           $current_coord1           $current_coord2           $last_shard\n";
	    push @coords, ($last_coord . ';' . $current_coord1 . ';' . $last_shard);
	    $last_coord = $current_coord2; 
	}
	elsif ($current_coord1 && $current_coord2 && $last_coord) { 
#	    print "2: $last_coord           $current_coord1           $current_coord2\n";
	    push @coords, ($last_coord . ';' . $current_coord1);
	    $last_coord = $current_coord2; 
	}
	elsif ($current_coord1 && $last_coord && $last_shard) { 
#	    print "3: $last_coord           $current_coord1           $last_shard\n";
	    push @coords, ($last_coord . ';' . $current_coord1 . ';' . $last_shard);
	    $last_coord = $current_coord1; 
	}
	elsif ($current_coord1 && $last_coord) { 
#	    print "4: $last_coord           $current_coord1\n";
	    push @coords, ($last_coord . ';' . $current_coord1);
	    $last_coord = $current_coord1; 
	}
	elsif ($current_coord1 && $current_coord2) { # should not happen?
#	    print "5: $current_coord1           $current_coord2\n";
#	    push @coords, ($current_coord1 . ';' . $current_coord2);
	    $last_coord = $current_coord2; 
	}
	elsif ($current_coord1) { 
	    $last_coord = $current_coord1; 
	}
	$last_shard = $shard; 
    }

    return(@coords);
}
#------------------------------------------------------------------------------------------------#
sub parse_element {
    my $element = shift;
    my ($chr, $strand, $p1, $p2, $shard);

    if ($element =~ /^\?+$/) {
	return('', '', '', '', '');
    }
    if ($element =~ /^[A-Za-z]$/) {
	return($chr, $strand, $p1, $p2, $element);
    }
    if ($element =~ /^([^+-]+)([+-]):(\d+)-(\d+)$/) {
	$chr = $1;
	$strand = $2;
	$p1 = $3;
	$p2 = $4;
	return($chr, $strand, $p1, $p2, '');
    }
    if ($element =~ /^([^+-]+)([+-]):(\d+)$/) {
	$chr = $1;
	$strand = $2;
	$p1 = $3;
	return($chr, $strand, $p1, '', '');
    }
}
#------------------------------------------------------------------------------------------------#
sub process_file_coords {
    my ($line, $name, $chr1, $strand1, $pos1_start, $pos1_end, $chr2, $strand2, $pos2_start, $pos2_end, $shard, $within, $species, $list_between, $show_biotype, $genome_data) = @_;
    my $out_string = '';

    # set up the dataEntry object 
    my $entry = new Grass::DataEntry(-name       => $name,
				     -chr1       => $chr1,
				     -strand1    => $strand1,
				     -pos1_start => $pos1_start,
				     -pos1_end   => $pos1_end,
				     -chr2       => $chr2,
				     -strand2    => $strand2,
				     -pos2_start => $pos2_start,
				     -shard      => $shard,
				     -pos2_end   => $pos2_end);
    # get fusion prediction
    my $dataset = [];
    push @$dataset, $entry;

    my $rgann = new Grass::Annotation::RGannotator(-within   => $within,
						   -dataset  => $dataset,
						   -species  => $species,
						   -list_between => $list_between,
						   -show_biotype => $show_biotype,
						   -genome_data => $genome_data);
    $rgann->getAnnotation();
    my $string = $rgann->format_for_printing();
    my @results = split "\n", $string; 

    # get the results string (may be more than one line) and prepend original line to it
    chomp $line;
    foreach (@results) { $out_string .= $line . "\t" . $_ . "\n"; }

    return($out_string);
}
#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
sub usage {

    print <<HERE;

grass.pl

Description - gets gene annotations across rearrangement breakpoints


options...

   -within        : consider a range around the given breakpoint, within n
   -species       : species (HUMAN, MOUSE, DOG etc)
   -coord         : breakpoint coordiate pair and optional shard (formats 1:-:123-456,2:+:345-678 or 1:-:123,2:+:345 or 1:-:123,2:+:345,atgatatat)
   -r_file        : filtered refract input file - format type: tab delimited: coord string in third column. 
                    coord format As for -coord or eg 10-:92877;13+:103483915 or eg 10-:92877;13+:1034-2000;13+:103483915  (double event)   
   -file          : input file - format type: tab delimited, coord string in first column. As for -coord or eg 10-:92877;13+:103483915 (refract file format)      
   -ensembl_api   : path to Ensembl api (eg /software/pubseq/PerlModules/Ensembl/www_58_1)
   -remote        : use the remote ensembl server, do not try to use the local server
   -list_between  : list genes lying between a coordinate pair if they are one the same chromosome, on the same strand, and the distance between them is < 1MB
   -show_biotype  : shows the biotype (eg protein_coding) for each gene
   -help          : Print this message


examples...

grass.pl -coord  3:-:129389225,3:-:129390103,AA

grass.pl -r_file PD4107a.AllDisruptions.txt

grass.pl -coord 2:+:188365485-188365837,2:+:188417763-188418155 -list_between


Author : las


HERE

exit(0);
}
