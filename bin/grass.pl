#!/usr/bin/perl

# adds gene and mutational consequence annotation to bedpe formatted rearrangement data

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use Sanger::CGP::Grass::GenomeData;
use Sanger::CGP::Grass::DataEntry;
use Sanger::CGP::Grass::FlankingBases;
use Sanger::CGP::Grass::Annotation::RGannotator;
use Sanger::CGP::Grass::VcfContigs;
use Sanger::CGP::Grass::VcfConverter;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Vcf::VcfProcessLog;

use Carp;
use Pod::Usage qw(pod2usage);
use Data::Dumper;
use Getopt::Long;

$| = 1;

my $within = 0; # how close to a range should the gene be before we flag it?
my $coord = 0;
my $file = '';
my $file2 = '';
my $outfile = '';
my $is_bedpe = 0;
my $list_between = 0;
my $show_biotype = 0;
my $use_all_biotypes = 0;
my $gene_id_required = 0;

# if using the genome cache
my $genome_cache_file = '';

# if using the ensembl api
my $ensembl_api = ''; # path to ensembl_api install eg /software/pubseq/PerlModules/Ensembl/www_58_1

# if using ensembl api and for vcf generation
my $species = ''; # eg HUMAN

# for generating vcf format
my $ref = '';
my $assembly = '';
my $platform = '';
my $protocol = '';
my $tumour = '';
my $acc_tumour = '';
my $acc_source_tumour = '';
my $study_tumour = '';
my $normal = '';
my $acc_normal = '';
my $acc_source_normal = '';
my $study_normal = '';

my $help = 0;
my $arg_count = scalar @ARGV;

GetOptions( 'within:s'      => \$within,
	    'coord:s'       => \$coord, # formats 1:-:123-456,2:+:345-678 or 1:-:123,2:+:345 or 1:-:123,2:+:345,atgatatat (shard at end)
	    'file:s'        => \$file,
	    'r_file:s'      => \$file2,
	    'outfile:s'     => \$outfile,
	    'genome_cache:s'=> \$genome_cache_file,
	    'ensembl_api:s' => \$ensembl_api,
	    'species:s'     => \$species,
	    'ref:s'         => \$ref,
	    'assembly:s'    => \$assembly,
	    'platform:s'    => \$platform,
	    'protocol:s'    => \$protocol,
	    'tumour:s'            => \$tumour,
	    'acc_tumour:s'        => \$acc_tumour,
	    'acc_source_tumour:s' => \$acc_source_tumour,
	    'study_tumour:s'      => \$study_tumour,
	    'normal:s'           => \$normal,
	    'acc_normal:s'       => \$acc_normal,
	    'acc_source_normal:s'=> \$acc_source_normal,
	    'study_normal:s'     => \$study_normal,
	    'list_between'  => \$list_between,
	    'show_biotype'  => \$show_biotype,
	    'use_all_biotypes' => \$use_all_biotypes,
	    'gene_id_required' => \$gene_id_required,
	    'help'          => \$help,
         );

# check inputs
if ($help || $arg_count == 0) { usage(); }

# set up access to genome data  from EnsemblDB (if species/ensembl_api supplied) or a cached flat file version (if genome_cache supplied)
my $genome_data = new Sanger::CGP::Grass::GenomeData(-species      => $species,
					-ensembl_api  => $ensembl_api,
					-genome_cache => $genome_cache_file,
                                        -use_all_biotypes => $use_all_biotypes,
					-gene_id_required => $gene_id_required);

my $field = 0; # field of file that contains the coordinate string

if ($coord) { do_coord($within, $species, $coord, $list_between, $show_biotype, $genome_data); }
elsif ($file)   { $field = 0; ($outfile, $is_bedpe) = do_file($within, $species, $file,    $outfile, $field, '', $list_between, $show_biotype, $genome_data); }
elsif ($file2)  { $field = 2; ($outfile, $is_bedpe) = do_file($within, $species, $file2,   $outfile, $field, 'refract', $list_between, $show_biotype, $genome_data); } # 2 is the field array index that contains the refract coordinate string
# note that do_file flips the orientation of the second strand of brassI entries to give brassII-like constructed strands rather than seq-read strands

unless ($is_bedpe) { print "Not bedpe format so no vcf file will be created\n";}
exit if !($is_bedpe);

# get the flanking bases for use in vcf conversion (appends the up and downstream bases to the end of each line)
if ($file || $file2) {
    my $FlankingBases = new Sanger::CGP::Grass::FlankingBases(-infile => $outfile,
							      -ref    => $ref );
    $FlankingBases->process();
}

# $outfile = 'HCC1395_191535.v1_ann.bedpe';

make_vcf_file($outfile, $ref, $species, $assembly, $platform, $protocol, $tumour, $acc_tumour, $acc_source_tumour, $study_tumour, $normal, $acc_normal, $acc_source_normal, $study_normal);

#------------------------------------------------------------------------------------------------#
sub do_coord {
    my ($within, $species, $coord, $list_between, $show_biotype, $genome_data) = @_;

    # set up the dataEntry object
    my $entry = new Sanger::CGP::Grass::DataEntry(-coord => $coord);

    # get fusion prediction
    my $dataset = [];
    push @$dataset, $entry;

    my $rgann = new Sanger::CGP::Grass::Annotation::RGannotator(-within   => $within,
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
    my $is_bedpe = 0;

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
	if (($line =~ /k?e?y?\s*EXPECTED/i) || ($line =~ /study\tsample/i) || ($line =~ /^#/)) {
	    do_header_line($fh_out, $show_biotype, $line);
	    next;
	}

	my $entry_is_bedpe = do_data_line($fh_out, $line, $field, $is_refract, $list_between, $show_biotype, $genome_data);
	if ($entry_is_bedpe) { $is_bedpe = 1; };
    }
    return($outfile, $is_bedpe);
}
#------------------------------------------------------------------------------------------------#
sub do_header_line {
    my ($fh_out, $show_biotype, $line) = @_;

    if ($line =~ /k?e?y?\s*EXPECTED/i) {
	if ($show_biotype) {
	    print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tfusion_flag\tUpstream_base\tDownstream_base(rev_strand)\n";
	}
	else {
	    print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tfusion_flag\tUpstream_base\tDownstream_base(rev_strand)\n";
	}
    }
    elsif ($line =~ /\tsample\t/i) {
	if ($show_biotype) {
	    print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tbiotype\tfusion_flag\tUpstream_base\tDownstream_base(rev_strand)\n";
	}
	else {
	    print $fh_out "$line\tgene\tgene_id\ttranscript_id\tstrand\tend_phase\tregion\tregion_number\ttotal_region_count\tfirst/last\tgene\tgene_id\ttranscript_id\tstrand\tphase\tregion\tregion_number\ttotal_region_count\tfirst/last\tfusion_flag\tUpstream_base\tDownstream_base(rev_strand)\n";
	}
    }
    elsif ($line =~ /^#/) {
	print $fh_out "$line\n";
    }
}
#------------------------------------------------------------------------------------------------#
sub do_data_line {
    my ($fh_out, $line, $field, $is_refract, $list_between, $show_biotype, $genome_data) = @_;

    return unless ($line);

    my $is_bedpe = 0;

    # process main data line
    my @line = split "\t", $line;
    my ($entry1, $entry2, $entry3);

    # decide what sort of file we have and get the coordinate details into a DataEntry object

    # brassII .tab format
    if (($line[0] =~ /^\S+$/) && ($line[1] =~ /^[+-]$/) && ($line[2] =~ /^\d+$/) && ($line[3] =~ /^\d+$/) && ($line[4] =~ /^\S+$/) && ($line[5] =~ /^[+-]$/) && ($line[6] =~ /^\d+$/) && ($line[7] =~ /^\d+$/) && ($line[8] =~ /^[\.ATGCNatgcn]+$/) && ($line[9] =~ /^[\.ATGCNatgcn]+$/)) {
	$entry1 = new Sanger::CGP::Grass::DataEntry(-chr1       => $line[0],
						    -strand1    => $line[1],
						    -pos1_start => $line[2],
						    -pos1_end   => $line[3],
						    -chr2       => $line[4],
						    -strand2    => $line[5],
						    -pos2_start => $line[6],
						    -pos2_end   => $line[7],
						    -shard      => $line[9]);
    }
    # brassI marked.rg format - need to flip the second strand orientation to give brassII type strands
    elsif (($line[0] =~ /^\S+$/) && ($line[1] =~ /^[+-]$/) && ($line[2] =~ /^\d+$/)&& ($line[3] =~ /^\d+$/)) {
	if ($line[5] eq '+') { $line[5] = '-'; }
	elsif ($line[5] eq '-') { $line[5] = '+'; }
	$entry1 = new Sanger::CGP::Grass::DataEntry(-chr1       => $line[0],
						    -strand1    => $line[1],
						    -pos1_start => $line[2],
						    -pos1_end   => $line[3],
						    -chr2       => $line[4],
						    -strand2    => $line[5],
						    -pos2_start => $line[6],
						    -pos2_end   => $line[7]);
    }
    # brassI filter bedpe format - need to flip the second strand orientation
    # also add 1 to the start of each range because bedpe is zero based and reference is 1 based.
    elsif (($line[0] =~ /^\S+$/) && ($line[1] =~ /^\d+$/)&& ($line[2] =~ /^\d+$/) && ($line[3] =~ /^\S+$/) && ($line[4] =~ /^\d+$/)&& ($line[5] =~ /^\d+$/) && ($line[6] =~ /^\S+$/) && ($line[7] =~ /^\S+$/) && ($line[8] =~ /^[+-]$/) && ($line[9] =~ /^[+-]$/)  && ($line[11] =~ /^\d+$/)  && ($line[12] =~ /^\d+$/) ) {
	if ($line[9] eq '+') { $line[9] = '-'; }
	elsif ($line[9] eq '-') { $line[9] = '+'; }
	$entry1 = new Sanger::CGP::Grass::DataEntry(-chr1       => $line[0],
						    -strand1    => $line[8],
						    -pos1_start => ($line[1] + 1),
						    -pos1_end   => ($line[2] + 1),
						    -chr2       => $line[3],
						    -strand2    => $line[9],
						    -pos2_start => $line[4],
						    -pos2_end   => $line[5]);
	$is_bedpe = 1;
    }
    # brassII filter bedpe format
    # add 1 to the start of each range because bedpe is zero based and reference is 1 based.
    elsif (($line[0] =~ /^\S+$/) && ($line[1] =~ /^\d+$/)&& ($line[2] =~ /^\d+$/) && ($line[3] =~ /^\S+$/) && ($line[4] =~ /^\d+$/)&& ($line[5] =~ /^\d+$/) && ($line[6] =~ /^\S+$/) && ($line[7] =~ /^\S+$/) && ($line[8] =~ /^[+-]$/) && ($line[9] =~ /^[+-]$/) ) {
	$entry1 = new Sanger::CGP::Grass::DataEntry(-chr1       => $line[0],
						    -strand1    => $line[8],
						    -pos1_start => ($line[1] + 1),
						    -pos1_end   => ($line[2] + 1),
						    -chr2       => $line[3],
						    -strand2    => $line[9],
						    -pos2_start => $line[4],
						    -pos2_end   => $line[5]);
	$is_bedpe = 1;
    }
    # refract format - coordinate usually in second field, possibly more than one pair of coordinates, ?? in place of shards or unknown coordinate
    elsif ($is_refract) {
	my ($coord1,$coord2,$coord3) = split_refract_string($line[$field]);
	if ($coord1) { $entry1 = new Sanger::CGP::Grass::DataEntry(-coord => $coord1); }
	if ($coord2) { $entry2 = new Sanger::CGP::Grass::DataEntry(-coord => $coord2); }
	if ($coord3) { $entry3 = new Sanger::CGP::Grass::DataEntry(-coord => $coord3); }
	unless ($coord1) {
	    print $fh_out "$line\n";
	    return;
	}
    }
    # standard format coordinate pair in specified field
    elsif ($field) { $entry1 = new Sanger::CGP::Grass::DataEntry(-coord => $line[$field]); }
    # standard format coordinate pair in first field
    else           { $entry1 = new Sanger::CGP::Grass::DataEntry(-coord => $line[0]); }

    # to make sure the 2nd strand (flipped from brassI to brassII layout) is printed, not the original line...
    $line = join "\t", @line;

    # skip it if we couldn't get an interpretable coord string
    unless ($entry1->chr1 && $entry1->strand1 && $entry1->pos1_start && $entry1->pos1_end && $entry1->chr2 && $entry1->strand2 && $entry1->pos2_start && $entry1->pos2_end) {
	print $fh_out "$line\n";
	return;
    }

    # get the annotation results string for each coordinate pair
    my ($out_string, $out_stringb, $out_stringc);
    $out_string = process_file_coords($line, $entry1, $within, $species, $list_between, $show_biotype, $genome_data);

    if ($is_refract) { # only get these extra coodinates with refract output.
	if ($entry2) { $out_stringb = process_file_coords($line, $entry2, $within, $species, 0, $show_biotype, $genome_data); }
	if ($entry3) { $out_stringc = process_file_coords($line, $entry3, $within, $species, 0, $show_biotype, $genome_data); }
	# merge the data from the compsite refract coordiate string. just use the first line of annotation
        $out_string = merge_refract($out_string, $out_stringb, $out_stringc);
    }

    print $fh_out $out_string;

    return($is_bedpe);
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
sub merge_refract {
    my ($out_string, $out_stringb, $out_stringc) = @_;

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
    return($out_string);
}
#------------------------------------------------------------------------------------------------#
sub process_file_coords {
    my ($line, $entry, $within, $species, $list_between, $show_biotype, $genome_data) = @_;
    my $out_string = '';

   # get fusion prediction
    my $dataset = [];
    push @$dataset, $entry;

    my $rgann = new Sanger::CGP::Grass::Annotation::RGannotator(-within   => $within,
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
sub make_vcf_file {
    my ($file, $ref, $species, $assembly, $platform, $protocol, $tumour, $acc_tumour, $acc_source_tumour, $study_tumour, $normal, $acc_normal, $acc_source_normal, $study_normal) = @_;

    # could get out the bam header and get the contig and sample vcf objects
    #$header = `samtools view -H /nfs/cancer_trk0003/00000009/191035.v1.brm.bam`;

    # ...or just create the contig objects from the fai file and the wt and tumour sample objects here...

    my $contigs = [];
    if ($ref && $species && $assembly) {
	my $contig_o = new Sanger::CGP::Grass::VcfContigs(-fai => ($ref . '.fai'),
							  -species => $species,
							  -assembly => $assembly);
	$contigs = $contig_o->generate();
    }

    my $mt_sample = new Sanger::CGP::Vcf::Sample( -name => $tumour,
						  -study => $study_tumour,
						  -platform => $platform,
						  -seq_protocol => $protocol,
						  -accession => $acc_tumour,
						  -accession_source => $acc_source_tumour,
						  -description => 'Mutant' );

    my $wt_sample = new Sanger::CGP::Vcf::Sample( -name => $normal,
						  -study => $study_normal,
						  -platform => $platform,
						  -seq_protocol => $protocol,
						  -accession => $acc_normal,
						  -accession_source => $acc_source_normal,
						  -description => 'Normal' );

    # make the process logs to put in the vcf header
    my $opts = join " ", $0, @ARGV;
    my @process_logs = ();
    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog( -input_vcf_source => basename($0),
							     -input_vcf_ver => Sanger::CGP::Grass->VERSION,
							     -input_vcf_param => $opts );

    my $source = basename($0). '_v'. Sanger::CGP::Grass->VERSION;

    if ($file) {
	my $VcfConverter = new Sanger::CGP::Grass::VcfConverter(-infile  => $file,
								-contigs => $contigs );
	$VcfConverter->convert($wt_sample, $mt_sample, \@process_logs, $ref, $source);
    }
}
#----------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
sub usage {

    print <<HERE;

grass.pl

Description - gets gene annotations and fusion predictions across rearrangement breakpoints


options...

   -within        : consider a range around the given breakpoint, within n
   -species       : species (HUMAN, MOUSE, DOG etc)
   -coord         : breakpoint coordiate pair and optional shard (formats 1:-:123-456,2:+:345-678 or 1:-:123,2:+:345 or 1:-:123,2:+:345,atgatatat)
   -r_file        : filtered refract input file - format type: tab delimited: coord string in third column.
                    coord format As for -coord or eg 10-:92877;13+:103483915 or eg 10-:92877;13+:1034-2000;13+:103483915  (double event)
   -file          : input file - format type: tab delimited, coord string in first column. As for -coord or eg 10-:92877;13+:103483915 (refract file format)
   -outfile       : what file to print output to. Default is input_file.bedpe
   -ensembl_api   : path to Ensembl api (eg /software/pubseq/PerlModules/Ensembl/www_58_1)
   -genome_cache  : Ensembl genome cache file - generated by Vagrent
   -list_between  : list genes lying between a coordinate pair if they are one the same chromosome, on the same strand, and the distance between them is < 1MB
   -show_biotype  : shows the biotype (eg protein_coding) for each gene
   -use_all_biotypes : 1 or 0. Use all biotypes, not just protein_coding (default).
                       Only relevant when using an Ensembl remote server since non-protein_coding entries are filtered out by Vagrent.
   -gene_id_required : get the gene_id not just its name.
                       Only relevant when using an Ensembl remote server since Vagrent currently does not supply this information.

   -ref           : fasta reference file (with associated fai file). For vcf out file generation.
   -assembly      : sequence assembly used (eg GRCh37). For vcf out file generation.
   -platform      : sequencing platform used (eg HiSeq). For vcf out file generation.
   -protocol      : sequencing experimental design (eg genomic, pulldown). For vcf out file generation.
   -tumour         : name of tumour sample. For vcf out file generation.
   -acc_tumour     : name of tumour sample accession id. For vcf out file generation.
   -acc_source_tumour : source of tumour sample accession id. For vcf out file generation.
   -study_tumour   : study id associated with tumour sample. For vcf out file generation.
   -normal         : name of normal sample. For vcf out file generation.
   -acc_normal     : name of normal sample accession id. For vcf out file generation.
   -acc_source_normal : source of normal sample accession id. For vcf out file generation.
   -study_normal   : study id associated with normal sample. For vcf out file generation.

   -help          : Print this message

examples...

grass.pl -genome_cache Homo_sapiens.GRCh37.74.vagrent.cache.gz -coord  3:-:129389225,3:-:129390103,AA

grass.pl -species HUMAN -ensembl_api www_74_1-r_file PD4107a.AllDisruptions.txt

grass.pl -genome_cache Homo_sapiens.GRCh37.74.vagrent.cache.gz -coord 2:+:188365485-188365837,2:+:188417763-188418155 -list_between

grass.pl -genome_cache Homo_sapiens.GRCh37.74.vagrent.cache.gz -species HUMAN -ref /nfs/cancer_ref01/human/37/genome.fa -assembly GRCh37 -platform HiSeq -protocol genomic -tumour HCC1395 -acc_tumour 1234 -acc_source_tumour COSMIC_SAMPLE_ID -study_tumour 111 -normal 1395BL -acc_normal 1235 -acc_source_normal COSMIC_SAMPLE_ID -study_normal 222 -file HCC1395_191535.v1.bedpe



Author : las


HERE

exit(0);
}
