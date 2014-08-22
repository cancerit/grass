LICENCE
=======
Copyright (c) 2014 Genome Research Ltd.

Author: Lucy Stebbings <cgpit@sanger.ac.uk>

This file is part of grass.

grass is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.



GRASS
=====

Takes pairs of coordinates describing potential rearrangement events and predicts the 'most interesting' fusion consequences that could occur.

# Input


Each end coordinate may be exact or can be represented as a range.

Either an Ensembl server or an Ensembl genome cache file (generated by Vagrent) must be available to supply the gene information used in the annotation process.


# Usage

    perl grass.pl -genome_cache /<path_to>/<species>.<ref>.<ensembl_ver>.vagrent.cache.gz  -coord <coord_str>

or

    perl grass.pl -species HUMAN -ensembl_api /<path_to>/<ensembl_api> -coord <coord_str>

or

    perl grass.pl -genome_cache /<path_to>/<species>.<ref>.<ensembl_ver>.vagrent.cache.gz  -file <file_of_coords>



examples...

    grass.pl  -genome_cache Homo_sapiens.GRCh37.74.vagrent.cache.gz -coord  3:-:129389225,3:-:129390103,AA
    grass.pl -species HUMAN -ensembl_api www_74_1 -r_file PD1234a.AllDisruptions.txt
    grass.pl  -genome_cache Homo_sapiens.GRCh37.74.vagrent.cache.gz -coord 2:+:188365485-188365837,2:+:188417763-188418155 -list_between

Note that the ensembl cache file isn't an exact copy of ensembl server db since gene names are applied after downloading the cache file precursor, using the most up to date CCDS names available.


## Options

Please see the help generated by ```grass.pl``` when executed with no parameters for up-to-date information.

# Additional detail

## Fusion Flag Value Descriptions

Fusion predictions with Higher value flags are considered more 'interesting' and are chosen over lower value flags.

* 910	In frame fusion gene, 2 different genes, broken in coding exons. Same orientation. reading frame (RF) preserved
* 900	In frame fusion gene, 2 different genes, broken in coding introns. Same ori. Reading frame preserved
* 890	In frame fusion gene, same gene, broken in coding exons. Same orientation. RF preserved
* 880	In frame fusion gene, same gene, broken in 2 different coding introns. Same ori. RF preserved
* 860	5UTR to 5UTR fusion, different genes. Same orientation
* 840	5UTR to 5UTR fusion, same gene. Same orientation
* 820	5UTR to coding fusion, different genes. 1 break 5UTR, 1 break coding region. Same ori. Ambiguous RF
* 800	5UTR to coding fusion, same gene. 1 break in 5UTR, 1 break in coding region. Same ori. Ambiguous RF
* 780	Fusion, 2 different genes, broken in coding exons. Same orientation. Ambiguous reading frame
* 760	Fusion, same gene, broken in coding exons. Same orientation. Ambiguous RF
* 740	Fusion, 2 different genes, broken in ambiguous coding region. Same orientation. Ambiguous RF
* 720	Fusion, same gene, broken in ambiguous coding region. Same orientation. Ambiguous reading frame
* 715	Truncated protein product. Stop codon formed at breakpoint junction.
* 710	Out of frame fusion gene, 2 different genes, broken in coding exons. Same orientation. RF differet
* 700	Out of frame fusion gene, 2 different genes, broken in coding introns. Same ori. RF different
* 690	Out of frame fusion gene, same gene, broken in coding exons. Same orientation. RF differet
* 680	Out of frame fusion gene, same gene, broken in 2 different coding introns. Same ori. RF different
* 660	5UTR to 3UTR fusion, different genes. Same orientation
* 640	5UTR to 3UTR fusion, same gene. Same orientation
* 620	3UTR to coding fusion, 2 different genes. 1 break 3UTR, 1 break coding region. Same ori. Ambiguous RF
* 600	3UTR to coding fusion, same gene. 1 break in 3UTR, 1 break in coding region. Same ori. Ambiguous RF
* 580	5UTR to coding fusion, 2 different genes. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
* 570	5UTR to coding fusion, same gene. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
* 560	3UTR to coding fusion, 2 different genes. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
* 550	3UTR to coding fusion, same gene. 1 break 3UTR, 1 break coding region. opposite ori. Ambiguous RF
* 540	fusion, 2 different genes, broken in coding introns. Opposite orientation. Ambiguous reading frame
* 530	fusion, same gene, broken in coding introns. Opposite orientation. Ambiguous reading frame
* 520	Fusion, 2 different genes, broken in coding exons. Opposite ori. Ambiguous RF
* 510	Fusion, same gene, broken in coding exons. Opposite ori. Ambiguous RF
* 500	Fusion, 2 different genes, broken in ambiguous coding region. Opposite ori. Ambiguous RF
* 490	Fusion, same gene, broken in ambiguous coding region. Opposite ori. Ambiguous RF
* 460	In frame fusion gene, same gene, broken in same coding intron. Same orientation. RF unchanged
* 450	Fusion, same gene, broken in same coding intron. Opposite ori. Ambiguous RF
* 440	3UTR to 5UTR fusion, different genes. Same orientation
* 420	3UTR to 5UTR fusion, same gene. Same orientation
* 400	3UTR to 3UTR fusion, different genes. Same orientation
* 380	3UTR to 3UTR fusion, same gene. Same orientation
* 360	5UTR to 5UTR fusion, different genes. Opposite orientation
* 350	5UTR to 5UTR fusion, same gene. Opposite orientation
* 320	5UTR to 3UTR fusion, different genes. Different orientation
* 310	5UTR to 3UTR fusion, same gene. Different orientation
* 280	3UTR to 5UTR fusion, different genes. Different orientation
* 270	3UTR to 5UTR fusion, same gene. Different orientation
* 240	3UTR to 3UTR fusion, different genes. Different orientation
* 230	3UTR to 3UTR fusion, same gene. Different orientation
* 200	Intron to something else fusion
* 180	something else to intron fusion
* 100	something else to something else fusion


## Processing logic

Start with breakpoint1, shard (optional), breakpoint2

Each breakpoint may be a defined base or a range.

Get Ensembl genomic 'slice' for each breakpoint (+- a distance specified by the 'within').
The start of the slice is numbered from 1 -> n.

Get genes/transcripts spanning each breakpoint -
get ccds transcripts if there are any.
Pick the longest coding one.
The coordinates for the genes/transcripts are relative to the slice so may go negative (can be confusing).

See if each breakpoint overlaps, is entirely within, or spans completely
-the coding region of the transcript
-the transcribed region of the gene

## resulting region types:

* exon - entirely in a single exon
* intron - entirely in a single intron
* CDS - entirely in coding region
* maybeCDS - overlapping coding region
* 5UTR - entirely in the 5' untranslated region
* 5UTRintron - entirely in the 5' untranslated region, in a single intron
* 5UTRexon - entirely in the 5' untranslated region, in a single exon
* 3UTR - entirely in the 3' untranslated region
* 3UTRintron - entirely in the 3' untranslated region, in a single intron
* 3UTRexon - entirely in the 3' untranslated region, in a single exon
* UTR - entirely in the untranslated region of non-coding gene
* UTRintron - entirely in the untranslated region of non-coding gene, in a single intron
* UTRexon - entirely in the untranslated region of non-coding gene, in a single exon
* nonCDS - overlapping the untranslated region
* 5nonCDS - overlapping the 5' untranslated region
* 3nonCDS - overlapping the 3' untranslated region
* upstream - entirely upstream of transcribed region
* downstream - entirely downstream of transcribed region
* For breakpoints in introns and exons, work out which number intron or exon it is in.

For exact exonic breakpoints, see if this is the very first or last base of an exon (common for RNAseq data).

## Phases

If a breakpoint is in an intron, get the phase (start phase) for the following exon.

* Phase/end phase is defined between bases.
* end phase = phase just after a base
* phase (start) = phase just before a base
* phase 0 means this base is the first of the codon and no preceding bases are needed to complete the codon.
* phase 1 means the preceding 1 base is needed to complete the codon.
* phase 2 means the preceding 2 bases are needed to complete the codon.
* If breakpoint 1 is exact, and in an exon, get the end phase
* If breakpoint 2 is exact, and in an exon, get the phase
* If breakpoint 1 end phase = breakpoint 2 start phase, an in frame product may be formed.

## Combining ends

If both breakpoint ends are in coding regions and the gene orientations are the same, a fusion product may form.
Whether an in frame fusion forms depends on the phase at each, if there is a shard, and whether a stop codon is created.

Regulatory gene region to different gene fusions are also very important.

Test each combination of transcript and check the fusion flag value.
Filter results and report only the highest 2 fusion flag value combinations.

## Strands

Must look for potential fusion products on both + and - strands.

Take into account the strand the breakpoint is on and the strand the gene is on to work out if both disrupted genes are in the same orientation.
There is a grid of possible outcomes.

## Stop codons

A stop codon may be formed when 2 coding regions fuse (only checking same orientation fusions for stops).
If the breakpoints are exact, get the 2 bases upstream and downstream of the breakpoint, and get any shard.
Check the phase at the upstream end and work through the codons across the join to see if any are stops (TGA,TAA,TAG).
A stop codon results in a fusion flag of 715.

# Comparing to older grass versions

To get the equivalent of old Grass's output (tied to Sanger architecture) run with the -use_all_biotypes and the -gene_id_required flags and supply ensembl connection details.

Comparing old and new grass. Example commands...

old Grass (held locally at the Sanger) command:

```
/.../CGP/projects/Grass/perl/scripts/grass.pl -remote -ensembl_api 74 -coord X:-:84500101-84500632,X:+:84510265-84510265
```

PanCan grass equivalent command (for comparison with old Grass, this also uses remote ensembl server):

```
perl -I /.../CGP/projects/vagrent_vcf/perl/modules/ ~/CGP_git/grass/bin/grass.pl -species HUMAN -ensembl_api /software/pubseq/PerlModules/Ensembl/www_74_1 -gene_id_required -use_all_biotypes -coord X:-:84500101-84500632,X:+:84510265-84510265
```

Comparing new grass using ensembl server, and new grass using Vagrent generated ensembl cache file. Example commands...

PanCan grass command - using remote ensembl server (for comparison with cache ensembl file):

```
perl -I /software/CGP/projects/vagrent_vcf/perl/modules/ ~/CGP_git/grass/bin/grass.pl -species HUMAN -ensembl_api /software/pubseq/PerlModules/Ensembl/www_74_1 -coord X:-:84500101-84500632,X:+:84510265-84510265
```

PanCan grass command - using cache ensembl file:

```
perl -I /software/CGP/projects/vagrent_vcf/perl/modules/ ~/CGP_git/grass/bin/grass.pl  -genome_cache /lustre/scratch104/sanger/am3/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz  -coord X:-:84500101-84500632,X:+:84510265-84510265
```

Note that the ensembl cache file isn't an exact copy of ensembl server db since gene names were applied after downloading the cache file precursor, using the most up to date CCDS names available.
