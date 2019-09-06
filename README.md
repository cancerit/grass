GRASS
=====

Gene Rearrangement AnalySiS

Takes pairs of coordinates describing potential rearrangement events and predicts the 'most interesting' fusion consequences that could occur.

| Master | Dev |
|---|---|
|  [![Build Status](https://travis-ci.org/cancerit/grass.svg?branch=master)](https://travis-ci.org/cancerit/grass) | [![Build Status](https://travis-ci.org/cancerit/grass.svg?branch=dev)](https://travis-ci.org/cancerit/grass) |

### Dependencies/Install
This package has dependencies on other software.  These __must be installed in advance__:

* [cgpVcf v2.0+](https://github.com/cancerit/cgpVcf/releases)
* [VAGrENT v3.0+](https://github.com/cancerit/VAGrENT/releases)

And various perl modules.

Please use `setup.sh` to install other minor dependancies (perl libraries).  Setting the environment variable `CGP_PERLLIBS` allows you to to append to `PERL5LIB` during install.  Without this all dependancies are installed into the target area.

Please be aware that this expects basic C compilation libraries and tools to be available.

---


# Input

Each end coordinate may be exact or can be represented as a range.

Either an Ensembl server or an Ensembl genome cache file (generated by [VAGrENT](http://cancerit.github.io/VAGrENT/)) must be available to supply the gene information used in the annotation process.

We recommend using a cache file as this is significantly quicker.

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


## Creating a release
#### Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

#### Cutting the release
1. Update `lib/Sanger/CGP/Grass.pm` to the correct version.
2. Update `CHANGES.md` to show major items.
3. Run `./prerelease.sh`
  * Set env variable `GRASS_GRCH37_FA` to point to a human GRCh37 reference fasta+fai.
  * Set env variable `GRASS_ENS_API` to point to `.../Ensembl/www_74_1` perl API.
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.

## Docker, Singularity and Dockstore

There is a pre-built image containing this codebase on quay.io.

* [dockstore-cgpwgs][ds-cgpwgs-git]: Contains additional tools for WGS analysis.

This was primarily designed for use with dockstore.org but can be used as normal containers.

The docker images are known to work correctly after import into a singularity image.

---

LICENCE
=======
Copyright (c) 2014-2019 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

This file is part of GRASS.

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

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
