[![Build Status](https://www.travis-ci.org/ba1/Vicinator.svg?branch=master)](https://www.travis-ci.org/ba1/Vicinator) 
[![codecov](https://codecov.io/gh/ba1/Vicinator/branch/master/graph/badge.svg)](https://codecov.io/gh/ba1/Vicinator) 
[![PyPI version](https://badge.fury.io/py/Vicinator.svg)](https://badge.fury.io/py/Vicinator) 
[![Requirements Status](https://requires.io/github/ba1/Vicinator/requirements.svg?branch=master)](https://requires.io/github/ba1/Vicinator/requirements/?branch=master) 
[![Documentation Status](https://readthedocs.org/projects/vicinator/badge/?version=latest)](https://vicinator.readthedocs.io/en/latest/?badge=latest) 
[![Code style:black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Vicinator

### What is Vicinator for?

Vicinator visualizes the microsynteny of grouped proteins (e.g. orthologs) across a large collection of genomes. 
As input, it requires a mapping of the genomes' proteins to the respective protein groups and a directory containing 
the genomes' feature files, i.e. files of the format *\*.gff* or *\*_feature_table.txt*.

![image](https://user-images.githubusercontent.com/8181764/104918766-86b5e980-5995-11eb-8a6b-9f2505c74973.png)


### What is Vicinator not for?

As stated above, Vicinator relies on a pre-computed grouping of proteins across genomes. It can not find these 
groups of genes for you.

### Installation

Vicinator is written for Python 3.6+

It is recommended to install Vicinator inside a virtual environment, e.g. with venv:

`python3 -m venv myenv`

This activates the new environment called *myenv*. While activated, you can install the latest version via pip. 
The following command installs the latest version and all unmet requirements automatically.

`pip install --upgrades vicinator`

Requirements:
  -    ansi2html>=1.5.2
  -    colorama>=0.4.4
  -    ete3>=3.1.2
  -    pandas>=1.1.3
  -    importlib-metadata>=3.1.1
  -    setuptools-scm>=5.0.1

### Options

```
python3 vicinator/vicinator.py --help
                                                                                                                                                                                                  
usage: Vicinator [-h] --tabular-ortholog-groups <orthology_table>
                 --feat-tables-dir <dir_path> --reference <file_path>
                 --centerprotein-accession <str> --extension-size <int>
                 [--tree <newick_tree_file_path>] [--outdir <dir_path>]
                 [--prefix <str>] [--outputlabel-map <file_path>]
                 [--nprocs <int>] [--force] [--version]

Track Microsynteny of target proteins and its orthologs across genomes.

required arguments:
  --tabular-ortholog-groups <orthology_table>
                        path to mapping file with format
                        ortholog_group_id<tab>genome_id<tab>protein_seq_id
  --feat-tables-dir <dir_path>
                        path to directory of *.feature_tables.txt or *.gff3
                        files that shall be screen

required arguments (neighborhood):
  --reference <file_path>
                        path to a ncbi style feature table file that acts as a
                        reference
  --centerprotein-accession <str>
                        unique identifier of the central gene of the window
  --extension-size <int>
                        defines the #features that are co-checked to the left
                        and right of the centerprotein

optional arguments (output):
  --tree <newick_tree_file_path>
                        path to newick tree that includes all taxa to be
                        screened
  --outdir <dir_path>   path to desired output directory
  --prefix <str>        if option is set, shows intergenic distances of genes
                        surrounding the center gene
  --outputlabel-map <file_path>
                        Attempts to replace genome accessions in the outputs
                        with a replacement string. Requires a two-column map
                        file formatted like so: 'genome file accession' <tab>
                        'replacement string'

optional arguments (run):
  --nprocs <int>        Number of CPUs for parallel processing of genomes.
                        Default: Number of CPUs-1
  --force               if option is set, existing ortholog databases in the
                        output dir are ignored and will be overwritten
```

### Input: Required Arguments

<br/>

`--tabular-ortholog-groups <orthology_table>`

>Vicinator requires a tab-separated three-column mapping of orthologs that is formatted like so:
>
> **group_id** &nbsp;&nbsp; \tab &nbsp;&nbsp;**genome_id** &nbsp;&nbsp; \tab &nbsp;&nbsp;**protein_id**
> ![example mapping file](https://user-images.githubusercontent.com/8181764/104924281-815c9d00-599d-11eb-9cb5-3e309f188bcd.png)

<br/>

`  --feat-tables-dir <dir_path>`

>Vicinator expects the path to a directory containing *.gff* format or *_feature_table.txt* 
> files of all the genomes you want to trace the microsynteny in.
>
> A recommended source for these files is NCBI RefSeq. In order for the mapping to work, the filenames 
> should correspond to the **genome_ids** specified in the mapping file:
> 
> E.g. line 7: **OG_2 &nbsp;&nbsp;  genomeB  &nbsp;&nbsp; protein_X011**
> <br/>
> triggers a search in a feature file named **genomeB.gff** or **genomeB_genomic.gff** or **genomeB_feature_table.txt** 
> in the directory specified with `--feat-tables-dir`. Effectively, it tries to locate the protein_X011 in this feature file. 

<br/>

`--reference <file_path>`
> the path to a reference genome feature file where the center-protein accession must be found

<br/>

`--centerprotein-accession` & `--extension-size <int>`

>Identifies the window of vicinity around a center-protein which is traced based on the findings in the reference 
> genome.  
> ![Vicinator Window in Reference Genome](https://user-images.githubusercontent.com/8181764/104915463-f83f6900-5990-11eb-9930-552b95109d16.png)

<br/>

## Example Basic Usage

`vicinator --tabular-ortholog-groups orthogenome_map.tsv --feat-tables-dir ./gff_dir --outdir ./results --reference gff_dir/MUSMU@10090@1.gff --centerprotein XP_006539605.1 --extension-size 3`

## Example Advanced Usage

When vicinator receives a phylogenetic tree (with genome_ids as leaf labels) it will trace the microsynteny in order of 
increasing phylogentic distance to the reference genome specified. 

`vicinator --tabular-ortholog-groups orthogenome_map.tsv --feat-tables-dir ./gff_dir --outdir ./results --reference gff_dir/MUSMU@10090@1.gff --centerprotein XP_006539605.1 --extension-size 3 --tree phylogeny.nwk`


