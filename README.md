# TMarSel

**TMarSel** is a tool for a Tailored Selection of gene families as Markers for microbial phylogenomics.

## Table of contents

- [Inputs](#inputs)
- [Outputs](#outputs)
- [Installation](#installation)
- [Basic usage](#basic-usage)
- [Examples](#examples)
- [Citation](#citation)

## Flags

* `-h` Display help message.
* `-i` [**required**] Either a single annotation of ORFs into gene families file OR a direcotry containing multiple annotation files. File(s) contains three columns: `orf|bit_score|gene_family`
* `-o` [**required**] Output directory to save the ORFs and statistics of each marker.
* `-raw` [**required**] IF input file(s) contain raw annotations from functional databases, which contains multiple columns depending on the database.
* `-db` [**required**] IF input contains the raw annotations set in `-raw`. Name of the database used for genome annotation (`eggnog` or `kegg`).
* `-k` Number of markers to select (default is 50).
* `-min_markers` Minimum number of markers per genome. Can be a percentage or a number (default is 1). Genomes with fewer markers than the indicated value are discarded.
* `-th` Threshold for filtering copies of each gene family per genome (default is 1.0) Retain the ORFs within `-th` of the maximum bit score for each gene family and genome. Lower values (e.g. 0.0) retains all ORFs, whereas higher values (e.g. 1.0) retains only the ORF with the highest bit score.
* `-p` Exponent of the power mean (cost function) used to select markers (default is 0.0). We recommend not changing this value unless you are familiar with the method. Default value yields the optimal combination of markers.

## Outputs

1. `k` files containing the ORFs, genome and file of origin for each marker (see below). Files are saved to `./output_dir/orfs`. 

| orf | genome | file |
| --- | --- | --- |
| G000006605_1748 | G000006605 | kofamscan_wol2_example.tsv |
| G000006725_378 | G000006725 | kofamscan_wol2_example.tsv |
| ... | ... | ... |

2. Statistics. Files are saved to `./output_dir/statistics`

* Number of markers per genome (see below). A given marker can contain more than one ORF per genome, therefore we provide the number of different markers (`k`) and the total number of markers. The `details` column is `;` separated. Each item indicates the marker name and the number of ORFs (i.e. copies) in the genome.

| genome | number_of_different_markers | total_number of markers | details |
| --- | --- | --- | --- |
| G000006605 | 10 | 10 | K01889:1;K01866:1;K01872:1; ... |
| G000006725 | 9 | 10 | K02358:2;K01872:1;K01866:1; ... |
| ... | ... | ... | ... |

* Number of genomes per marker (see below). We provide the number of genomes containing the marker. The `details` column contains the genome name and the number of ORFs of the marker.

| marker | number_of_genomes | details |
| --- | --- | --- |
| K01409 | 1509 | G000093065:2;G900097235:2;G002074035:2; ... |
| K01866 | 1508 | G900097235:2;G001941465:2;G000006605:1; ... |

## Installation

* **Anaconda**
* **pip**
* **GitHub**

## Basic usage

```bash
 python TMarSel.py -i input_file -o output_dir
```

After installation, type `python TMarSel.py -h` to learn all the options.

## Examples

We provide multiple examples to showcase the usage of **TMarSel**. Data can be downloaded as explained in [files](data/files.md).

### 1\. Annotations of 1,510 genomes from the Web of Life 2 database

* EggNOG annotations contained in a single file with three columns `orf|bit_score|gene_family`. **See** [annotation](doc/genome_annotation.md) for formating the raw annotation files.

```bash
python TMarSel/TMarSel.py \
    --input_file    data/wol2/emapper_wol2_example.tsv \
    --output_dir    out/wol2 
```

* KEGG annotations contained in a single file with three columns `orf|bit_score|gene_family`. **See** [annotation](doc/genome_annotation.md) for formating the raw annotation files.

```bash
python TMarSel/TMarSel.py \
    --input_file    data/wol2/kofamscan_wol2_example.tsv \
    --output_dir    out/wol2 
```

### 2\. Annotations of 793 metagenome-assembled genomes (MAGs) from the Earth Microbiome Project

* EggNOG annotations contained multiple (not compressed) files with three columns `orf|bit_score|gene_family`. **See** [annotation](doc/genome_annotation.md) for formating the raw annotation files.

```bash
python TMarSel/TMarSel.py \
    --input_file_or_dir data/emp/eggnog_format \
    --output_dir        out/emp
```

* EggNOG annotations contained in multiple (not compressed) files with raw annotations.

```bash
python TMarSel/TMarSel.py \
    --input_file_or_dir data/emp/eggnog \
    --output_dir        out/emp \
    --database          eggnog \
    --raw_annotations
```

* KEGG annotations contained in multiple files with three columns `orf|bit_score|gene_family`. **See** [annotation](doc/genome_annotation.md) for formating the raw annotation files.

```bash
python TMarSel/TMarSel.py \
    --input_file_or_dir data/emp/kegg_format \
    --output_dir        out/emp
```

* KEGG annotations contained in multiple files with raw annotations.

```bash
python TMarSel/TMarSel.py \
    --input_file_or_dir data/emp/kegg \
    --output_dir        out/emp \
    --database          kegg \
    --raw_annotations
```

## Citation

The current version of **TMarSel** is described in 

* x.x
