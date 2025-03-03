# TMarSel

**TMarSel** is a tool for a Tailored Selection of gene families as Markers for microbial phylogenomics.

## Table of contents

* Inputs[#inputs]
* Outputs[#outputs]
* Installation[#installation]
* Basic usage[#basic-usage]
* Examples[#examples]
* Citation[#citation]

## Inputs

* `-h` Display help message.
* -i` [**required**] Either a single annotation of ORFs into gene families file OR a file containing a list of annotation file names (one per line).\nFile(s) contains three columns: `orf|bit_score|gene_family`
* `-o` [**required**] Output directory to save the ORFs and statistics of each marker.
* `-input_dir` [**required**] IF providing a list of files in `-i`.\nDirectory containing the input files.
* `-raw` [**required**] IF input file(s) contain raw annotations from functional databases, which contains multiple columns depending on the database.
* `-db` [**required**] IF input contains the raw annotations set in `-raw`.\nName of the database used for genome annotation (`eggnog` or `kegg`).
* `-k` Number of markers to select (default is 50).
* `-min_markers` Minimum number of markers per genome. Can be a percentage or a number (default is 1).\nGenomes with fewer markers than the indicated value are discarded.
* `-th` Threshold for filtering copies of each gene family per genome (default is 1.0)\nRetain the ORFs within `-th` of the maximum bit score for each gene family and genome. Lower values (e.g. 0.0) retains all ORFs, whereas higher values (e.g. 1.0) retains only the ORF with the highest bit score.
* `-p` Exponent of the power mean (cost function) used to select markers (default is 0.0).\nWe recommend not changing this value unless you are familiar with the method\nDefault value yields the optimal combination of markers..
* `-c` Indicate whether the input file(s) is (.xz) compressed.

## Outputs

* `k` files containing the ORFs, genome and file of origin for each marker (see below). Files are saved to `./output_dir/orfs`. 

| orf | genome | file |
| --- | --- | --- |
| G000006605_1748 | G000006605 | kofamscan_wol2_example.tsv |
| G000006725_378 | G000006725 | kofamscan_wol2_example.tsv |
| ... | ... | ... |

* Statistics. Files are saved to `./output_dir/statistics`

    * Number of markers per genome (see below). A given marker can contain more than one ORF per genome, therefore we provide the number of different markers (`k`) and the total number of markers. The `details` column is `;` separatend. Each item indicates the marker name and the number of ORFs (i.e. copies) in the genome.

| genome | number_of_different_markers | total_number of markers | details |
| --- | --- | --- | --- |
| G000006605 | 10 | 10 | K01889:1;K01866:1;K01872:1;K02600:1;K01867:1;K07497:1;K01409:1;K02982:1;K02358:1;K02867:1 |
| G000006725 | 9 | 10 | K02358:2;K01872:1;K01866:1;K02600:1;K01867:1;K01409:1;K01889:1;K02982:1;K02867:1
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

After installation, type `python TMarSel.py -h` to see all the options.

## Examples

We provide multiple examples to showcase the usage of **TMarSel**.

1. **Annotations of 1,510 genomes from the Web of Life 2 database**

    * EggNog annotations contained in a single file with three columns `orf|bit_score|gene_family`.

```bash
python TMarSel/TMarSel.py \
--input_file    data/wol2/emapper_wol2_example.tsv \
--output_dir out/wol2 
```

    * EggNog annotations contained in a single (xz compressed) file with raw annotations.

```bash
python TMarSel/TMarSel.py \
--input_file    data/wol2/emapper_wol2.annotations.xz \
--output_dir    out/wol2 \
--compressed \
--raw_annotations \
--database  eggnog
```
    * KEGG annotations contained in a single file with three columns `orf|bit_score|gene_family`.

```bash
python TMarSel/TMarSel.py \
--input_file    data/wol2/kofamscan_wol2_example.tsv \
--output_dir out/wol2 
```
    * KEGG annotations contained in a single (xz compressed) file with raw annotations.

```bash
python TMarSel/TMarSel.py \
--input_file    data/wol2/kofamscan_wol2.tsv.xz \
--output_dir    out/wol2 \
--compressed \
--raw_annotations \
--database  kegg
```

2. **Annotations of 793 metagenome-assembled genomes (MAGs) from the Earth Microbiome Project**

    * EggNog annotations contained in multiple files with raw annotations.

```bash
python TMarSel/TMarSel.py \
--input_file    data/emp/mags.txt \
--input_dir_files data/emp/eggnog \
--output_dir    out/emp \
--raw_annotations \
--database  eggnog
```

    * KEGG annotations contained in multiple (not compressed) files.

```bash
python TMarSel/TMarSel.py \
--input_file    data/emp/mags.txt \
--input_dir_files data/emp/kegg \
--output_dir    out/emp \
--raw_annotations \
--database  kegg
```


## Citation

The current version of **TMarSel** is described in 

    * x.x
