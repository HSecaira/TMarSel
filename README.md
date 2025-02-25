# TMarSel

**TMarSel** is a tool for a Tailored Selection of gene families as Markers for microbial phylogenomics.

## Inputs

* Mandatory: `-i` Either a single annotation of ORFs into gene families file OR a file containing a list of annotation file names (one per line).
* Mandatory: `-db` Name of the database used for genome annotation (`kegg` or `eggnog`).
* Mandatory: `-o` Output directory to save the ORFs of each marker.
* Optional: `-input_dir` Directory containing the input files. **Required** if providing a list of file names in `-i`
* Optional: `-k` Number of markers to select. Default is 10.
* Optional: `-th` Threshold for filtering copies of each gene family per genome. Default is `1.0`, retaining only the ORFs with the maximum bit score for each gene family and genome (i.e. less copies). `0.0` reatains all the ORFs (i.e. all copies).
* Optional: `-p` Exponent of the power mean (cost function) used to select markers. Default is 0.0 which yields the optimal combination of markers.
* Optional: `-c` Indicate whether the input file(s) are (.xz) compressed.
* Optional: `-min_markers` Minimum number of markers per genome. Can be a percentage or a number. Genomes with fewer markers than the indicated value are discarded. Default is 1.
* `-h` Display help message.

## Outputs

* `k` files containing the ORFs of each marker gene.
* Protein sequences of each marker (**To implement**)
* Statistics (e.g. number of markers per genome, number of genomes per marker, description of genes; **To implement**)

## Installation

* **Anaconda**
* **pip**
* **GitHub**

## Usage

After installation, type `python TMarSel.py -h` to see all the options.

## Examples

We provide two four examples to showcase the usage of **TMarSel**.

1. **Annotations of 1,510 genomes from the Web of Life 2 database**
    * EggNog annotations contained in a single (compressed) file.
        `python TMarSel.py -i data/wol2/emapper_wol2.annotations.xz -db eggnog -o out/wol2 -c`
    * KEGG annotations contained in a single (compressed) file.
        `python TMarSel.py -i data/wol2/kofamscan_wol2.tsv.xz -db kegg -o out/wol2 -c`
2. **Annotations of 793 metagenome-assembled genomes (MAGs) from the Earth Microbiome Project**
    * EggNog annotations contained in multiple (not compressed) files.
        `python TMarSel.py -i data/emp/mags.txt -input_dir data/emp/eggnog -db eggnog -o out/emp`
    * KEGG annotations contained in multiple (not compressed) files.
        `python TMarSel.py -i data/emp/mags.txt -input_dir data/emp/kegg -db kegg -o out/emp`

## Publication

The current version of **TMarSel** is described in 
    * x.x
