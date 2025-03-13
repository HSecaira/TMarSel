#!/usr/bin/env python3
'''
Select gene families as markers for microbial phylogenomics
Version: 0.1.0
Authors: Henry Secaira and Qiyun Zhu
March 2025
'''

# Import the necessary libraries
import numpy as np
import pandas as pd
import os, sys, argparse
import gzip, bz2, lzma
from scipy.stats import pmean
from tqdm import tqdm


##################################################
# Functions

def load_genome_annotations_single_file(input_data, database, raw_annotations):
    """
    Load genome annotation data from EggNOG or KEGG format.

    Parameters
    ----------
        input_data : str
            Input file
        database : str (optional)
            Either 'eggnog' or 'kegg' to specify the format
        compressed : bool (optional)
            Whether the file is compressed (default is True)
        raw_annotations : bool (optional)
            Whether the file contains raw annotations (default is False)
    
    Returns
    -------
        df : pandas.DataFrame
            Dataframe with columns ['orf', 'bit_score', 'gene_family', 'genome'], indexed by 'orf'.
    """

    # Check if file exists
    if not os.path.exists(input_data):
        raise FileNotFoundError(f"File '{input_data}' does not exist. Please check the filename and path.")

    # Determine the correct open function
    if input_data.endswith('.xz'):
            open_func = lzma.open
    elif input_data.endswith('.gz'):
        open_func = gzip.open
    elif input_data.endswith('.bz2'):
        open_func = bz2.open
    else:
        open_func = open

    # Get file name
    file_name = os.path.basename(input_data)
    # Load data
    tmp = []
    with open_func(input_data, mode = 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            row = line.strip().split('\t')
            # If raw annotations, read the raw data
            if raw_annotations:
                if database == 'eggnog':
                    orf, score, gene = row[0], row[3], row[4].split('|')[0].split('@')[0]
                elif database == 'kegg':
                    if not line.startswith('*'):
                        continue
                    orf, score, gene = row[1], row[4], row[2]
                else:
                    raise ValueError("Invalid database type. Choose 'eggnog' or 'kegg'.")
            else:
                orf, score, gene = row[0], row[1], row[2]

            genome = orf.split('_')[0]
            # Add data to tmp list
            tmp.append([orf, float(score), gene, genome, file_name])

    df = pd.DataFrame(tmp, columns=['orf', 'bit_score', 'gene_family', 'genome', 'file_name'])
    df.set_index('orf', inplace=True)

    return df

def load_genome_annotations_multiple_files(input_data, database, raw_annotations):
    """
    Load genome annotation data from multiple EggNOG or KEGG files.
    
    Parameters
    ----------
        input_data : str
            Directory containing the input files.
        database : str (optional)
            Either 'eggnog' or 'kegg' to specify the format.
        compressed : bool (optional)
            Whether the file is compressed (default is True)
        raw_annotations : bool (optional)
            Whether the file contains raw annotations (default is False)
    
    Returns
    -------
        df : pandas.DataFrame
            Dataframe with columns ['orf', 'bit_score', 'gene_family', 'genome'], indexed by 'orf'.
    """

    # List to store data
    tmp = []
    # Check input directory is not empty
    if not os.path.exists(input_data):
        raise FileNotFoundError(f"Directory '{input_data}' does not exist. Please check the directory name and path.")
    # Iterate over input files
    files = os.listdir(input_data)
    for file_name in files:
        fIn = os.path.join(input_data, file_name)
        # Determine the correct open function
        if file_name.endswith('.xz'):
            open_func = lzma.open
        elif file_name.endswith('.gz'):
            open_func = gzip.open
        elif file_name.endswith('.bz2'):
            open_func = bz2.open
        else:
            open_func = open

        # Load data
        with open_func(fIn, mode = 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                row = line.strip().split('\t')
                # If raw annotations, read the raw data
                if raw_annotations:
                    if database == 'eggnog':
                        orf, score, gene = row[0], row[3], row[4].split('|')[0].split('@')[0]
                    elif database == 'kegg':
                        if not line.startswith('*'):
                            continue
                        orf, score, gene = row[1], row[4], row[2]
                    else:
                        raise ValueError("Invalid database type. Choose 'eggnog' or 'kegg'.")
                else:
                    orf, score, gene = row[0], row[1], row[2]

                genome = orf.split('_')[0]
                # Add data to tmp list
                tmp.append([orf, float(score), gene, genome, file_name])

    # Create a single DataFrame after all files are processed
    df = pd.DataFrame(tmp, columns=['orf', 'bit_score', 'gene_family', 'genome', 'file_name'])
    df.set_index('orf', inplace=True)

    return df

def filter_copies(df, threshold):
    '''
    Filter out copies of gene families per genome based on a threshold.

    Parameters
    -----------
        df : pandas.DataFrame
            DataFrame with columns 'bit_score', 'gene_family', and 'genome', indexed by 'orf'.
        threshold : float
            Float value between 0 and 1 indicating the closeness to the maximum bit score for each gene family and genome.
            Lower values (e.g. 0.0) retains all ORFs, whereas higher values (e.g. 1.0) retains only the ORFs with the highest bit score.
    Returns
    --------
        filtered_df : pandas.DataFrame
            DataFrame with columns 'bit_score', 'gene_family', and 'genome', indexed by 'orf'.
    '''
    # Check if threshold is a valid value
    if not 0 <= threshold <= 1:
        raise ValueError('Threshold should be a value between 0 and 1.')
    # Calculate the maximum score for each genome and gene family combination
    max_scores = df.groupby(['genome', 'gene_family'])['bit_score'].transform('max')
    # Keep rows where the score is greater than or equal to max_score * threshold
    filtered_df = df[df['bit_score'] >= max_scores * threshold]

    return filtered_df

def get_edges(filtered_df):
    '''
    Get the mapping between genome and gene families.

    Parameters
    -----------
        filtered_df : pandas.DataFrame
            DataFrame with columns 'bit_score', 'gene_family', and 'genome', indexed by 'orf'.
    Returns
    --------
        edges_genes : numpy.ndarray
            Array of gene families.
        edges_genomes : numpy.ndarray
            Array of genomes.
    '''
    edges_genomes = filtered_df['genome'].values
    edges_genes = filtered_df['gene_family'].values

    return edges_genes, edges_genomes

def build_copy_number_matrix(edges_genes, edges_genomes):
    '''
    Build a copy number matrix from the arrays of gene families and genomes.

    Parameters
    -----------
        edges_genes : numpy.ndarray
            Array of gene families.
        edges_genomes : numpy.ndarray
            Array of genomes.
    Returns
    --------
        adj : numpy.ndarray
            A matrix where rows represent genes and columns represent genomes. 
            Each entry represents the number of copies of a gene family in a genome.
        genomes : numpy.ndarray
            Array of unique genomes.
        genes : numpy.ndarray
            Array of unique gene families.
    '''
    # Get unique elements
    genes, genes_indices = np.unique(edges_genes, return_inverse = True)
    genomes, genomes_indices = np.unique(edges_genomes, return_inverse = True)
    # Calculate bin counts for each combination of indices
    counts = np.bincount(genes_indices * len(genomes) + genomes_indices, minlength = len(genes) * len(genomes))
    # Reshape counts as adjacency matrix
    adj = counts.reshape(len(genes), len(genomes))
    return adj, genomes, genes


def remove_genes(adj):
    """
    Remove genes that are present in less than 4 genomes.
    
    Parameters
    -----------
        adj : numpy.ndarray
            A matrix where rows represent genes and columns represent genomes.    
    Returns
    --------
        remove : numpy.ndarray
            Indices of rows (genes) to be removed.
    """
    # remove = np.array([i for i in range(len(adj)) if np.count_nonzero(adj[i]) < 4])
    remove = np.where(np.count_nonzero(adj, axis = 1) < 4)[0]
    return remove

def get_genomes_to_keep(adj, k, markers_index, min_markers, genomes):
    """
    Filter genomes based on a minimum number of markers.

    Parameters
    -----------
        adj : numpy.ndarray 
            A matrix where rows represent genes and columns represent genomes.
        k : int
            The total number of markers selected.
        markers_index : numpy.ndarray
            Indices of rows corresponding to selected markers
        min_markers : float or int
            Threshold for filtering genomes. 
            If a float (0-1), it is treated as a percentage of k.
            If an int, it is used as an absolute threshold.
    Returns
    --------
        numpy.ndarray
            Array of genomes that meet the threshold of minimum number of markers per genome
    """
    threshold = min_markers * k if isinstance(min_markers, float) and 0 <= min_markers <= 1 else min_markers
    return genomes[adj[markers_index].sum(axis = 0) >= threshold]

def greedy_power_mean_sample_final(data, k, p, pseudocount):
    """
    Select k rows from a matrix such that the selection criterion by column is maximized.

    Parameters
    ----------
    data : numpy.ndarray
        A matrix where rows represent genes and columns represent genomes.
    k : int
        Number of rows (markers) to select.
    p : int
        Exponent of generalized power mean
    Returns
    -------
    numpy.ndarray
        Indices of rows corresponding to selected markers
    """

    n, m = data.shape

    # Matrix is empty
    if n == 0 or m == 0:
        raise ValueError(f'Matrix is empty!')

    # Matrix contains only zeroes
    if (data == 0).all():
        raise ValueError(f'Matrix only contains 0\'s')

    if k >= n:
        raise ValueError(f'k should be smaller than {n}')
    
    # Add pseudocount
    data = data + pseudocount

    # Cumulative gene counts
    counts = np.zeros(m, dtype = int)

    # Gene indices in original data matrix
    indices = np.arange(n)

    # Indices of selected genes
    selected = []

    # Select k genes iteratively and display progress bar
    with tqdm(total = k, desc = f'Selection progress') as pbar:
        for i in range(k):
            # calculate counts after adding each gene
            sums_ = counts + data

            # Select a gene that maximizes the power mean gene count per genome, using the cumulative matrix
            if isinstance(p, int) or isinstance(p, np.int64): 
                choice = pmean(sums_, int(p), axis = 1).argmax()
            elif p == 'min':
                choice = sums_.min(axis = 1).argmax()
            elif p == 'max':
                choice = sums_.max(axis = 1).argmax()
            else:
                raise ValueError(f'Invalid p: {p}.')

            # Append index of selected gene
            selected.append(indices[choice])

            # Update per-species gene counts
            counts = sums_[choice]

            # Remove selected gene from data matrix
            data = np.delete(data, choice, axis = 0)

            # Remove selected gene from indices
            indices = np.delete(indices, choice)

            # Update progress bar
            pbar.update(1)

    return np.array(selected)

def reformat_column_counts(columns_series):
    columns_counts = columns_series.value_counts().items()
    return ';'.join(f"{item}:{count}" for item, count in columns_counts)


def save_marker_orfs(markers_index, genes_mod, filtered_df, genomes_to_keep, output_dir):
    """
    Save statistics ORFs for each marker gene to individual text files.

    Parameters
    ----------
        markers_index : numpy.ndarray
            Indices of marker genes.
        genes_mod : numpy.ndarray
            Gene family names
        filtered_df : pandas.DataFrame)
            DataFrame with columns 'bit_score', 'gene_family', and 'genome', indexed by 'orf'.
        genomes_to_keep : numpy.ndarray
            Genomes to keep based on the minimum number of markers per genome.
        output_dir : str
            Directory to save the ORF files.
    
    Returns:
        None
    """
    
    # Ensure the output directories exists
    os.makedirs(f'{output_dir}/orfs', exist_ok = True)
    os.makedirs(f'{output_dir}/statistics', exist_ok = True)
    # Get marker names
    markers_names = genes_mod[markers_index]
    # Get ORFs of markers and genomes
    orfs_markers = filtered_df[filtered_df['gene_family'].isin(markers_names) & filtered_df['genome'].isin(genomes_to_keep)].copy() #make a copy to avoid SettingWithCopyWarning
    
    # Save number of markers per genome
    fOut_markers_per_genome = os.path.join(output_dir, 'statistics/number_of_markers_per_genome.tsv')
    number_markers_per_genome = orfs_markers.groupby('genome').agg(
        number_of_different_markers=('gene_family', 'nunique'),
        total_number_of_markers=('gene_family', 'count'),
        details = ('gene_family', reformat_column_counts)
    )
    number_markers_per_genome.to_csv(fOut_markers_per_genome, sep = '\t', index = True)

    # Save number of genomes per marker
    fOut_genomes_per_marker = os.path.join(output_dir, 'statistics/number_of_genomes_per_marker.tsv')
    number_genomes_per_marker = orfs_markers.groupby('gene_family').agg(
        number_of_genomes=('genome', 'nunique'),
        details=('genome', reformat_column_counts)
    )
    number_genomes_per_marker.index.name = 'marker'
    number_genomes_per_marker.to_csv(fOut_genomes_per_marker, sep = '\t', index = True)

    # Save ORFs for each marker gene and show progress bar
    with tqdm(total = len(markers_names), desc = 'Saving progress') as pbar:
        for marker in markers_names:
            file_path_orfs = os.path.join(output_dir, f'orfs/{marker}.txt')
            # Get ORFs for the marker gene
            orfs = orfs_markers[orfs_markers['gene_family'] == marker]
            # Write ORFs, genome, and file_name to file
            orfs.to_csv(file_path_orfs, sep = '\t', columns = ['genome', 'file_name'], header = False, index = True, mode = 'w')
            # Update progress bar
            pbar.update(1)

def print_tmarsel():
    tmarsel = [
        r" _____ __  __            ____       _ ",
        r"|_   _|  \/  | __ _ _ __/ ___|  ___| |",
        r"  | | | |\/| |/ _` | '__\___ \ / _ \ |",
        r"  | | | |  | | (_| | |   ___) |  __/ |",
        r"  |_| |_|  |_|\__,_|_|  |____/ \___|_|"
    ]
    
    for line in tmarsel:
        print(line)


def main(argv = None):

    # Print banner
    print_tmarsel()

    args_parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
    description = f'TMarSel: Tailored Marker Selection of gene families for microbial phylogenomics\nVersion: 0.1.0\nBasic usage: python TMarSel.py -i input_file -k markers -o output_dir\nType python TMarSel.py -h for help')
    args_parser.add_argument('-i', '--input_file_or_dir', type = str, required = True,
     help = '[required] File containing the genome annotations of ORFs into gene families.\nEither a single annotation file OR a file containing a list of annotation file names (one per line)')
    args_parser.add_argument('-o', '--output_dir', type = str, required = True, 
     help = 'Output directory for the ORFs and statistics of each marker')
    args_parser.add_argument('-raw', '--raw_annotations', action = 'store_true', 
    help = '[required] IF input file contains raw annotations')
    args_parser.add_argument('-db', '--database', type = str, 
    help = '[required] IF input contains the raw annotations set in "-raw".\nDatabase used for genome annotation: eggnog or kegg')
    args_parser.add_argument('-k', '--markers', type = int, default = 50, 
    help = 'Number of markers to select (default is 50)')
    args_parser.add_argument('-min_markers', '--min_number_markers_per_genome', type = lambda x: float(x) if '.' in x else int(x), default = 1,
     help = 'Minimum number of markers per genome. Can be a percentage or a number (default is 1).\nGenomes with fewer markers than the indicated value are discarded')
    args_parser.add_argument('-th', '--threshold', type = float, default = 1.0, 
     help = 'Threshold for filtering copies of each gene family per genome (default is 1.0)\nRetain the ORFs within "threshold" of the maximum bit score for each gene family and genome.\nLower values (e.g. 0.0) retains all ORFs, whereas higher values (e.g. 1.0) retains only the ORF with the highest bit score')
    args_parser.add_argument('-p', '--exponent', type = int, default = 0, 
     help = 'Exponent of the power mean (cost function) used to select markers (default is 0).\nWe recommend not changing this value unless you are familiar with the method.')
    args = args_parser.parse_args(argv)

    # Get names
    input_data = args.input_file_or_dir
    database = args.database
    threshold = args.threshold
    k = args.markers
    p = args.exponent
    output_dir = args.output_dir
    min_markers = args.min_number_markers_per_genome
    raw_annotations = args.raw_annotations

    # Run TMarSel
    run_TMarSel(input_data, database, threshold, k, p, output_dir, min_markers, raw_annotations)
    
    return 0

def run_TMarSel(input_data, database, threshold, k, p, output_dir, min_markers, raw_annotations):

    # Load data
    if not os.path.exists(input_data):
        raise FileNotFoundError(f"'{input_data}' does not exist. Please check the filename and/or path of directory")
    else:
        if os.path.isfile(input_data):
            print(f'Loading genome annotation data from a single file with "-raw" (annotations) set to {raw_annotations}...')
            df = load_genome_annotations_single_file(input_data, database, raw_annotations)
        elif os.path.isdir(input_data):
            file_names = os.listdir(input_data)
            print(f'Loading genome annotation data from {len(file_names)} files with "-raw" (annotations) set to {raw_annotations}...')
            df = load_genome_annotations_multiple_files(input_data, database, raw_annotations)

    print(f'\tAnnotation file has : {df.shape[0]} ORFs assigned to a gene family\n')
    # Filter copies
    print(f'Filtering copies with threshold: {threshold}...')
    filtered_df = filter_copies(df, threshold)
    print(f'\t{df.shape[0] - filtered_df.shape[0]} ORFs were filtered out, leaving {filtered_df.shape[0]} ORFs for further analysis\n')
    # Get edges
    edges_genes, edges_genomes = get_edges(filtered_df)
    # Build copy number matrix
    print('Building matrix for marker selection...')
    adj, genomes, genes = build_copy_number_matrix(edges_genes, edges_genomes)
    print(f'\tMatrix has {adj.shape[0]} gene families from {adj.shape[1]} genomes/MAGs\n')
    # Remove genes
    print(f'Removing gene families present in less than 4 genomes/MAGs...')
    remove = remove_genes(adj)
    genes_mod = np.delete(genes, remove, axis = 0)
    adj_mod = np.delete(adj, remove, axis = 0)
    print(f'\tMatrix now has {adj_mod.shape[0]} gene families from {adj_mod.shape[1]} genomes/MAGs\n')
    # Select markers
    print(f'Selecting {k} markers with parameter p = {p}...')
    markers_index = greedy_power_mean_sample_final(data = adj_mod, k = k, p = p, pseudocount = 0.1)
    # Remove genomes with less than min_markers markers
    print(f'\nRemoving genomes/MAGs with less than {min_markers} markers...')
    genomes_to_keep = get_genomes_to_keep(adj_mod, k, markers_index, min_markers, genomes)
    print(f'\t{len(genomes_to_keep)} fit the criteria above\n')
    # Save ORFs of each marker
    print('Saving statistics and ORFs for each marker gene...')
    save_marker_orfs(markers_index, genes_mod, filtered_df, genomes_to_keep, output_dir)
    print('Done!')


##################################################

if __name__ == '__main__':
    status = main()
    sys.exit(status)











