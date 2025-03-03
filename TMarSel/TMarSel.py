#!/usr/bin/env python3
'''
Select gene families as markers for microbial phylogenomics
Version: 0.1.0
Authors: Henry Secaira and Qiyun Zhu
February 2025
'''

# Import the necessary libraries
import numpy as np
import pandas as pd
import lzma, argparse, os, sys, pyfiglet
from scipy.stats import pmean
from tqdm import tqdm


##################################################
# Functions

def load_genome_annotations_single_file(fIn, database, compressed, raw_annotations):
    """
    Load genome annotation data from EggNOG or KEGG format.
    
    Parameters:
        fIn (str): Input file path.
        database (str): Either 'eggnog' or 'kegg' to specify the format.
        compressed (bool, optional): Whether the file is compressed (default is True).
    
    Returns:
        pd.DataFrame: A dataframe with columns ['orf', 'bit_score', 'gene_family', 'genome'], indexed by 'orf'.
    """

    # Check if file exists
    if not os.path.exists(fIn):
        raise FileNotFoundError(f"File '{fIn}' does not exist. Please check the filename and path.")

    # Determine the correct open function
    try:
        with lzma.open(fIn, mode = 'rt') as f:
            f.read(1)  # Try reading a character to check compression
        if not compressed:
            raise ValueError(
                f'File "{fIn}" appears to be compressed.\nTry running the command with "-c" to parse compressed files.'
            )
        open_func = lzma.open
    except lzma.LZMAError:
        if compressed:
            raise ValueError(f'File "{fIn}" is not a valid LZMA-compressed file.')
        open_func = open

    # # Check if file is a single annotation file or a list of annotation files
    # with open_func(fIn, mode = 'rt') as f:
    #     line = f.readline()
    #     if len(line.split('\t')) == 1:
    #         raise ValueError('Please provide a single annotation file, not a list of files.')

    # Load data
    tmp = []
    with open_func(fIn, mode='rt') as f:
        # Get file name
        file_name = os.path.basename(fIn)
        # Iterate over lines in the file
        for line in f:
            if not line.startswith('#'):
                row = line.strip().split('\t')

                # If raw annotations, read the raw dara
                if raw_annotations:
                    orf, score = row[0], row[3]
                    genome = orf.split('_')[0]

                    if database == 'eggnog':
                        # Extract most basal OG
                        gene = row[4].split('|')[0].split('@')[0]
                    elif database == 'kegg':
                        # Extract KO identifier
                        gene = row[1]
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

def load_genome_annotations_multiple_files(file_names, database, compressed, raw_annotations):
    """
    Load genome annotation data from multiple EggNOG or KEGG files.
    
    Parameters:
        file_names (list of str): List of input file paths.
        database (str): Either 'eggnog' or 'kegg' to specify the format.
        compressed (bool): Whether the files are compressed.
    
    Returns:
        pd.DataFrame: A dataframe with columns ['orf', 'bit_score', 'gene_family', 'genome'], indexed by 'orf'.
    """

    tmp = []
    for fIn in file_names:
        if database == 'eggnog':
            fIn = fIn + '.emapper.annotations'
        elif database == 'kegg':
            fIn = fIn + '.tsv'
        else:
            raise ValueError("Invalid database type. Choose 'eggnog' or 'kegg'.")
        # Check if file exists
        if not os.path.exists(fIn):
            raise FileNotFoundError(f"File '{fIn}' does not exist. Please check the filename and path.")

        # Determine the correct open function
        try:
            with lzma.open(fIn, mode = 'rt') as f:
                f.read(1)  # Try reading a character to check compression
            if not compressed:
                raise ValueError(
                    f'File "{fIn}" appears to be compressed.\nTry running the command with "-c" to parse compressed files.'
                )
            open_func = lzma.open
        except lzma.LZMAError:
            if compressed:
                raise ValueError(f'File "{fIn}" is not a valid LZMA-compressed file.')
            open_func = open

        # Get file name
        file_name = os.path.basename(fIn)

        # Load data into tmp list
        with open_func(fIn, mode = 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    row = line.strip().split('\t')
                    # If raw annotations, read the raw dara
                    if raw_annotations:
                        orf, score = row[0], row[3]
                        genome = orf.split('_')[0]

                        if database == 'eggnog':
                            gene = row[4].split('|')[0].split('@')[0]  # Extract most basal OG
                        elif database == 'kegg':
                            gene = row[1]  # Extract KO identifier
                        else:
                            raise ValueError("Invalid database type. Choose '-db eggnog' or '-db kegg'.")
                    else:
                        orf, score, gene = row[0], row[1], row[2]
                        genome = orf.split('_')[0]

                    tmp.append([orf, float(score), gene, genome, file_name])

    # Create a single DataFrame after all files are processed
    df = pd.DataFrame(tmp, columns=['orf', 'bit_score', 'gene_family', 'genome', 'file_name'])
    df.set_index('orf', inplace=True)

    return df

def filter_copies(df, threshold):
    # Check if threshold is a valid value
    if not 0 <= threshold <= 1:
        raise ValueError('Threshold should be a value between 0 and 1.')
    # Calculate the maximum score for each genome and gene family combination
    max_scores = df.groupby(['genome', 'gene_family'])['bit_score'].transform('max')
    # Keep rows where the score is greater than or equal to max_score * threshold
    filtered_df = df[df['bit_score'] >= max_scores * threshold]

    return filtered_df

def get_edges(filtered_df):
    edges_genomes = filtered_df['genome'].values
    edges_genes = filtered_df['gene_family'].values

    return edges_genes, edges_genomes

def build_copy_number_matrix(edges_genes, edges_genomes):
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
    
    Parameters:
        adj (numpy.ndarray): A matrix where rows represent genes and columns represent genomes.
    
    Returns:
        numpy.ndarray: Indices of rows (genes) to be removed.
    """
    # remove = np.array([i for i in range(len(adj)) if np.count_nonzero(adj[i]) < 4])
    remove = np.where(np.count_nonzero(adj, axis = 1) < 4)[0]
    return remove

def get_genomes_to_keep(adj, k, markers_index, min_markers, genomes):
    """
    Filter genomes based on the number of markers.

    Parameters:
        adj (numpy.ndarray): A 2D matrix where rows represent genes and columns represent genomes.
        k (int): The total number of markers.
        markers_index (array-like): Indices of rows corresponding to marker genes.
        min_markers (float or int): Threshold for filtering genomes. 
                                             If a float (0-1), it is treated as a percentage of k.
                                             If an int, it is used as an absolute threshold.
    Returns:
        numpy.ndarray: Array of genomes that meet the threshold.
    """
    threshold = min_markers * k if isinstance(min_markers, float) and 0 <= min_markers <= 1 else min_markers
    return genomes[adj[markers_index].sum(axis = 0) >= threshold]

def greedy_power_mean_sample_final(data, k, p, pseudocount):
    """Select k rows from a matrix such that the selection criterion by column is maximized.

    Parameters
    ----------
    data : ndarray (2D)
        Input data matrix.
    k : int
        Number of rows to select.
    p : int
        Exponent.

    Returns
    -------
    ndarray (1D)
        Row indices selected in order.
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

    Parameters:
        markers_index (array-like): Indices of marker genes.
        genes_mod (array-like): List of gene names corresponding to indices.
        filtered_df (pd.DataFrame): DataFrame containing gene_family and genome information.
        genomes_to_keep (array-like): List of genomes to keep.
        output_dir (str, optional): Directory to save the ORF files. Default is './orfs'.
    
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


def save_statistics(data, genomes, genes_mod, markers_index, filtered_df, output_dir):
    '''
    '''

    # Ensure the output directory exists
    os.makedirs(f'{output_dir}/statistics', exist_ok = True)

    # Save number of markers per genome


def main(argv = None):

    # Print welcome message
    ascii_banner = pyfiglet.figlet_format("TMarSel")
    print(ascii_banner)

    args_parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
    description = f'TMarSel: Tailored Selection gene families as Markers for microbial phylogenomics\nVersion: 0.1.0\nAuthors: Henry Secaira and Qiyun Zhu\nFebruary 2025\nBasic usage: python TMarSel.py -i input_file -k markers -o output_dir\nType python TMarSel.py -h for help')
    args_parser.add_argument('-i', '--input_file', type = str, required = True,
     help = '[required] File containing the genome annotations of ORFs into gene families.\nEither a single annotation file OR a file containing a list of annotation file names (one per line)')
    args_parser.add_argument('-o', '--output_dir', type = str, required = True, 
     help = 'Output directory for the ORFs and statistics of each marker')
    args_parser.add_argument('-input_dir', '--input_dir_files', type = str,
     help = '[required] IF providing a list of file names in "-i".\nDirectory containing the input files')
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
    args_parser.add_argument('-c', '--compressed', action = 'store_true', 
     help = 'Input file(s) is (.xz) compressed')
    args = args_parser.parse_args(argv)

    # Get names
    fIn = args.input_file
    database = args.database
    threshold = args.threshold
    k = args.markers
    p = args.exponent
    output_dir = args.output_dir
    compressed = args.compressed
    min_markers = args.min_number_markers_per_genome
    input_dir = args.input_dir_files
    raw_annotations = args.raw_annotations

    # Run TMarSel
    run_TMarSel(fIn, database, threshold, k, p, output_dir, compressed, min_markers, input_dir, raw_annotations)
    
    return 0

def run_TMarSel(fIn, database, threshold, k, p, output_dir, compressed, min_markers, input_dir, raw_annotations):
    # Load data
    # In case of multiple annotation files
    if input_dir:
        with open(fIn, 'r') as f:
            file_names = [f'{input_dir}/{line.strip()}' for line in f]
        print(f'Loading genome annotation data from {len(file_names)} files with "-raw" (annotations) set to {raw_annotations}...')
        df = load_genome_annotations_multiple_files(file_names, database, compressed, raw_annotations)
    else:
        print(f'Loading genome annotation data from a single file with "-raw" (annotations) set to {raw_annotations}...')
        df = load_genome_annotations_single_file(fIn, database, compressed, raw_annotations)

    # df = load_genome_annotations(fIn, database, compressed)    
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











