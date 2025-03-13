# This file contains the functions to be used in the notebooks for the project

import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import anndata as ad
import pandas as pd
import scipy.sparse
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from matplotlib_venn import venn2
from typing import List, Optional, Union
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

homeo_genes = ['Sall1', 'Cx3cr1', 'P2ry12', 'Olfml3', 'Tmem119', 'Cd68', 'Itgam', 'Cst3']
DAM_genes = ['Apoe', 'B2m', 'Ctsb', 'Lpl', 'Cst7', 'Tyrobp', 'Trem2', 'Cd9', 'Itgax', 'Cd63', 'Fth1', 'Fabp5', 'Spp1', 'Axl', 'Mertk', 'Mdk', 'Ptn', 'Lgals3']
Mapk_genes = ['Syk', 'Prkca', 'Mef2c', 'Elk4', 'Trp53', 'Mertk', 'Axl', 'Mapk14', 'Mapk1', 'Gadd45a',]
Nfkb_genes = ['Nfkb1','Nfkb2', 'Tlr4', 'Nfkbia', 'Cd40', 'Tlr4', 'Rela', 'Relb', 'Irf1', 'Irf3', 'Irf4']
Macroph_genes = ['Mrc1', 'Mrc2', 'Ccr2', 'Ly6c2', 'Lyz2', 'Vim', 'Ifi204', 'S100a10', 'Msrb1']
DC_genes = ['Cd14', 'Cd1a', 'Cd209', 'Itgam', 'Cd207', 'Fcer1', '33D1', 'eyfp', 'Flt3', 'Csf1r', 'Id2', 'Clec4k', 'Cd8a']
TIM_genes = ['Malat1', 'Tmx4', 'Jund', 'Btg2', 'Fosb', 'Fos','Jun', 'Klf4', 'Klf2', 'Gm34455']



def summarize_adata(adata, mt_gene_prefix="mt-", ribo_gene_prefixes=("Rps", "Rpl"), min_counts=500, min_genes=200, min_cells=3):
    # Calculate the total number of cells and genes
    total_cells = adata.n_obs
    total_genes = adata.n_vars

    # Calculate number of cells with fewer than 500 UMI counts
    cells_with_500_or_more = (adata.X.sum(axis=1) >= min_counts).sum()
    cells_with_less_than_500 = total_cells - cells_with_500_or_more
    print(f"Number of cells with fewer than {min_counts} UMI counts: {cells_with_less_than_500}")

    # Calculate number of cells with fewer than 200 genes
    cells_with_200_or_more_genes = (adata.X > 0).sum(axis=1) >= min_genes
    cells_with_less_than_200_genes = total_cells - cells_with_200_or_more_genes.sum()
    print(f"Number of cells with fewer than {min_genes} genes: {cells_with_less_than_200_genes}")

    # Calculate number of genes expressed in fewer than 3 cells
    genes_with_3_or_more_cells = (adata.X > 0).sum(axis=0) >= min_cells
    genes_with_less_than_3_cells = total_genes - genes_with_3_or_more_cells.sum()
    print(f"Number of genes expressed in fewer than {min_cells} cells: {genes_with_less_than_3_cells}")

    # Annotate mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(mt_gene_prefix)

    # Annotate ribosomal genes with both prefixes
    adata.var["ribo"] = adata.var_names.str.startswith(ribo_gene_prefixes[0]) | adata.var_names.str.startswith(ribo_gene_prefixes[1])

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    # Summarize mitochondrial and ribosomal genes
    mt_genes = adata.var["mt"].sum()
    ribo_genes = adata.var["ribo"].sum()

    pct_counts_mt = adata.obs["pct_counts_mt"].mean()
    pct_counts_ribo = adata.obs["pct_counts_ribo"].mean()

    total_counts = adata.obs["total_counts"].sum()
    total_counts_mt = (adata[:, adata.var["mt"]].X.sum())
    total_counts_ribo = (adata[:, adata.var["ribo"]].X.sum())

    avg_counts_mt_per_cell = total_counts_mt / total_cells
    avg_counts_ribo_per_cell = total_counts_ribo / total_cells

    avg_pct_counts_mt_per_cell = adata.obs["pct_counts_mt"].mean()
    avg_pct_counts_ribo_per_cell = adata.obs["pct_counts_ribo"].mean()

    print(f"Number of mitochondrial genes: {mt_genes}")
    print(f"Number of ribosomal genes: {ribo_genes}")
    print(f"Percentage of total reads from mitochondrial genes: {pct_counts_mt:.2f}%")
    print(f"Percentage of total reads from ribosomal genes: {pct_counts_ribo:.2f}%")
    print(f"Average reads per cell from mitochondrial genes: {avg_counts_mt_per_cell:.2f}")
    print(f"Average reads per cell from ribosomal genes: {avg_counts_ribo_per_cell:.2f}")
    print(f"Average percentage of mitochondrial reads per cell: {avg_pct_counts_mt_per_cell:.2f}%")
    print(f"Average percentage of ribosomal reads per cell: {avg_pct_counts_ribo_per_cell:.2f}%")

    # Print current status of the AnnData object
    print(f"Current Anndata has {adata.n_obs} cells and {adata.n_vars} genes, with a total amount of {total_counts} UMI counts")


def QC_filter_adata(adata, mt_threshold=5, ribo_threshold=5, min_counts=500, min_genes=200, min_cells=3):
    """
    Perform quality control filtering on AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    mt_threshold (float): Threshold for the percentage of mitochondrial genes. Default is 5.
    ribo_threshold (float): Threshold for the percentage of ribosomal genes. Default is 5.
    min_counts (int): Minimum number of counts required for a cell. Default is 500.
    min_genes (int): Minimum number of genes required for a cell. Default is 200.
    min_cells (int): Minimum number of cells required for a gene to be kept. Default is 3.

    Returns:
    AnnData: The filtered AnnData object.
    """
    
    # Annotate the group of mitochondrial and ribosomal genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith("Rps") | adata.var_names.str.startswith("Rpl")

    #Print number of cells and genes before filtering
    print(f"Initial AnnData has {adata.n_obs} cells and {adata.n_vars} genes")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    # Print the number of cells with more than mt_threshold percent mitochondrial genes
    cells_with_more_than_mt_threshold = (adata.obs["pct_counts_mt"] > mt_threshold).sum()
    print(f"{cells_with_more_than_mt_threshold} of {adata.n_obs} cells contain more than {mt_threshold}% mitochondrial genes. Will be filtered")

    # Print the number of cells with less than ribo_threshold percent ribosomal genes
    cells_with_less_than_ribo_threshold = (adata.obs["pct_counts_ribo"] < ribo_threshold).sum()
    print(f"{cells_with_less_than_ribo_threshold} of {adata.n_obs} cells contain less than {ribo_threshold}% ribosomal genes. Will be filtered")

    # Filter cells with more than mt_threshold percent mitochondrial genes
    adata = adata[adata.obs.pct_counts_mt < mt_threshold, :]

    # Filter cells with less than ribo_threshold percent ribosomal genes
    adata = adata[adata.obs.pct_counts_ribo >= ribo_threshold, :]

    # Print the number of cells with fewer than min_counts counts to be filtered
    cells_with_less_than_min_counts = (adata.X.sum(axis=1) < min_counts).sum()
    print(f"{cells_with_less_than_min_counts} of {adata.n_obs} cells have fewer than {min_counts} counts. Will be filtered")

    # Apply additional filtering based on counts, genes, and cells
    sc.pp.filter_cells(adata, min_counts=min_counts)

    # Print the number of cells with fewer than min_genes genes to be filtered
    cells_with_less_than_min_genes = (adata.X > 0).sum(axis=1) < min_genes
    print(f"{cells_with_less_than_min_genes.sum()} of {adata.n_obs} cells have fewer than {min_genes} genes. Will be filtered")

    sc.pp.filter_cells(adata, min_genes=min_genes)

    # Print the number of genes expressed in fewer than min_cells cells to be filtered
    genes_with_less_than_min_cells = (adata.X > 0).sum(axis=0) < min_cells
    print(f"{genes_with_less_than_min_cells.sum()} of {adata.n_vars} genes are expressed in fewer than {min_cells} cells. Will be filtered")

    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Print the total number of cells and genes after filtering
    print(f"Filtered AnnData has {adata.n_obs} cells and {adata.n_vars} genes")

    # Recalculate QC metrics for the filtered data
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    return adata

def visualize_QC_measures(adata, title="QC Measures"):
    """
    Visualize QC measures in violin plots.

    Parameters:
    adata (AnnData or dict of AnnData): The input AnnData object or a dictionary of AnnData objects.
    title (str): The title of the plot. Default is "QC Measures".
    """
    if isinstance(adata, dict):
        # Handling dictionary of AnnData objects
        print(f"Visualizations {title} for multiple datasets:")
        
        # Combined violin plots
        combined_data = adata[list(adata.keys())[0]].concatenate(*[adata[key] for key in list(adata.keys())[1:]], batch_key="dataset", batch_categories=list(adata.keys()))
        sc.pl.violin(combined_data, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"], groupby="dataset", jitter=0.4, multi_panel=True, stripplot=False)
        
        # Scatter plots for each dataset
        fig, axes = plt.subplots(1, len(adata), figsize=(15, 5))
        if len(adata) == 1:
            axes = [axes]
            
            sc.pl.scatter(single_adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", ax=ax, show=False)
            ax.set_xlim(0, 40000)
            ax.set_ylim(0, 6000)
            ax.grid(True)
            ax.set_title(f"{adata_key} - pct_counts_mt")
            ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1))
        plt.tight_layout()
        plt.show()

        fig, axes = plt.subplots(1, len(adata), figsize=(15, 5))
        if len(adata) == 1:
            axes = [axes]
        for ax, (adata_key, single_adata) in zip(axes, adata.items()):
            sc.pl.scatter(single_adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_ribo", ax=ax, show=False)
            ax.set_xlim(0, 40000)
            ax.set_ylim(0, 6000)
            ax.grid(True)
            ax.set_title(f"{adata_key} - pct_counts_ribo")
            ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1))
        plt.tight_layout()
        plt.show()
    else:
        # Handling single AnnData object
        print(f"Visualizations {title}:")

        # Violin plots
        sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"], jitter=0.4, multi_panel=True)
        
        # Adding gridlines to the mitochondrial and ribosomal plots
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.4, ax=axes[0], show=False)
        sc.pl.violin(adata, ["pct_counts_ribo"], jitter=0.4, ax=axes[1], show=False)

        for ax in axes:
            ax.grid(True)

        plt.tight_layout()
        plt.show()

        # Scatter plots
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", ax=ax[0], show=False)
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_ribo", ax=ax[1], show=False)

        for a in ax:
            a.set_xlim(0, 40000)
            a.set_ylim(0, 6000)
            a.grid(True)
            a.legend(loc='upper right', bbox_to_anchor=(1.2, 1))

        plt.tight_layout()
        plt.show()


def hash_demulitplex(adata, hashtag_input, number_of_noise_barcodes=None):
    """
    Function to demultiplex droplets using the Hashsolo package.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object.
    hashtag_input : str or list-like
        - If a str, interpret it as a prefix for hashtag genes in `var['gene_ids']`.
        - If a list or Index, interpret it as explicit gene_ids/var_names to use.
    number_of_noise_barcodes : int, optional
        Number of noise barcodes to pass to `hashsolo`.

    Returns
    -------
    AnnData
        The updated AnnData object with demultiplexing results.
    """
    import numpy as np
    import scanpy.external as sce

    # 1. Figure out if `hashtag_input` is a prefix (string) or list of genes.
    if isinstance(hashtag_input, str):
        # We assume it's a prefix in the 'gene_ids' column.
        hashtag_mask = adata.var['gene_ids'].str.startswith(hashtag_input)
    else:
        # Otherwise, assume it's a list/Index of actual gene names or gene_ids.
        # We can check both var.index (the var_names) and var['gene_ids']
        # in case your data is keyed either way:
        if hasattr(hashtag_input, "__iter__"):
            # Convert to a set for faster membership checking
            hashtag_input_set = set(hashtag_input)
            # Create a mask that is True if var.index OR var['gene_ids'] is in the set
            hashtag_mask = adata.var_names.isin(hashtag_input_set) | adata.var['gene_ids'].isin(hashtag_input_set)
        else:
            raise TypeError(
                "hashtag_input must be either a string prefix or a list of explicit gene names."
            )
    
    # 2. Check if we found any genes
    if hashtag_mask.sum() == 0:
        print(f"No matching hashtag genes found for input: {hashtag_input}")
        return adata

    # 3. Subset adata to those hashtag columns
    hto = adata[:, hashtag_mask].copy()

    # 4. Transfer counts from hto.X to hto.obs (one column per hashtag)
    # NOTE: For large data, consider using a sparse-friendly approach
    hto.obs[hto.var_names] = hto.X.toarray()

    # 5. Run HashSolo
    if number_of_noise_barcodes is not None:
        sce.pp.hashsolo(hto, hashtag_input, number_of_noise_barcodes=number_of_noise_barcodes)
    else:
        sce.pp.hashsolo(hto, hashtag_input)

    # 6. Print summary stats
    singlets = (hto.obs['most_likely_hypothesis'] == 1).sum()
    doublets = (hto.obs['most_likely_hypothesis'] == 2).sum()
    negatives = (hto.obs['most_likely_hypothesis'] == 0).sum()
    print(f'Number of predicted singlets: {singlets}')
    print(f'Number of predicted doublets: {doublets}')
    print(f'Number of predicted negatives: {negatives}')

    classifications = hto.obs['Classification'].value_counts()
    print('Number of cells for each hashtag:')
    for classification, count in classifications.items():
        print(f'{classification}: {count}')

    # 7. Map classifications back to original AnnData object
    #    Note: This mapping may need adjusting if your classification labels
    #    are different from what is in var['gene_ids'] or var_names.
    def map_to_var_index(label):
        # If the label is directly in adata.var['gene_ids'], get the var index (like 'K257')
        if label in adata.var['gene_ids'].values:
            return adata.var.index[adata.var['gene_ids'] == label][0]
        else:
            # Return the label as-is if not found
            return label

    adata.obs['Classification'] = hto.obs['Classification'].map(map_to_var_index)
    adata.obs['most_likely_hypothesis'] = hto.obs['most_likely_hypothesis']

    return adata


def plot_total_counts_vs_cells(
    adata, 
    bins=250, 
    zoom_x_range=(0, 3000), 
    zoom_y_range=(0, 600), 
    xlabel="Total counts", 
    ylabel="Number of cells", 
    figure_size=(10, 6), 
    inset_size=("30%", "30%"), 
    inset_location="upper right", 
    border_pad=2, 
    rotation=45
):
    """
    Plots the total counts vs the number of cells for a given anndata object, with an optional zoomed inset.
    
    Parameters:
        adata (AnnData): The annotated data matrix.
        bins (int): Number of bins for the histogram.
        zoom_x_range (tuple): Range for the x-axis in the zoomed inset (default: (0, 3000)).
        zoom_y_range (tuple): Range for the y-axis in the zoomed inset (default: (0, 600)).
        xlabel (str): Label for the x-axis (default: "Total counts").
        ylabel (str): Label for the y-axis (default: "Number of cells").
        figure_size (tuple): Size of the main figure (default: (10, 6)).
        inset_size (tuple): Size of the inset plot as a percentage of the main plot (default: ("30%", "30%")).
        inset_location (str): Location of the inset in the plot (default: "upper right").
        border_pad (int): Padding around the inset plot (default: 2).
        rotation (int): Rotation for the x-axis tick labels (default: 45).
    """
    # Main plot
    fig, ax = plt.subplots(figsize=figure_size)
    
    # Plot the histogram of total counts
    ax.hist(adata.obs["total_counts"], bins=bins, log=False, label="Total counts in bins")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    # Set x-ticks and rotate labels
    ax.set_xticks(range(0, int(adata.obs["total_counts"].max()) + 1000, 1000))
    plt.xticks(rotation=rotation)
    
    # Create inset plot for zoomed area
    axins = inset_axes(ax, width=inset_size[0], height=inset_size[1], loc=inset_location, borderpad=border_pad)
    axins.hist(adata.obs["total_counts"], bins=bins, log=False)
    
    # Set the zoom range in the inset plot
    x1, x2 = zoom_x_range
    y1, y2 = zoom_y_range
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    
    # Add a rectangle to highlight the zoomed area in the main plot
    rect = plt.Rectangle((x1, y1), x2-x1, y2-y1, edgecolor='r', facecolor='none', linestyle='--')
    ax.add_patch(rect)
    
    # Show the plot
    plt.show()

# Example usage:
# plot_total_counts_vs_cells(adata, bins=200, zoom_x_range=(0, 5000), zoom_y_range=(0, 1000))



def annotate_cellcycle_mouse(adata):
    # Get cell cycle genes from the regev lab data file. sort/split into appropriate lists
    cell_cycle_genes = [x.strip().lower().title() for x in open('../../regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x.lower().title() for x in cell_cycle_genes if x.lower().title() in adata.var_names]
    # Cell cycle scoring function from scanpy
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    # subset adata with only necessary genes
    adata_cc_genes = adata[:, cell_cycle_genes]
    

    adata.obs["cellcycle"] = adata_cc_genes.obs["phase"]
    adata.obs["S_score"] = adata_cc_genes.obs["S_score"]
    adata.obs["G2M_score"] = adata_cc_genes.obs["G2M_score"]
    
    return adata

# Function to generate a dictionary of lists of genes relevant to analysis which are present in the dataset
def find_genes(adata, marker_genes: dict) -> dict:

    #use supplied marker_genes dictionary to find genes of interest
    Genes_of_interest = []
    found_genes = {key: [] for key in marker_genes.keys()}

    # loop through gene lists and genes, see if it is present in the passed adata, if so append it to a list+dictionary 
    # so we have a dicitonary of lists of only genes in the dataset, as well as a list of all genes from the different types
    for gene_list_name, gene_list in marker_genes.items():
        for gene in gene_list:
            if gene in adata.var.index:
                print(f" ✓ {gene} is in adata.var")
                found_genes[gene_list_name].append(gene)
            else:
                print(f"{gene} is not in adata.var")

    return found_genes



def plot_umap_with_top_genes(adata, umap_color=None, n_top_genes=3):
    """
    Plots a square UMAP embedding with clusters colored by 'umap_color' on the left,
    and displays individual tables for each cluster stacked vertically on the right.
    Each table shows the top 'n_top_genes' genes for the cluster, with gene names in the top row,
    and their log fold changes, -10log(p-values), 'pts', and mean expression beneath.
    The first column contains the labels 'Gene', 'LFC', '-10logp', 'pts', and 'MeanExpr'.
    Gene names are colored red if the log fold change is positive and blue if negative.

    Parameters:
    - adata: AnnData object containing your single-cell data.
    - umap_color: str, the key in adata.obs to color the UMAP plot (e.g., 'leiden_0.5').
    - n_top_genes: int, number of top genes to display for each cluster.

    This function assumes that:
    - Differential expression analysis has been performed and the results are stored in adata.uns['rank_genes_groups'].
    - The name of the adata object is stored in adata.uns['name'].
    - 'pts' is stored in adata.uns['rank_genes_groups']['pts'] as a pandas DataFrame.
    """
    # Helper function to compute mean expression per group
    def grouped_obs_mean(adata, group_key, layer=None):
        """
        Computes the mean expression of genes for each group and the rest of the cells.

        Parameters:
        - adata: AnnData object containing your single-cell data.
        - group_key: str, the key in adata.obs to group the cells (e.g., 'leiden_0.5').
        - layer: str or None, the layer of adata to use for expression values.

        Returns:
        - mean_expr_df: DataFrame of mean expressions for each group.
        - mean_expr_rest_df: DataFrame of mean expressions for the rest of the cells.
        """
        if layer is not None:
            getX = lambda x: x.layers[layer]
        else:
            getX = lambda x: x.X

        groups = adata.obs[group_key].unique()
        mean_expr_df = pd.DataFrame(index=adata.var_names)
        mean_expr_rest_df = pd.DataFrame(index=adata.var_names)

        for group in groups:
            # Cells in the group
            idx_in_group = adata.obs[group_key] == group
            # Cells not in the group
            idx_rest = ~idx_in_group

            # Mean expression in the group
            X_in_group = getX(adata[idx_in_group])

            if isinstance(X_in_group, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)):
                mean_expr_in_group = X_in_group.mean(axis=0).A1  # Convert to 1D array
            else:
                mean_expr_in_group = np.asarray(X_in_group.mean(axis=0)).ravel()

            mean_expr_df[group] = mean_expr_in_group

            # Mean expression in the rest of the cells
            X_rest = getX(adata[idx_rest])

            if isinstance(X_rest, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)):
                mean_expr_rest = X_rest.mean(axis=0).A1  # Convert to 1D array
            else:
                mean_expr_rest = np.asarray(X_rest.mean(axis=0)).ravel()

            mean_expr_rest_df[group] = mean_expr_rest

        return mean_expr_df, mean_expr_rest_df

    # Retrieve adata name from adata.uns['name']
    adata_name = adata.uns.get('name', 'adata')

    # Access the rank genes groups data
    rank_genes_groups = adata.uns['rank_genes_groups']
    gene_names = rank_genes_groups['names']  # Structured array
    logfoldchanges = rank_genes_groups['logfoldchanges']
    pvals_adj = rank_genes_groups['pvals_adj']

    # Access 'pts' as a DataFrame
    pts = adata.uns['rank_genes_groups']['pts']
    pts_rest = adata.uns['rank_genes_groups']['pts_rest']

    # Compute mean expression per group
    mean_expr_df, mean_expr_rest_df = grouped_obs_mean(adata, group_key=umap_color)

    # Get group names from rank_genes_groups
    group_names = gene_names.dtype.names

    num_groups = len(group_names)

    # Calculate figure height based on number of groups
    fig_height = max(6, num_groups * 2)  # Ensure minimum height
    fig_width = 14  # Wider figure to accommodate both UMAP and tables

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(nrows=num_groups, ncols=2, width_ratios=[1, 1.5])

    # Create the UMAP plot axis occupying the left side
    ax_umap = fig.add_subplot(gs[:, 0])

    # Generate the UMAP plot on the ax_umap
    sc.pl.umap(adata, color=umap_color, legend_loc="on data", ax=ax_umap, show=False)

    # Ensure the UMAP plot is square
    ax_umap.set_aspect('equal')
    

    # Set the plot title to include adata_name and umap_color
    ax_umap.set_title(f"UMAP of {adata_name} colored by '{umap_color}'")

    # Add definitions text under the UMAP plot
    definitions_text = (
        r"$\mathbf{Genes\ ordered\ by\ Z-score.}$""\n\n"
        r"$\mathbf{LFC}$: Fold change for each gene for each group.""\n\n"
    r"$\mathbf{-log10p}$: -log10 Benjamini-Hochberg corrected p-values.""\n\n"
    r"$\mathbf{MeanExpr}$: Mean expression of the gene in the group.""\n\n"
    r"$\mathbf{MeanExprRest}$: Mean expression of the gene in the rest of the groups""\n""(excluding current group).""\n\n"
    r"$\mathbf{pts}$: Fraction of cells expressing the genes for each group.""\n\n"
    r"$\mathbf{pts\_rest}$: Fraction of cells from the union of the rest of each group""\n"" expressing the genes."
)
    # Create an axis for the definitions text
    ax_definitions = fig.add_subplot(gs[0, 0])
    # Position the text under the UMAP plot
    ax_definitions.axis('off')  # Hide the axis
    ax_definitions.text(
        -0.1, -0.1, definitions_text,
        ha='left', va='top',
        fontsize=10,
        transform=ax_umap.transAxes,
        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=5)
    )

    # Prepare data for tables
    tables_data = {}
    for group_name in group_names:
        top_genes = gene_names[group_name][:n_top_genes]
        top_logfoldchanges = logfoldchanges[group_name][:n_top_genes]
        top_pvals_adj = pvals_adj[group_name][:n_top_genes]
        top_neg_log_pvals = -np.log10(top_pvals_adj)

        # Retrieve 'pts' values for the top genes and current group
        top_pts = pts.loc[top_genes, group_name].values
        top_pts_rest = pts_rest.loc[top_genes, group_name].values

        # Retrieve mean expression values for the top genes in the current group
        mean_expr = mean_expr_df.loc[top_genes, group_name].values
        mean_expr_rest = mean_expr_rest_df.loc[top_genes, group_name].values

        # Store data in a list of lists
        table = [
            ['Gene'] + list(top_genes),
            ['LFC'] + ["{:.2f}".format(lfc) for lfc in top_logfoldchanges],
            ['-10logp'] + ["{:.2f}".format(pval) for pval in top_neg_log_pvals],
            ['MeanExpr'] + ["{:.2f}".format(me) for me in mean_expr],
            ['MeanExpr_rest'] + ["{:.2f}".format(me) for me in mean_expr_rest],
            ['pts'] + ["{:.2f}".format(pt) for pt in top_pts],
            ['pts_rest'] + ["{:.2f}".format(pt) for pt in top_pts_rest]            
        ]

        # Store the colors for gene names
        gene_colors = ['black'] + ['red' if lfc > 0 else 'blue' for lfc in top_logfoldchanges]

        tables_data[group_name] = (table, gene_colors)

    # Now, create a table for each group, stacked vertically on the right
    for i, group_name in enumerate(group_names):
        table_data, gene_colors = tables_data[group_name]

        # Create a new axis for the table on the right
        ax_subtable = fig.add_subplot(gs[i, 1])
        ax_subtable.axis('off')

        # Create the table
        the_table = ax_subtable.table(cellText=table_data,
                                      cellLoc='center',
                                      loc='center')

        # Adjust the table properties
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 2)  # Adjust table scale if necessary

        # Set cell colors and properties
        for (row, col), cell in the_table.get_celld().items():
            cell.set_linewidth(0.5)
            if row == 0 and col > 0:
                # Gene names in the top row
                cell._text.set_color(gene_colors[col])
                cell.set_text_props(weight='bold')
                #truncate gene names if they are too long
                cell._text.set_text(cell._text.get_text()[:7])
            if col == 0 and row > 0:
                # Labels in the first column
                cell.set_text_props(weight='bold')
                cell._text.width = 0.2
            if row == 0 or col == 0:
                cell.set_facecolor('#f1f1f1')  # Light grey background for headers

        # Set the table title
        ax_subtable.set_title(f"Cluster {group_name} vs Rest", fontsize=10, pad=25)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(hspace=1)
    

    # Display the plot
    plt.show()

def plot_umap_with_gene_list(adata, umap_color, n_top_genes=3):
    """
    Plots a square UMAP embedding with clusters colored by 'umap_color' and displays
    a list of top 'n_top_genes' for each cluster in a horizontal row. Cluster names are 
    displayed in a single row, with corresponding genes listed below each name and 
    colored by up-regulation (red) or down-regulation (blue) based on log fold change.

    Parameters:
    - adata: AnnData object containing your single-cell data.
    - umap_color: str, the key in adata.obs to color the UMAP plot (e.g., 'leiden_0.5').
    - n_top_genes: int, number of top genes to display for each cluster.

    This function assumes that:
    - Differential expression analysis has been performed and the results are stored in adata.uns['rank_genes_groups'].
    - The name of the adata object is stored in adata.uns['name'].
    """
    # Retrieve adata name from adata.uns['name']
    adata_name = adata.uns.get('name', 'adata')

    # Ensure using log_norm prior to running rank_genes_groups
    adata.X = adata.layers['log_norm'] if 'log_norm' in adata.layers.keys() else print("No log_norm layer found in adata.layers")

    # Run rank_genes with umap_color parameter as groupby to ensure correct gene ranking if previous run with different groupby
    sc.tl.rank_genes_groups(adata, groupby=umap_color, method='wilcoxon', rankby_abs = True)

    # Access the rank genes groups data
    rank_genes_groups = adata.uns['rank_genes_groups']
    gene_names = rank_genes_groups['names']  # Structured array
    logfoldchanges = rank_genes_groups['logfoldchanges']

    # Get group names from rank_genes_groups
    group_names = gene_names.dtype.names

    # Set up figure
    fig, ax_umap = plt.subplots(figsize=(8, 8))

    # Generate the UMAP plot
    sc.pl.umap(adata, color=umap_color, legend_loc="on data", ax=ax_umap, show=False)
    ax_umap.set_aspect('equal')
    ax_umap.set_title(f"UMAP of {adata_name} colored by '{umap_color}'")

    # Define starting position for the cluster names (horizontal row)
    y_pos = 0.05
    x_start = 0.1
    x_spacing = 0.15  # Horizontal spacing between clusters
    gene_y_spacing = 0.03  # Vertical spacing between genes

    # Loop through clusters and add cluster name and gene lists below
    for i, group_name in enumerate(group_names):
        # Position the cluster name at the top of each column
        x_pos = x_start + i * x_spacing
        fig.text(x_pos, y_pos, f"Cluster {group_name}", ha='center', va='top', fontsize=10, weight='bold')

        # List top genes below each cluster name
        for j, (gene, lfc) in enumerate(zip(gene_names[group_name][:n_top_genes], logfoldchanges[group_name][:n_top_genes])):
            color = 'red' if lfc > 0 else 'blue'
            gene_y_pos = y_pos - (j + 1) * gene_y_spacing
            fig.text(x_pos, gene_y_pos, gene, ha='center', va='top', fontsize=10, color=color)

    # Adjust layout and show plot
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.3)  # Add space at the bottom for gene text
    plt.show()

def gene_expression_percentage(adata, gene_id, category, logfc_threshold=0.5):
    """
    Prints the percentage of cells in which the specified gene is expressed above 0 and above a given logFC threshold,
    both in total and within each group of the specified category.

    Parameters:
    adata (AnnData): The AnnData object containing single-cell RNA sequencing data.
    gene_id (str): The ID of the gene to analyze.
    category (str): The name of the category in adata.obs to group the cells.
    logfc_threshold (float): The log fold-change threshold for considering a gene as expressed. Default is 0.5.
    """
    if gene_id not in adata.var_names:
        print(f"Gene {gene_id} not found in the dataset.")
        return

    if category not in adata.obs:
        print(f"Category {category} not found in adata.obs.")
        return

    # Get the expression values for the gene across all cells
    gene_expression = adata[:, gene_id].X

    # Flatten the array to ensure it's a 1D array
    gene_expression = gene_expression.A.flatten() if hasattr(gene_expression, "A") else gene_expression

    # Calculate the percentage of cells with expression > 0 in total
    total_cells = gene_expression.shape[0]
    expressed_cells = (gene_expression > 0).sum()
    percentage_expressed = (expressed_cells / total_cells) * 100

    # Calculate the percentage of cells with expression > logfc_threshold in total
    highly_expressed_cells = (gene_expression > logfc_threshold).sum()
    percentage_highly_expressed = (highly_expressed_cells / total_cells) * 100

    print(f"Gene {gene_id}:")
    print(f"Total percentage of cells with expression > 0: {percentage_expressed:.2f}%")
    print(f"Total percentage of cells with expression > {logfc_threshold} logFC: {percentage_highly_expressed:.2f}%")

    # Get the unique groups in the specified category
    groups = adata.obs[category].unique()

    # Calculate the percentages within each group
    for group in groups:
        group_mask = adata.obs[category] == group
        group_gene_expression = gene_expression[group_mask]

        group_total_cells = group_gene_expression.shape[0]
        group_expressed_cells = (group_gene_expression > 0).sum()
        group_percentage_expressed = (group_expressed_cells / group_total_cells) * 100

        group_highly_expressed_cells = (group_gene_expression > logfc_threshold).sum()
        group_percentage_highly_expressed = (group_highly_expressed_cells / group_total_cells) * 100

        print(f"\nGroup '{group}':")
        print(f"Percentage of cells with expression > 0: {group_percentage_expressed:.2f}%")
        print(f"Percentage of cells with expression > {logfc_threshold} logFC: {group_percentage_highly_expressed:.2f}%")

# Example usage
# Assuming 'adata' is your pre-processed AnnData object, 'gene_id' is the gene you want to analyze,
# and 'category' is the annotation in adata.obs (e.g., 'leiden' or any other category)
# gene_expression_percentage_by_category(adata, 'gene_id', 'category')



def plot_proportion_heatmap_total(df, group_col, feature_col, adata, umap_color):
    # Calculate the total number of cells for each genotype
    total_genotype_cells = df.groupby(group_col).size()
    
    # Group by group_col and feature_col, then calculate size
    count_df = df.groupby([group_col, feature_col]).size().unstack(fill_value=0)
    
    # Normalize by the total number of cells for each genotype
    proportion_df = count_df.div(total_genotype_cells, axis=0) * 100
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Plotting the heatmap
    sns.heatmap(proportion_df.T, annot=True, cmap='crest', fmt=".2f", linewidths=.5, ax=ax1)
    ax1.set_title(f'Proportion of {group_col} by {feature_col}')

    # Plotting the UMAP
    sc.pl.umap(adata, color=umap_color, legend_loc="on data", ax=ax2, show=False)
    ax2.set_title('UMAP Plot')

    plt.tight_layout()
    plt.show()

def plot_proportion_heatmap(df, group_col, feature_col, adata, umap_color):
    # Group by group_col and feature_col, then calculate size
    count_df = df.groupby([group_col, feature_col]).size().unstack(fill_value=0)
    
    # Normalize to get the percentage of each feature type within each group
    proportion_df = count_df.div(count_df.sum(axis=0), axis=1) * 100

    # Debugging: Print the proportion DataFrame to verify
    print(proportion_df)
    print(proportion_df.sum(axis=1))  # This should show 100 for each row

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Plotting the heatmap
    sns.heatmap(proportion_df, annot=True, cmap='crest', fmt=".2f", linewidths=.5, ax=ax1)
    ax1.set_title(f'Proportion of {feature_col} by {group_col}')

    # Plotting the UMAP
    sc.pl.umap(adata, color=umap_color, legend_loc="on data", ax=ax2, show=False)
    ax2.set_title('UMAP Plot')

    plt.tight_layout()
    plt.show()

# def plot_volcano(adata, group='LPS', logfc_threshold=0.3, pval_threshold=0.05):
# ## based on the code from the following link:
# ## https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html

#     # Extract log fold changes and adjusted p-values for the specified group
#     logfc = adata.uns['rank_genes_groups']['logfoldchanges'][group]
#     pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][group]
    
#     # Convert adjusted p-values to -log10 scale
#     neg_log_pvals = -np.log10(pvals_adj)
    
#     # Plot all genes
#     plt.figure(figsize=(10, 7))
#     plt.scatter(logfc, neg_log_pvals, s=1, label="Not significant", color='grey')

#     # Identify up- and down-regulated genes based on thresholds
#     down = (logfc <= -logfc_threshold) & (pvals_adj <= pval_threshold)
#     up = (logfc >= logfc_threshold) & (pvals_adj <= pval_threshold)

#     # Plot down-regulated genes
#     plt.scatter(logfc[down], neg_log_pvals[down], s=3, label="Down-regulated", color="blue")

#     # Plot up-regulated genes
#     plt.scatter(logfc[up], neg_log_pvals[up], s=3, label="Up-regulated", color="red")

#     # Annotate significant points
#     texts = []
#     for i in np.where(up)[0]:
#         texts.append(plt.text(logfc[i], neg_log_pvals[i], adata.uns['rank_genes_groups']['names'][group][i]))


#     adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.1))

#     # Add labels and lines for thresholds
#     plt.xlabel("Log Fold Change (logFC)")
#     plt.ylabel("-Log10 Adjusted P-Value (-logFDR)")
#     plt.axvline(-logfc_threshold, color="grey", linestyle="--")
#     plt.axvline(logfc_threshold, color="grey", linestyle="--")
#     plt.axhline(-np.log10(pval_threshold), color="grey", linestyle="--")
#     plt.xlim(-10,10)
#     plt.ylim(0, 20)
#     plt.legend()
#     plt.title(f"Volcano Plot for {group}")
#     plt.show()


def plot_volcano(adata, group='LPS', logfc_threshold=0.3, pval_threshold=0.05, top_n=20):
    # Get Name of data
    data_name = adata.uns['name']

    # Get reference group
    ref_group = adata.uns['rank_genes_groups']['params']['reference']

    # Create a DataFrame with necessary data
    df = pd.DataFrame({
        'logfc': adata.uns['rank_genes_groups']['logfoldchanges'][group],
        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][group],
        'gene_names': adata.uns['rank_genes_groups']['names'][group]
    })

    # Replace zero p-values with the smallest positive float
    df['pvals_adj'] = df['pvals_adj'].replace(0, np.nextafter(0, 1))

    # Calculate -log10 adjusted p-values
    df['neg_log_pvals'] = -np.log10(df['pvals_adj'])

    # Check neg_log_pvals for inf and nan values
    neg_log_pvals = np.array(df['neg_log_pvals'])
    print(f"Number of inf values in neg_log_pvals: {np.isinf(neg_log_pvals).sum()}")
    print(f"Number of nan values in neg_log_pvals: {np.isnan(neg_log_pvals).sum()}")

    # Replace inf and nan values with 0
    df['neg_log_pvals'] = df['neg_log_pvals'].replace([np.inf, -np.inf], 0)

    # Assign significance categories
    df['significance'] = 'Not significant'
    df.loc[(df['logfc'] >= logfc_threshold) & (df['pvals_adj'] <= pval_threshold), 'significance'] = 'Up-regulated'
    df.loc[(df['logfc'] <= -logfc_threshold) & (df['pvals_adj'] <= pval_threshold), 'significance'] = 'Down-regulated'

    # Define color mapping
    colors = {'Not significant': 'grey', 'Up-regulated': 'red', 'Down-regulated': 'blue'}

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot each category
    for key, group_df in df.groupby('significance'):
        ax.scatter(group_df['logfc'], group_df['neg_log_pvals'], s=3, color=colors[key], label=key, alpha=0.7)

    # Collect text annotations
    down_texts = []
    up_texts = []

    # Annotate top N up-regulated genes
    up_genes = df[df['significance'] == 'Up-regulated']
    top_up_genes = up_genes.nsmallest(top_n, 'pvals_adj')

    for _, row in top_up_genes.iterrows():
        up_texts.append(
            ax.text(row['logfc'], row['neg_log_pvals'], row['gene_names'], fontsize=12)
        )

    # Annotate top N down-regulated genes
    down_genes = df[df['significance'] == 'Down-regulated']
    top_down_genes = down_genes.nsmallest(top_n, 'pvals_adj')

    for _, row in top_down_genes.iterrows():
        down_texts.append(
            ax.text(row['logfc'], row['neg_log_pvals'], row['gene_names'], fontsize=12)
        )



    # Create a dictionary of top N up- and down-regulated genes
    annotated_genes = {
        'Up-regulated': top_up_genes['gene_names'].tolist(),
        'Down-regulated': top_down_genes['gene_names'].tolist()
    }
    # Add threshold lines
    neg_vline = ax.axvline(-logfc_threshold, color="grey", linestyle="--")
    pos_vline = ax.axvline(logfc_threshold, color="grey", linestyle="--")
    hline = ax.axhline(-np.log10(pval_threshold), color="grey", linestyle="--")



    # Set labels and title
    ax.set_xlabel("Log Fold Change (logFC)")
    ax.set_ylabel("-Log10 Adjusted P-Value (-logFDR)")
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, df['neg_log_pvals'].max() + 1)
    ax.legend()
    title_obj = ax.set_title(f"Volcano Plot for {group} vs {ref_group} in {data_name}")

    mpl_objs = [neg_vline, pos_vline, hline, title_obj]

    adjust_text(
        up_texts,
        ax=ax,
        objects=mpl_objs,
        ensure_inside_axes=True,
        explode_radius=6,
        pull_threshold=3,
        force_text=[5, 5],
        force_static=[2, 1],
        force_pull=[0.1, 0.1],
        only_move={'static': 'x+', 'text': 'x+'},
        arrowprops=dict(arrowstyle='->', color='black', shrinkA=8),
        max_move=25,
        min_arrow_len=1,
        expand_axes=True
    )

    adjust_text(
        down_texts,
        ax=ax,
        objects=mpl_objs,
        ensure_inside_axes=True,
        explode_radius=6,
        pull_threshold=3,
        force_text=[5, 5],
        force_static=[2, 1],
        force_pull=[0.1, 0.1],
        only_move={'static': 'x-', 'text': 'x-'},
        arrowprops=dict(arrowstyle='->', color='black', shrinkA=8),
        max_move=25,
        min_arrow_len=1,
        expand_axes=True
    )
    # Improve layout
    plt.subplots_adjust(top=0.85)
    plt.tight_layout()
    plt.show()

    return annotated_genes

def plot_volcano_from_df(
    df,
    logfc_col='log2FoldChange',
    pval_col='padj',
    gene_col=None,
    logfc_threshold=0.3,
    pval_threshold=0.05,
    top_n=20,
    title='Volcano Plot',
    figsize=(10, 7)
):
    """
    Plots a volcano plot from a DataFrame containing differential expression results.

    Parameters:
    - df: pandas DataFrame containing the differential expression results.
    - logfc_col: str, name of the column containing log2 fold changes.
    - pval_col: str, name of the column containing adjusted p-values.
    - gene_col: str or None, name of the column containing gene names. If None, the DataFrame's index is used.
    - logfc_threshold: float, threshold for log2 fold change to consider as significant.
    - pval_threshold: float, threshold for adjusted p-value to consider as significant.
    - top_n: int, number of top genes to annotate for upregulated and downregulated genes.
    - title: str, title of the plot.
    - figsize: tuple, size of the figure.
    """
    # Create a copy of the DataFrame to avoid modifying the original
    df = df.copy()

    # If gene_col is specified, use it; otherwise, use the index as gene names
    if gene_col:
        df['gene_names'] = df[gene_col]
    else:
        df['gene_names'] = df.index

    # Replace zero or negative p-values with the smallest positive float
    df[pval_col] = df[pval_col].replace(0, np.nextafter(0, 1e-10))

    # Calculate -log10 adjusted p-values
    df['neg_log_pvals'] = -np.log10(df[pval_col])

    # Check for infinite or NaN values
    neg_log_pvals = np.array(df['neg_log_pvals'])
    num_inf = np.isinf(neg_log_pvals).sum()
    num_nan = np.isnan(neg_log_pvals).sum()
    print(f"Number of inf values in neg_log_pvals: {num_inf}")
    print(f"Number of nan values in neg_log_pvals: {num_nan}")

    # Replace inf and NaN values with zero
    df['neg_log_pvals'] = df['neg_log_pvals'].replace([np.inf, -np.inf], 0)
    df['neg_log_pvals'] = df['neg_log_pvals'].fillna(0)

    # Assign significance categories
    df['significance'] = 'Not significant'
    df.loc[
        (df[logfc_col] >= logfc_threshold) & (df[pval_col] <= pval_threshold),
        'significance'
    ] = 'Up-regulated'
    df.loc[
        (df[logfc_col] <= -logfc_threshold) & (df[pval_col] <= pval_threshold),
        'significance'
    ] = 'Down-regulated'

    # Define color mapping
    colors = {'Not significant': 'grey', 'Up-regulated': 'red', 'Down-regulated': 'blue'}

    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)

    # Plot each category
    for key, group_df in df.groupby('significance'):
        ax.scatter(
            group_df[logfc_col],
            group_df['neg_log_pvals'],
            s=10,
            color=colors[key],
            label=key,
            alpha=0.7
        )

    # Collect text annotations
    up_texts = []
    down_texts = []

    # Annotate top N upregulated genes
    up_genes = df[df['significance'] == 'Up-regulated']
    top_up_genes = up_genes.nsmallest(top_n, pval_col)

    for _, row in top_up_genes.iterrows():
        up_texts.append(
            ax.text(row[logfc_col], row['neg_log_pvals'], row['gene_names'], fontsize=12)
        )

    # Annotate top N downregulated genes
    down_genes = df[df['significance'] == 'Down-regulated']
    top_down_genes = down_genes.nsmallest(top_n, pval_col)

    for _, row in top_down_genes.iterrows():
        down_texts.append(
            ax.text(row[logfc_col], row['neg_log_pvals'], row['gene_names'], fontsize=12)
        )

    # Add threshold lines
    neg_vline = ax.axvline(-logfc_threshold, color="grey", linestyle="--")
    pos_vline = ax.axvline(logfc_threshold, color="grey", linestyle="--")
    hline = ax.axhline(-np.log10(pval_threshold), color="grey", linestyle="--")

    # Set labels and title
    ax.set_xlabel("Log2 Fold Change (log2FC)")
    ax.set_ylabel("-Log10 Adjusted P-Value (-log10 padj)")
    ax.set_xlim(df[logfc_col].min() - 1, df[logfc_col].max() + 1)
    ax.set_ylim(0, df['neg_log_pvals'].max() + 1)
    ax.legend()
    title_obj = ax.set_title(title)

    mpl_objs = [neg_vline, pos_vline, hline, title_obj]

    # Adjust text annotations to prevent overlap
    adjust_text(
        up_texts + down_texts,
        ax=ax,
        objects=mpl_objs,
        arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
        expand_text=(1.2, 1.2),
        expand_points=(1.2, 1.2),
        expand_objects=(1.2, 1.2),
        force_text=(0.5, 0.5),
        force_points=(0.2, 0.2),
        lim=1000,
    )

    # Improve layout
    plt.tight_layout()
    plt.show()

    # Create a dictionary of top N up- and down-regulated genes
    annotated_genes = {
        'Up-regulated': top_up_genes['gene_names'].tolist(),
        'Down-regulated': top_down_genes['gene_names'].tolist()
    }

    return annotated_genes

def get_significant_genes(adata, group='', logfc_threshold=0.3, pval_threshold=0.05, top_n=None):
    # Extract log fold changes and adjusted p-values for the specified group
    logfc = adata.uns['rank_genes_groups']['logfoldchanges'][group]
    pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][group]
    
    # Identify up- and down-regulated genes based on thresholds
    down = (logfc <= -logfc_threshold) & (pvals_adj <= pval_threshold)
    up = (logfc >= logfc_threshold) & (pvals_adj <= pval_threshold)
    
    # Get the names of the significant genes
    down_genes = adata.uns['rank_genes_groups']['names'][group][down]
    up_genes = adata.uns['rank_genes_groups']['names'][group][up]

    # sort the genes by pvals_adj
    down_genes = [x for _, x in sorted(zip(adata.uns['rank_genes_groups']['pvals_adj'][group][down], down_genes))]
    up_genes = [x for _, x in sorted(zip(adata.uns['rank_genes_groups']['pvals_adj'][group][up], up_genes))]

    # If top_n is provided, take only the first top_n entries
    if top_n is not None:
        down_genes = down_genes[:top_n]
        up_genes = up_genes[:top_n]
    
    
    return down_genes, up_genes

def venn_gene(adata1, adata2, group, logfc_threshold=0.3, pval_threshold=0.05, top_n=None):
    
    
    # Get the significant genes for each group

    if top_n is not None:
        down_genes1, up_genes1 = get_significant_genes(adata1, group, logfc_threshold, pval_threshold, top_n)
        down_genes2, up_genes2 = get_significant_genes(adata2, group, logfc_threshold, pval_threshold, top_n)
    else:
        down_genes1, up_genes1 = get_significant_genes(adata1, group, logfc_threshold, pval_threshold)
        down_genes2, up_genes2 = get_significant_genes(adata2, group, logfc_threshold, pval_threshold)
    
    # Create the Venn diagram
    plt.figure(figsize=(12, 8))
    venn2([set(up_genes1), set(up_genes2)])
    plt.title("Up-regulated Genes")
    plt.show()

    plt.figure(figsize=(12, 8))
    venn2([set(down_genes1), set(down_genes2)], set_labels=[adata1, adata2])
    plt.title("Down-regulated Genes")
    plt.show()
    

def plot_dual_contrast(adata, gene, constant_contrast, contrast_variable, level1, level2):
    """
    Plot a bar plot comparing two contrast groups (e.g. treatment: VEHICLE vs LPS)
    while holding other experimental conditions constant (e.g. disease and apoe).

    For each contrast group:
      - Cells are filtered by the constant conditions and the contrast variable is set 
        to the specified level.
      - The function groups cells by 'Classification' (i.e. sample) and computes the 
        mean normalized count for the gene of interest.
      - A bar is drawn with height equal to the overall mean (across samples), with 
        individual sample means overlaid as dots.
    
    In addition, the function performs a two-sample t-test on the aggregated sample-level 
    means (comparing the two contrast groups). The p-value is BH-corrected, and both the 
    corrected p-value and significance stars are printed and annotated on the plot.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object with expression data and observations.
    gene : str
        The gene of interest.
    constant_contrast : dict
        Dictionary of conditions that remain constant (e.g. {'disease': 'E4', 'apoe': 'WT'}).
    contrast_variable : str
        The variable being contrasted (e.g. 'treatment').
    level1 : str
        First level for the contrast variable (e.g. 'VEHICLE').
    level2 : str
        Second level for the contrast variable (e.g. 'LPS').
    
    Returns:
    --------
    None
        Displays the annotated bar plot.
    """
    
    def get_sample_means(adata_subset, gene):
        # Use log-normalized counts (assumed in adata.layers['log_norm'])
        adata_subset.X = adata_subset.layers['log_norm'].copy()
        
        if gene not in adata_subset.var_names:
            raise ValueError(f"Gene '{gene}' not found in adata.var_names.")
        
        # Extract gene expression values as a 1D array
        gene_expr = adata_subset[:, gene].X
        if hasattr(gene_expr, "toarray"):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = np.array(gene_expr).flatten()
        
        # Create DataFrame with expression and sample identity (Classification)
        df = pd.DataFrame({
            'normalized_counts': gene_expr,
            'Classification': adata_subset.obs['Classification'].values
        })
        # Compute mean expression per sample
        sample_means = df.groupby('Classification')['normalized_counts'].mean()
        return sample_means

    def subset_for_level(level):
        mask = np.ones(adata.n_obs, dtype=bool)
        for key, value in constant_contrast.items():
            mask &= (adata.obs[key] == value)
        mask &= (adata.obs[contrast_variable] == level)
        return adata[mask].copy()
    
    # Subset the data for the two contrast levels
    adata_level1 = subset_for_level(level1)
    adata_level2 = subset_for_level(level2)
    
    # Compute sample-level means
    sample_means1 = get_sample_means(adata_level1, gene)
    sample_means2 = get_sample_means(adata_level2, gene)
    
    # Compute overall means for plotting
    overall_mean1 = sample_means1.mean() if not sample_means1.empty else np.nan
    overall_mean2 = sample_means2.mean() if not sample_means2.empty else np.nan
    
    # Prepare bar plot
    fig, ax = plt.subplots(figsize=(8, 6))
    conditions = [level1, level2]
    overall_means = [overall_mean1, overall_mean2]
    bar_positions = np.arange(len(conditions))
    
    bars = ax.bar(bar_positions, overall_means, color=['lightblue', 'lightgreen'],
                  alpha=0.6, width=0.6, label='Overall Mean Expression')
    
    # Overlay individual sample means as dots with jitter for clarity
    jitter_std = 0.08
    for i, sample_means in enumerate([sample_means1, sample_means2]):
        if not sample_means.empty:
            x_positions = np.random.normal(loc=i, scale=jitter_std, size=len(sample_means))
            ax.scatter(x_positions, sample_means.values, color='black', alpha=0.8,
                       zorder=10, label='Sample Expression' if i == 0 else "")
    
    # ----------------- Two-sample t-test on Sample-Level Means -----------------
    # Check if both groups have at least two samples
    n1 = sample_means1.shape[0]
    n2 = sample_means2.shape[0]
    if n1 < 2 or n2 < 2:
        print("Not enough samples in one or both groups to perform a reliable t-test.")
        p_corrected = np.nan
    else:
        # Perform two-sample t-test (unequal variance)
        t_stat, p_val = ttest_ind(sample_means1.values, sample_means2.values, equal_var=False)
        # Apply BH correction (trivial here with one test)
        reject, p_adj, _, _ = multipletests([p_val], method='fdr_bh')
        p_corrected = p_adj[0]
    
    # Determine significance stars based on corrected p-value
    if not np.isnan(p_corrected):
        if p_corrected < 0.001:
            significance = '***'
        elif p_corrected < 0.01:
            significance = '**'
        elif p_corrected < 0.05:
            significance = '*'
        else:
            significance = 'ns'
    else:
        significance = 'NA'
    
    # Print the p-value and significance
    print(f"BH-corrected p-value for {contrast_variable}: {p_corrected:.3g} {significance}")
    
    # Annotate the plot with both the p-value and significance stars
    x_center = np.mean(bar_positions)
    y_max = max(overall_means) * 1.1
    annotation = f"{significance}\np = {p_corrected:.3g}"
    ax.text(x_center, y_max, annotation, ha='center', va='bottom', fontsize=14)
    ax.set_ylim(0, y_max * 1.2)
    
    ax.set_xlabel("Contrast Group")
    ax.set_ylabel(f"Normalized counts for {gene}")
    ax.set_title(f"Expression of {gene}\nConstant: {constant_contrast}, Contrast: {contrast_variable}")
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(conditions)
    ax.legend()
    plt.tight_layout()
    plt.show()


# def plot_volcano2(
#     adata,
#     group,
#     key="rank_genes_groups",
#     title=None,
#     adjusted_pvals=True,
#     show=True,
#     logfc_threshold=0.2,
#     filter_kwargs=None,
#     **kwargs
# ):
#     """
#     Plots the results of :func:`scanpy.tl.rank_genes_groups` in the form of a volcano plot.

#     Parameters
#     ----------
#     adata
#         Annotated data matrix.
#     group
#         Which group (as in :func:`scanpy.tl.rank_genes_groups`’s groupby argument) to return
#         results from. Can be a list. All groups are returned if `groups` is None.
#     key
#         Key differential expression groups were stored under.
#     title
#         Title of the resulting plot.
#     adjusted_pvals
#         Use adjusted p-values instead of raw p-values.
#     show
#         Whether to show the plot or return it.
#     logfc_threshold
#         Threshold for log fold change to determine vertical lines and coloring of points.
#     filter_kwargs
#         Keyword arguments to pass into :func:`scanpy.get.rank_genes_groups_df`.
#     kwargs
#         Keyword arguments to pass into :func:`seaborn.scatterplot`.

#     Returns
#     -------
#     If `show = False` returns the current axes. If `show = True` returns nothing.
#     """
#     if filter_kwargs is None:
#         filter_kwargs = {}
#     if title is None:
#         title = f"{key} {group}"

#     # Pull dataframe from adata object, and select columns of interest
#     de_df = sc.get.rank_genes_groups_df(adata, group=group, key=key, **filter_kwargs)
#     logfold = de_df["logfoldchanges"]
#     pvals = de_df["pvals_adj" if adjusted_pvals else "pvals"]
#     neg_log_pvals = np.negative(np.log10(pvals))

#     # Determine color based on thresholds
#     colors = np.where(
#         (logfold > logfc_threshold) & (neg_log_pvals > -np.log10(0.05)), "red",
#         np.where(
#             (logfold < -logfc_threshold) & (neg_log_pvals > -np.log10(0.05)), "blue", "lightgray"
#         )
#     )

#     # Plot logfold and -log pvals
#     ax = sns.scatterplot(
#         x=logfold,
#         y=neg_log_pvals,
#         hue=colors,  # Use the colors array for coloring
#         palette={"red": "red", "blue": "blue", "lightgray": "lightgray"},
#         legend=None,
#         size=0.1,
#         **kwargs
#     )

#     # Add vertical lines at logfc_threshold and -logfc_threshold
#     ax.axvline(logfc_threshold, color="gray", linestyle="--")
#     ax.axvline(-logfc_threshold, color="gray", linestyle="--")

#     # Add significance line at p = 0.05 and set title and axis labels
#     ax.axhline(-np.log10(0.05), 0, 1, color="lightgray", zorder=-10)
#     ax.set(ylabel="-log10(pval)", xlabel="logfoldchange", title=title)
#     ax.set_xlim(-10, 10)
#     ax.set_ylim(0, 20)

#     if show:
#         plt.show()
#     else:
#         return ax