# This file contains the functions to be used in the notebooks for the project

import re
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.sparse
from adjustText import adjust_text
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
    Perform quality control filtering on an AnnData object.

    Parameters:
    ----------
    adata : AnnData
        The input AnnData object containing single-cell RNA-seq data.
    mt_threshold : float, optional
        Maximum allowed percentage of mitochondrial gene counts per cell (default: 5%).
    ribo_threshold : float, optional
        Minimum required percentage of ribosomal gene counts per cell (default: 5%).
    min_counts : int, optional
        Minimum total UMI counts per cell required to retain the cell (default: 500).
    min_genes : int, optional
        Minimum number of genes detected per cell required to retain the cell (default: 200).
    min_cells : int, optional
        Minimum number of cells a gene must be expressed in to retain the gene (default: 3).

    Returns:
    -------
    AnnData
        Filtered AnnData object after applying all quality control thresholds.
    """

    # Annotate mitochondrial and ribosomal genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith("Rps") | adata.var_names.str.startswith("Rpl")

    print(f"Initial AnnData has {adata.n_obs} cells and {adata.n_vars} genes")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    # Identify cells failing each criterion
    cells_mt = adata.obs["pct_counts_mt"] > mt_threshold
    cells_ribo = adata.obs["pct_counts_ribo"] < ribo_threshold
    cells_counts = adata.obs["total_counts"] < min_counts
    cells_genes = adata.obs["n_genes_by_counts"] < min_genes

    # Print number of cells failing each individual filter
    print(f"{cells_mt.sum()} cells have >{mt_threshold}% mitochondrial genes.")
    print(f"{cells_ribo.sum()} cells have <{ribo_threshold}% ribosomal genes.")
    print(f"{cells_counts.sum()} cells have <{min_counts} total counts.")
    print(f"{cells_genes.sum()} cells have <{min_genes} expressed genes.")

    # Combine all filters into a single mask (remove if failing any)
    combined_failed_cells = cells_mt | cells_ribo | cells_counts | cells_genes
    print(f"{combined_failed_cells.sum()} unique cells will be removed based on combined QC thresholds.")

    # Apply filtering
    adata = adata[~combined_failed_cells, :].copy()

    # Filter genes expressed in too few cells
    genes_cells = (adata.X > 0).sum(axis=0) < min_cells
    print(f"{genes_cells.sum()} genes are expressed in fewer than {min_cells} cells. Will be filtered.")

    sc.pp.filter_genes(adata, min_cells=min_cells)

    print(f"Filtered AnnData has {adata.n_obs} cells and {adata.n_vars} genes")

    # Recalculate QC metrics for the filtered data
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    return adata


def visualize_QC_measures(adata, title="QC Measures"):
    """
    Visualize QC measures in violin and scatter plots.
    """
    qc_keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]

    def plot_scatter(single_adata, color_key, title):
        sc.pl.scatter(
            single_adata,
            x="total_counts",
            y="n_genes_by_counts",
            color=color_key,
            title=title,
            show=True
        )

    if isinstance(adata, dict):
        print(f"Visualizations {title} for multiple datasets:")

        combined_data = adata[list(adata.keys())[0]].concatenate(
            *[adata[key] for key in list(adata.keys())[1:]],
            batch_key="dataset",
            batch_categories=list(adata.keys())
        )

        # Check for missing QC fields
        missing_keys = [key for key in qc_keys if key not in combined_data.obs.columns]
        if missing_keys:
            print(f"Warning: Missing the following QC fields in combined_data.obs: {missing_keys}")

        valid_keys = [key for key in qc_keys if key in combined_data.obs.columns]
        sc.pl.violin(combined_data, valid_keys, groupby="dataset", jitter=0.4, multi_panel=True, stripplot=False)

        # Scatter plots for each dataset and QC metric
        for key in ["pct_counts_mt", "pct_counts_ribo"]:
            if key in combined_data.obs.columns:
                for name, ad in adata.items():
                    if key in ad.obs.columns:
                        plot_scatter(ad, key, f"{name} - {key}")

    else:
        print(f"Visualizations {title}:")
        valid_keys = [key for key in qc_keys if key in adata.obs.columns]
        sc.pl.violin(adata, valid_keys, jitter=0.4, multi_panel=True)

        for key in ["pct_counts_mt", "pct_counts_ribo"]:
            if key in adata.obs.columns:
                plot_scatter(adata, key, key)

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
    norm=False,
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
    
    # If norm is True, use [norm_total_counts] instead of [total_counts]
    if norm:
        if "norm_total_counts" not in adata.obs.columns:
            raise ValueError("norm_total_counts not found in adata.obs. Please save this data normalize the data first.")
        ax.hist(adata.obs["norm_total_counts"], bins=bins, log=False, label="Normalized total counts in bins")

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
    cell_cycle_genes = [x.strip().lower().title() for x in open('/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/regev_lab_cell_cycle_genes.txt')]
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

def custom_volcano_plot(
    df,
    logfc_col='log2FoldChange',
    pval_col='padj',
    lfc_thresh=0.25,
    pval_thresh=0.05,
    label_genes=None,
    top_n=10,
    x_lim=None,
    y_lim=None,
    figsize=(8, 6),
    dpi=300,
    title=None,
    use_adjust_text=True,
):
    """Volcano plot with optional auto‑labelling.

    If *label_genes* is **None**, the *top_n* most significant (–log₁₀ p)
    red/blue genes are labelled automatically.
    Axis limits can be overridden via *x_lim* and *y_lim*.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from adjustText import adjust_text

    df = df.copy()
    df['-log10_pval'] = -np.log10(df[pval_col])

    # colour assignment
    df['color'] = 'gray'
    df.loc[(df[logfc_col] >  lfc_thresh) & (df[pval_col] < pval_thresh), 'color'] = 'red'
    df.loc[(df[logfc_col] < -lfc_thresh) & (df[pval_col] < pval_thresh), 'color'] = 'blue'

    # auto‑select genes if none supplied
    if label_genes is None:
        sig = df['color'] != 'gray'
        label_genes = df[sig].nlargest(top_n, '-log10_pval').index.tolist()

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.scatter(df[logfc_col], df['-log10_pval'], c=df['color'], alpha=0.7)

    # significance thresholds
    y_thr = -np.log10(pval_thresh)
    x_thr = lfc_thresh

    # optional axis limits
    if x_lim is not None:
        ax.set_xlim(*x_lim)
    if y_lim is not None:
        ax.set_ylim(*y_lim)

    # draw threshold lines spanning full axes
    ax.axhline(y=y_thr, linestyle='--', color='black')
    ax.axvline(x=x_thr,  linestyle='--', color='black')
    ax.axvline(x=-x_thr, linestyle='--', color='black')

    ax.set_xlabel('log₂ fold‑change')
    ax.set_ylabel('–log₁₀ adjusted p')
    if title:
        ax.set_title(title)

    # label genes
    up_texts, down_texts = [], []
    for gene in label_genes:
        if gene not in df.index:
            continue
        x, y = df.at[gene, logfc_col], df.at[gene, '-log10_pval']
        t = ax.text(x, y, gene, fontsize=9,
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
        if x > 0:
            up_texts.append(t)
        else:
            down_texts.append(t)

    # directional nudging using adjustText
    if use_adjust_text:
        if up_texts:
            adjust_text(up_texts, ax=ax,
                        only_move={'static': 'yx+', 'text': 'yx+'},
                        ensure_inside_axes=True,
                        explode_radius=6,
                        pull_threshold=3,
                        force_text=[5, 5],
                        force_static=[2, 1],
                        force_pull=[0.1, 0.1],
                        arrowprops=dict(arrowstyle='->', color='black', shrinkA=4))
        if down_texts:
            adjust_text(down_texts, ax=ax,
                        only_move={'static': 'yx-', 'text': 'yx-'},
                        ensure_inside_axes=True,
                        explode_radius=6,
                        pull_threshold=3,
                        force_text=[5, 5],
                        force_static=[2, 1],
                        force_pull=[0.1, 0.1],
                        arrowprops=dict(arrowstyle='->', color='black', shrinkA=4))

    plt.tight_layout()
    plt.show()
    return fig


def summarize_celltype_results(cell_type, model, anova_table, tukey, output_dir):
    import pandas as pd
    import re

    # --- 3-Way ANOVA p ---
    try:
        interaction_p = anova_table.loc['C(apoe):C(disease):C(treatment)', 'PR(>F)']
    except KeyError:
        interaction_p = None

    # --- Largest effect from OLS model ---
    model_summary = model.summary2().tables[1]
    model_summary = model_summary.drop('Intercept', errors='ignore')
    largest_term_raw = model_summary['Coef.'].abs().idxmax()
    effect_size = model_summary.loc[largest_term_raw, 'Coef.']
    term_pval = model_summary.loc[largest_term_raw, 'P>|t|']

    # Clean model term for readability
    term = largest_term_raw.replace('C(', '').replace(')', '')
    term = term.replace(':', ' × ')
    term = re.sub(r'\[T\.(.*?)\]', r'=\1', term)

    # --- Tukey HSD ---
    res = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    res['meandiff'] = res['meandiff'].astype(float)
    res['p-adj'] = res['p-adj'].astype(float)
    sig_res = res[res['p-adj'] < 0.05]

    # --- Create DataFrames ---
    top_df = pd.DataFrame({
        'Statistic': ['3-way ANOVA p', 'Largest Effect Term', 'Effect Size (%)', 'p-value'],
        'Value': [
            round(interaction_p, 4) if interaction_p is not None else 'NA',
            term,
            round(effect_size, 2),
            round(term_pval, 4)
        ]
    })

    tukey_df = pd.DataFrame(columns=['Group 1', 'Group 2', 'Δ (%)', 'Tukey p-adj'])
    if not sig_res.empty:
        tukey_df['Group 1'] = sig_res['group1']
        tukey_df['Group 2'] = sig_res['group2']
        tukey_df['Δ (%)'] = sig_res['meandiff'].round(1)
        tukey_df['Tukey p-adj'] = sig_res['p-adj'].round(4)
        tukey_df = tukey_df.sort_values(by='Δ (%)', ascending=False)

    # --- Write two separate CSVs ---
    safe_ct = re.sub(r'\W+', '_', cell_type)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{safe_ct}_summary.csv").write_text(top_df.to_csv(index=False))
    (output_dir / f"{safe_ct}_tukey.csv").write_text(tukey_df.to_csv(index=False))

    print(f"Saved summary CSV for {cell_type} to {output_dir}/{safe_ct}_summary.csv")
    print(f"Saved Tukey CSV for {cell_type} to {output_dir}/{safe_ct}_tukey.csv")

def bivariate_quadrant_plot(
    df1, df2,
    label='log2FoldChange',
    label_genes=None,
    lfc_thresh=0.25,
    x_label=None, y_label=None,
    x_lim=None, y_lim=None,
    title=None,
    figsize=(6, 6), dpi=300,
    use_adjust_text=True
):
    """
    Bivariate plot of two log2FC vectors, coloring by quadrant based on log2FC thresholds only.
    Inclusion criteria: gene must be present in both df1 and df2 and have no NaNs in logFC or p-value.
    Coloring:
      - Red   = logFC_x >  lfc_thresh AND logFC_y >  lfc_thresh
      - Blue  = logFC_x < -lfc_thresh AND logFC_y < -lfc_thresh
      - Purple= one logFC > lfc_thresh and the other < -lfc_thresh
      - Gray  = any gene not meeting both |logFC| > lfc_thresh
    If label_genes is None, automatically label up to 10 genes per colored quadrant
    with the largest |logFC_x + logFC_y|.
    Optionally, pass x_label and y_label to override axis titles.
    Optionally, pass x_lim and y_lim as tuples to set axis limits.
    Adds Pearson correlation annotation and a dotted regression line,
    calculated only on the colored points.
    """
    # Handle axis labels
    if x_label is None or y_label is None:
        import inspect
        frame = inspect.currentframe().f_back
        names = {id(val): name for name, val in frame.f_locals.items()}
        if x_label is None:
            x_label = names.get(id(df1), 'Comparison A')
        if y_label is None:
            y_label = names.get(id(df2), 'Comparison B')

    # Intersect indices and build working DataFrame
    common = df1.index.intersection(df2.index)
    df = pd.DataFrame({
        'logFC_x': df1.loc[common, label],
        'logFC_y': df2.loc[common, label],
        'pval_x':  df1.loc[common, 'padj'],
        'pval_y':  df2.loc[common, 'padj']
    }).dropna()

    # Mask of genes passing logFC thresholds in BOTH comparisons
    mask = (
        (df['logFC_x'].abs() > lfc_thresh)
        & (df['logFC_y'].abs() > lfc_thresh)
    )

    # Assign colors based only on logFC quadrant
    def quadrant_color(r):
        if not mask.loc[r.name]:
            return 'gray'
        if r['logFC_x'] >  lfc_thresh and r['logFC_y'] >  lfc_thresh:
            return 'red'
        if r['logFC_x'] < -lfc_thresh and r['logFC_y'] < -lfc_thresh:
            return 'blue'
        return 'purple'

    df['color'] = df.apply(quadrant_color, axis=1)

    # If no explicit labels, pick top genes by |logFC_x + logFC_y| in each colored quadrant
    if label_genes is None:
        scores = (df['logFC_x'] + df['logFC_y']).abs()
        computed_labels = []
        for col in ['red', 'blue', 'purple']:
            genes_in_col = df.index[df['color'] == col]
            if not genes_in_col.empty:
                top_genes = scores.loc[genes_in_col].nlargest(10).index.tolist()
                computed_labels.extend(top_genes)
        label_genes = computed_labels

    # Filter only colored points for correlation and regression
    df_colored = df[df['color'] != 'gray']

    # Compute Pearson correlation and regression line on colored points
    from scipy.stats import pearsonr, linregress
    if not df_colored.empty:
        r_val, p_val = pearsonr(df_colored['logFC_x'], df_colored['logFC_y'])
        slope, intercept, _, _, _ = linregress(df_colored['logFC_x'], df_colored['logFC_y'])
    else:
        r_val, p_val, slope, intercept = [float('nan')]*4

    # Plot
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.scatter(df['logFC_x'], df['logFC_y'], c=df['color'], alpha=0.7)

    # Plot infinite regression line
    if not df_colored.empty:
        ax.axline((0, intercept), slope=slope, linestyle=':', color='black', linewidth=1)

    # Annotate correlation
    ax.text(
        0.05, 0.95,
        f"r={r_val:.2f}\np={p_val:.2g}",
        transform=ax.transAxes,
        ha='left', va='top', fontsize=10,
        bbox=dict(facecolor='white', edgecolor='none', pad=2)
    )

    # Draw threshold lines
    ax.axvline( lfc_thresh, color='black', linestyle='--')
    ax.axvline(-lfc_thresh, color='black', linestyle='--')
    ax.axhline( lfc_thresh, color='black', linestyle='--')
    ax.axhline(-lfc_thresh, color='black', linestyle='--')

    # Apply optional axis limits
    if x_lim is not None:
        ax.set_xlim(*x_lim)
    if y_lim is not None:
        ax.set_ylim(*y_lim)

    ax.set_xlabel(f'{x_label} (log₂FC)')
    ax.set_ylabel(f'{y_label} (log₂FC)')
    if title:
        ax.set_title(title)

    # Gene labeling
    if label_genes:
        texts = []
        for g in label_genes:
            if g in df.index and df.at[g, 'color'] != 'gray':
                x, y = df.at[g, 'logFC_x'], df.at[g, 'logFC_y']
                texts.append(ax.text(x, y, g, fontsize=9))
        if use_adjust_text and texts:
            from adjustText import adjust_text
            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(arrowstyle='->', color='black', shrinkA=5),
                expand_axes=True
            )
        for t in texts:
            t.set_bbox(dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    plt.tight_layout()
    return fig


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
        Dictionary of conditions that remain constant (e.g. {'apoe': 'E4', 'diesase': 'WT'}).
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
    
    bars = ax.bar(
        bar_positions,
        overall_means,
        color=['lightblue', 'lightgreen'],
        alpha=0.6,
        width=0.6,
        label='Overall Mean Expression'
    )
    
    # Overlay individual sample means as dots with jitter for clarity
    jitter_std = 0.08
    for i, sample_means in enumerate([sample_means1, sample_means2]):
        if not sample_means.empty:
            x_positions = np.random.normal(loc=i, scale=jitter_std, size=len(sample_means))
            ax.scatter(
                x_positions,
                sample_means.values,
                color='black',
                alpha=0.8,
                zorder=10,
                label='Sample Expression' if i == 0 else ""
            )
    
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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_quad_contrast(adata, gene, disease, dc_vehlps_results):
    """
    Plot a bar plot comparing gene expression under four contrast groups:
      - Two APOE conditions (E3 and E4) and two treatment levels (VEHICLE and LPS),
        while holding disease state constant (e.g. 'WT' or 'FAD').
      
      For each group:
        - Cells are filtered by APOE, treatment, and the constant disease state.
        - Cells are grouped by 'Classification' (i.e. sample) and the mean normalized count 
          for the gene (assumed to be stored in adata.layers['log_norm']) is computed.
        - A bar is drawn with height equal to the overall mean (across samples) and individual
          sample means are overlaid as jittered diamond markers.
      
      The function uses a custom color palette for each group.
      
      Parameters:
      -----------
      adata : AnnData
          The AnnData object with expression data and observations.
      gene : str
          The gene of interest.
      disease : str
          The disease state (e.g. 'WT' or 'FAD') to hold constant.
      dc_vehlps_results : dict
          A dictionary of DESeq2 results with keys of the form f"{apoe}{disease}_vehlps".
      
      Returns:
      --------
      None
          Displays the annotated bar plot.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Define color palettes for E3 and E4.
    e3_wt_palette = sns.color_palette("Blues", 2)    # two shades for WT in E3
    e3_fad_palette = sns.color_palette("Greens", 2)   # two shades for FAD in E3
    e4_wt_palette = sns.color_palette("Oranges", 2)   # two shades for WT in E4
    e4_fad_palette = sns.color_palette("Purples", 2)   # two shades for FAD in E4


    if disease == 'WT':
        e3_palette = e3_wt_palette
        e4_palette = e4_wt_palette
    elif disease == 'FAD':
        e3_palette = e3_fad_palette
        e4_palette = e4_fad_palette

    # For this function, we assume the desired groups are:
    # E3_{disease}_VEHICLE, E3_{disease}_LPS, E4_{disease}_VEHICLE, E4_{disease}_LPS
    desired_order = [
        f"E3_{disease}_VEHICLE",
        f"E3_{disease}_LPS",
        f"E4_{disease}_VEHICLE",
        f"E4_{disease}_LPS"
    ]
    


    # Map groups to colors using the appropriate palette.
    # (Here, we assume for simplicity that E3 uses the WT palette and E4 uses the WT palette.
    # Adjust if you want to incorporate FAD as a separate color family.)
    color_dict = {
        f"E3_{disease}_VEHICLE": e3_palette[0],
        f"E3_{disease}_LPS": e3_palette[1],
        f"E4_{disease}_VEHICLE": e4_palette[0],
        f"E4_{disease}_LPS": e4_palette[1],
    }
    
    # For positioning, we place the E3 groups close together and the E4 groups with a gap.
    positions = {
        f"E3_{disease}_VEHICLE": 0,
        f"E3_{disease}_LPS": 1,
        f"E4_{disease}_VEHICLE": 3,
        f"E4_{disease}_LPS": 4,
    }
# Define a helper to subset adata for a given APOE and treatment, keeping disease constant.
    def subset_for_conditions(apoe, treatment, disease):
        mask = (
            (adata.obs['apoe'] == apoe) &
            (adata.obs['treatment'] == treatment) &
            (adata.obs['disease'] == disease)
        )
        return adata[mask].copy()
    
    # Helper function to compute sample means for a given subset.
    def get_sample_means(adata_subset, gene):
        # Use log-normalized counts (assumed stored in adata.layers['log_norm'])
        adata_subset.X = adata_subset.layers['log_norm'].copy()
        if gene not in adata_subset.var_names:
            raise ValueError(f"Gene '{gene}' not found in adata.var_names.")
        gene_expr = adata_subset[:, gene].X
        if hasattr(gene_expr, "toarray"):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = np.array(gene_expr).flatten()
        # Create a DataFrame with expression and sample identity (Classification)
        df = pd.DataFrame({
            'normalized_counts': gene_expr,
            'Classification': adata_subset.obs['Classification'].values
        })
        # Compute mean expression per sample
        sample_means = df.groupby('Classification')['normalized_counts'].mean()
        return sample_means

    # Define the two APOE values and treatment levels.
    apoes = ['E3', 'E4']
    treatments = ['VEHICLE', 'LPS']

    overall_means = []         # To store the overall mean for each group.
    overall_stds = []          # To store the standard deviation.
    all_sample_data = []       # To store individual sample means with group labels.
    
    # Loop over the groups in desired order.
    for group in desired_order:
        # Parse group label (expected format: "E3_{disease}_VEHICLE")
        parts = group.split("_")
        apoe = parts[0]
        treatment = parts[2]
        # Subset data and compute sample means.
        adata_subset = subset_for_conditions(apoe, treatment, disease)
        sample_means = get_sample_means(adata_subset, gene)
        mean_val = sample_means.mean() if not sample_means.empty else np.nan
        std_val = sample_means.std() if not sample_means.empty else np.nan
        overall_means.append(mean_val)
        overall_stds.append(std_val)
        # Record individual sample means.
        for val in sample_means.values:
            all_sample_data.append({'group': group, 'percent': val})
    
    # Create DataFrames for group statistics and individual sample data.
    group_stats = pd.DataFrame({
        'group': desired_order,
        'x': [positions[g] for g in desired_order],
        'mean': overall_means,
        'std': overall_stds
    })
    df_long = pd.DataFrame(all_sample_data)
    
    # --- Plotting using matplotlib's bar plot style ---
    plt.figure(figsize=(8, 4))
    ax = plt.gca()
    bar_width = 0.35

    # Plot the bars with error bars.
    ax.bar(
        group_stats['x'],
        group_stats['mean'],
        yerr=group_stats['std'],
        width=bar_width,
        color=[color_dict[g] for g in desired_order],
        capsize=5,
        edgecolor='black'
    )

    # Overlay individual sample points as diamonds with slight x jitter.
    for g in desired_order:
        subset = df_long[df_long['group'] == g]
        x_pos = positions[g]
        jitter = np.random.uniform(-0.05, 0.05, size=len(subset))
        ax.plot(
            x_pos + jitter,
            subset['percent'],
            'D',
            markersize=7,
            markeredgecolor='black',
            markerfacecolor=color_dict[g],
            linestyle='none'
        )

    # Set custom x-axis ticks and labels.
    ax.set_xticks([positions[g] for g in desired_order])
    ax.set_xticklabels([g.replace("_", "-") for g in desired_order], rotation=45)
    
    ax.set_xlabel("Group (apoe-disease-treatment)")
    ax.set_ylabel(f"Mean Log-Normalized expression for {gene}")
    ax.set_title(f"Expression of {gene}\nDisease: {disease}")
    
    # --- Annotate with p-values and significance stars ---
    def get_significance(pval):
        if np.isnan(pval):
            return 'NA'
        if pval < 0.001:
            return '***'
        elif pval < 0.01:
            return '**'
        elif pval < 0.05:
            return '*'
        else:
            return 'ns'

    # Compute annotations for each APOE contrast.
    annotations = {}
    for apoe in apoes:
        key = f"{apoe}{disease}_vehlps"  # expected key format: e.g. "E3WT_vehlps"
        if key in dc_vehlps_results and gene in dc_vehlps_results[key].index:
            p_val = dc_vehlps_results[key].loc[gene, 'padj']
        else:
            p_val = np.nan
        sig = get_significance(p_val)
        # Group center is the mean of the two positions for that APOE.
        group_positions = [positions[f"{apoe}_{disease}_{t}"] for t in treatments]
        group_center = np.mean(group_positions)
        annotations[apoe] = (group_center, p_val, sig)

    # Determine a y-value for annotation (a little above the highest bar).
    valid_means = [m for m in overall_means if not np.isnan(m)]
    y_max = max(valid_means) * 1.1 if valid_means else 1

    # Only annotate if the p-value is neither 1 nor N/A.
    for apoe in apoes:
        group_center, p_val, sig = annotations[apoe]
        if np.isnan(p_val) or p_val == 1:
            continue  # Do not show annotation if p-val is N/A or 1
        annotation_text = f"{sig}\np = {p_val:.3g}"
        ax.text(group_center, y_max, annotation_text, ha='center', va='bottom', fontsize=14)

    ax.set_ylim(0, y_max * 1.5)

    
    plt.tight_layout()

    return plt.gcf()
