# This file contains the functions to be used in the notebooks for the project

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import anndata as ad
import pandas as pd
from matplotlib_venn import venn2
from adjustText import adjust_text
from typing import List, Optional, Union

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


def hash_demulitplex(adata, hashtag_prefix='Hashtag'):
    """
    Function to demultiplex droplets using the Hashsolo package.
    
    Parameters:
    adata (AnnData): The input AnnData object.
    hashtag_prefix (str): The prefix for hashtag genes in the var gene_ids.
    
    Returns:
    AnnData: The updated AnnData object with demultiplexing results.
    """
    
    from solo import hashsolo

    # Check if hashtag genes are present
    hashtag_genes = adata.var.gene_ids.str.startswith(hashtag_prefix)
    if sum(hashtag_genes) == 0:
        print(f"No genes with prefix '{hashtag_prefix}' found.")
        return adata

    # Extract hashtag data
    hto = adata[:, hashtag_genes].copy()
    
    # Run HashSolo
    hashsolo.hashsolo(hto)

    # Print the number of predicted singlets, doublets, and negatives
    singlets = sum(hto.obs['most_likely_hypothesis'] == 1)
    doublets = sum(hto.obs['most_likely_hypothesis'] == 2)
    negatives = sum(hto.obs['most_likely_hypothesis'] == 0)
    print(f'Number of predicted singlets: {singlets}')
    print(f'Number of predicted doublets: {doublets}')
    print(f'Number of predicted negatives: {negatives}')
    
    # Print the number of cells classified to each hashtag
    classifications = hto.obs['Classification'].value_counts()
    print('Number of cells for each hashtag:')
    for classification, count in classifications.items():
        print(f'{classification}: {count}')
    
    # Map classifications back to original AnnData object
    adata.obs['Classification'] = hto.obs['Classification'].map(lambda x: adata.var.index[adata.var['gene_ids'] == x][0] if x in adata.var['gene_ids'].values else x)

    # Update the original AnnData object with the most likely hypothesis and classification
    adata.obs['most_likely_hypothesis'] = hto.obs['most_likely_hypothesis']
    
    return adata

# Example usage:
# adata = sc.read_h5ad('path_to_your_anndata.h5ad')
# adata = demultiplex_droplets(adata)

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
def find_genes(adata):
    #define some lists of genes of interest for the analysis.
    Genes_of_interest = []
    homeo_genes = ['Sall1', 'Cx3cr1', 'P2ry12', 'P2ry13', 'Olfml3', 'Tmem119', 'Cd68', 'Itgam', 'Cst3']
    DAM_genes = ['Apoe', 'B2m', 'Ctsb', 'Lpl', 'Cst7', 'Tyrobp', 'Trem2', 'Cd9', 'Itgax', 'Cd63', 'Fth1', 'Spp1', 'Axl', 'Mertk', 'Mdk', 'Ptn', 'Lgals3']
    Mapk_genes = ['Syk', 'Prkca', 'Mef2c', 'Elk4', 'Trp53', 'Mertk', 'Axl', 'Mapk14', 'Mapk1', 'Gadd45a',]
    Nfkb_genes = ['Nfkb1','Nfkb2', 'Tlr4', 'Nfkbia', 'Cd40', 'Tlr4', 'Rela', 'Relb', 'Irf1', 'Irf3', 'Irf4']
    Macroph_genes = ['Mrc1', 'Mrc2', 'Ccr2', 'Ly6c2', 'Lyz2', 'Vim', 'Ifi204', 'S100a10', 'Msrb1']
    DC_genes = ['Cd14', 'Cd1a', 'Cd209', 'Itgam', 'Cd207', 'Fcer1', '33D1', 'eyfp', 'Flt3', 'Csf1r', 'Id2', 'Clec4k', 'Cd8a']
    TIM_genes = ['Malat1', 'Tmx4', 'Jund', 'Btg2', 'Fosb', 'Fos','Jun', 'Klf4', 'Klf2', 'Gm34455']

    # make a dictionary of the lists to iterate through
    gene_lists = {"Homeo": homeo_genes, "DAM": DAM_genes, "Mapk": Mapk_genes, "Nfkb": Nfkb_genes, "Macroph": Macroph_genes, "DC": DC_genes, "TIM": TIM_genes}
    found_genes = {key: [] for key in gene_lists.keys()}

    # loop through gene lists and genes, see if it is present in the passed adata, if so append it to a list+dictionary 
    # so we have a dicitonary of lists of only genes in the dataset, as well as a list of all genes from the different types
    for gene_list_name, gene_list in gene_lists.items():
        for gene in gene_list:
            if gene in adata.var.index:
                print(f" ✓ {gene} is in adata.var")
                Genes_of_interest.append(gene)
                found_genes[gene_list_name].append(gene)
            else:
                print(f"{gene} is not in adata.var")

    return Genes_of_interest, found_genes

import scanpy as sc

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
    # Create a DataFrame with necessary data
    df = pd.DataFrame({
        'logfc': adata.uns['rank_genes_groups']['logfoldchanges'][group],
        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][group],
        'gene_names': adata.uns['rank_genes_groups']['names'][group]
    })

    # Calculate -log10 adjusted p-values
    df['neg_log_pvals'] = -np.log10(df['pvals_adj'])



    #Check neg_log_pvals for inf and nan values
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

    # Annotate top N up-regulated genes
    up_genes = df[df['significance'] == 'Up-regulated']
    top_up_genes = up_genes.nsmallest(top_n, 'pvals_adj')

    for _, row in top_up_genes.iterrows():
        ax.text(row['logfc'], row['neg_log_pvals'], row['gene_names'], fontsize=8)

    # Annotate top N down-regulated genes
    down_genes = df[df['significance'] == 'Down-regulated']
    top_down_genes = down_genes.nsmallest(top_n, 'pvals_adj')

    for _, row in top_down_genes.iterrows():
        ax.text(row['logfc'], row['neg_log_pvals'], row['gene_names'], fontsize=8)

    # Create a dictionary of top N up- and down-regulated genes
    annotated_genes = {
        'Up-regulated': top_up_genes['gene_names'].tolist(),
        'Down-regulated': top_down_genes['gene_names'].tolist()
    }


    # Add threshold lines
    ax.axvline(-logfc_threshold, color="grey", linestyle="--")
    ax.axvline(logfc_threshold, color="grey", linestyle="--")
    ax.axhline(-np.log10(pval_threshold), color="grey", linestyle="--")

    # Set labels and title
    ax.set_xlabel("Log Fold Change (logFC)")
    ax.set_ylabel("-Log10 Adjusted P-Value (-logFDR)")
    ax.set_xlim(-10, 10)
    ax.set_ylim(0, df['neg_log_pvals'].max() + 1)
    ax.legend()
    ax.set_title(f"Volcano Plot for {group}")

    # Improve layout
    plt.tight_layout()
    plt.show()

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