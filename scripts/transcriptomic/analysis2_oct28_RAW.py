# %% [markdown]
# # Analysis
# 
# This notebook will constitute the analysis of the single cell data however the non-microglia will be removed from the beginning allowing re-normalization, re-QC, re-Feature selection to ideally prevent any skew from different cell types. 
# 
# the main steps in this notebook will be:
# 	1.	Demultiplexing and Barcode Removal (on individual datasets)
# 	2.	QC (on individual datasets)
# 	3.	Remove Non-Microglia (on individual datasets)
# 	4.	Merge Datasets
# 	5.	Normalization
# 	6.	Regress MitoCell cycle
# 	7.	Cell Cycle
# 	8.	Dimensionality reduction
# 			- PCA
# 			- nieghbor
# 			- UMAP
# 	9.	Cluster
# 	10.	Decoupler(DEseq)
# 	11.	Single variable comparison for DEGs; volcano plots.
# 			- Overlap between gene lists
# 	12.	GO (on decoupler list and DEG lists?)
# 	13. GSEA 
# 
# 
# These steps will be turned into a script which will generate reports for different variables at some point. 
# 

# %% [markdown]
# ## Initializaton

# %%
import numpy as np
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
from bbknn import bbknn
import matplotlib.pyplot as plt
import seaborn as sns
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import importlib
import gc

# Modify the sys.path to include the notebook directory to import functions
import sys
import os

notebook_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
if notebook_dir not in sys.path:
    sys.path.append(notebook_dir)

import func_lib as f

importlib.reload(f)

# ignore pd performance warning
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# Set the random seed for reproducibility
seed = 2

# %% [markdown]
# ## Data read in 

# %%

adatas = {
    "D1" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/collect_tube_1_batch_3_June_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    "D2" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/collect_tube_2_batch_3_june_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    "D3" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/Isol_Microglia_EFAD_TD_august_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    "D4" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/MICROGLIA_E3E4FAD_TD_23_05_2023_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    "D5" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/Tube_1_July_20_TD_YYRFC_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    "D6" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/Tube2_July_20_TD_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    "D7" : sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/Summer_proj_2024/Results_HT/Tube3_july_20_TD_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False),
    }

# genotype and treatment map for samples
genotype_map = {"K191" : "E3FAD", "K192" : "E3FAD", "K190" : "E3WT", "K193" : "E3WT", 
                "G118" : "E4FAD", "G119" : "E4FAD", "G120" : "E4FAD", "G122" : "E4WT",
                "K231" : "E3FAD", "K233" : "E3FAD", "K234" : "E3FAD", "K219" : "E4FAD",
                "G143" : "E4WT", "G146" : "E4WT", "K232" : "E3WT", "G180" : "E4WT", 
                "G183" : "E4WT", "G184" : "E4WT", "G185" : "E4WT", "G186" : "E4FAD",
                "K248" : "E3WT", "K249" : "E3WT", "K250" : "E3FAD", "K252" : "E3FAD", 
                "G177" : "E4FAD", "G178" : "E4FAD", "G179" : "E4FAD", "G145" : "E4FAD",
                "K257" : "E3WT", "K261" : "E3WT", "K262" : "E3WT", "K268" : "E3FAD",
                "K283" : "E3WT", "K284" : "E3FAD"}
treatment_map = {"K191" : "LPS", "K192" : "LPS", "K190" : "LPS", "K193" : "LPS", 
                "G118" : "VEHICLE", "G119" : "VEHICLE", "G120" : "VEHICLE", "G122" : "VEHICLE",
                "K231" : "VEHICLE", "K233" : "VEHICLE", "K234" : "VEHICLE", "K219" : "LPS",
                "G143" : "LPS", "G146" : "LPS", "K232" : "VEHICLE", "G180" : "LPS", 
                "G183" : "VEHICLE", "G184" : "VEHICLE", "G185" : "VEHICLE", "G186" : "VEHICLE",
                "K248" : "LPS", "K249" : "LPS", "K250" : "LPS", "K252" : "LPS", 
                "G177" : "LPS", "G178" : "LPS", "G179" : "LPS", "G145" : "LPS",
                "K257" : "VEHICLE", "K261" : "VEHICLE", "K262" : "VEHICLE", "K268" : "VEHICLE",
                "K283" : "LPS", "K284" : "LPS"}


# Load the CSV file containing manual annotations
manual_annotations = pd.read_csv('manual_annotation.csv', index_col=0)




# %% [markdown]
# ## Demux

# %%
#Demultiplex hashtags from datasets
for adata_key, adata in adatas.items():
    print(f"{adata_key}:")
    f.hash_demulitplex(adata)
    print("\n")

# %%
#remove doublets and empty droplets from datasets 
for adata_key, adata in adatas.items():
    print(f"Removing doublets and empty droplets from {adata_key}:")
    print(f"number of cells before: {adata.obs.shape[0]}")
    doublets = adata.obs['most_likely_hypothesis'] == 2
    negatives = adata.obs['most_likely_hypothesis'] == 0
    adata = adata[~doublets & ~negatives]
    adatas[adata_key] = adata
    print(f"number of cells after: {adata.obs.shape[0]}")
    print("\n")

# %%
# Apply genotype and treatment mapping for each adata separately
for adata_key, adata in adatas.items():
    print(f"Processing {adata_key}:")

    # Apply genotype mapping
    adata.obs["genotype"] = adata.obs["Classification"].map(lambda x: genotype_map.get(x, None))

    # Apply APOE mapping (first two characters of the genotype value)
    adata.obs["apoe"] = adata.obs["Classification"].map(lambda x: genotype_map.get(x, None)[:2] if genotype_map.get(x, None) else None)

    # Apply disease mapping (characters after the first two characters of the genotype value)
    adata.obs["disease"] = adata.obs["Classification"].map(lambda x: genotype_map.get(x, None)[2:] if genotype_map.get(x, None) else None)

    # Apply treatment mapping
    adata.obs["treatment"] = adata.obs["Classification"].map(lambda x: treatment_map.get(x, None))

    # Apply genotype_treatment mapping
    adata.obs['genotype_treatment'] = adata.obs['genotype'] + '_' + adata.obs['treatment']


    # Print statistics for each `adata`
    print(f"Number of cells of each genotype in {adata_key}:")
    print(adata.obs["genotype"].value_counts())
    print("\n")
    print(f"Number of cells of each treatment in {adata_key}:")
    print(adata.obs["treatment"].value_counts())
    print("\n")
    print(f"Number of cells of each genotype and treatment in {adata_key}:")
    print(adata.obs.groupby(["genotype", "treatment"]).size())
    print("\n")
    print(adata.obs.groupby(["Classification", "genotype"]).size().sort_values(ascending=False))
    print("\n")

# %%
# Step 1: Group and calculate group sizes for each 'Classification' and 'disease', then compute mean sizes
for adata_key, adata in adatas.items():
    print(f"Processing {adata_key}:")
    
    # Group by 'Classification' and 'disease', then calculate the size of each group
    group_sizes = adata.obs.groupby(["Classification", "disease"]).size().reset_index(name='size')

    # Group by 'disease' and calculate the mean size for each classification within each disease group
    mean_sizes = group_sizes.groupby("disease")["size"].mean().reset_index()

    # Print the mean sizes
    print(f"Mean sizes for {adata_key}:")
    print(mean_sizes)
    print("\n")

# Step 2: Remove hashtag genes from each adata separately
for adata_key, adata in adatas.items():
    print(f"Removing hashtag genes from {adata_key}:")
    
    # Find hashtag genes
    hashtag_genes = [gene for gene in adata.var_names if any(ht in gene for ht in genotype_map.keys())]

    # Remove hashtag genes
    adata = adata[:, ~adata.var_names.isin(hashtag_genes)]

    # Update the adatas dictionary with the modified adata
    adatas[adata_key] = adata

    # Check if hashtag genes are present in the .var DataFrame
    print(f"Checking for presence of hashtag genes in {adata_key}:")
    for hashtag in genotype_map.keys():
        if hashtag in adata.var_names:
            print(f"{hashtag} is present in the .var DataFrame.")
        else:
            print(f"{hashtag} is not present in the .var DataFrame.")
    print("\n")

# %% [markdown]
# ## QC

# %%
# Print Quality Control metric for each dataset
for adata_key, adata in adatas.items():
    print(f"{adata_key}:")
    f.summarize_adata(adata, mt_gene_prefix="mt-", ribo_gene_prefixes=("Rps", "Rpl"), min_counts=2000, min_genes=300, min_cells=3)
    print("\n")

# %%
# Visualize pre-filtering QC measures
f.visualize_QC_measures(adatas)

# %%
# Basic QC thresholds
QC_filtered_adatas = {}
for adata_key, adata in adatas.items():
    print(f"{adata_key}:")
    QC_filtered_adatas[adata_key] =f.QC_filter_adata(adata, mt_threshold=5, ribo_threshold=5, min_counts=2000, min_genes=300, min_cells=3)
    print("\n")

# %%
f.visualize_QC_measures(QC_filtered_adatas)

# %% [markdown]
# ## Remove Non-Microglia

# %%
# Initialize counters for total cells before and after removal
total_cells_before = 0
total_cells_after = 0

# Loop through each adata in the QC_filtered_adatas dictionary
for i, (adata_key, adata) in enumerate(QC_filtered_adatas.items()):
    print(f"Processing {adata_key}:")

    # Get the cells that are in the adata
    adata_cell_ids = adata.obs.index

    # Add the number of cells in the current adata to the total before removal
    total_cells_before += len(adata_cell_ids)

    # Define the expected suffix based on the current index (e.g., _0 for first, _1 for second)
    expected_suffix = f'_{i}'

    # Filter manual_annotations to include only the cells with the matching suffix
    manual_annotation_cells_with_suffix = manual_annotations[manual_annotations.index.str.endswith(expected_suffix)]

    # Remove the suffix from manual_annotations to match the format in adata_cell_ids
    manual_annotation_cells = manual_annotation_cells_with_suffix.index.str.replace(expected_suffix, '')

    # Find the cells present in both adata and manual_annotations (after suffix removal)
    matching_cells = adata_cell_ids.intersection(manual_annotation_cells)

    # Find cells annotated as 'Other' in the matching manual_annotations
    other_cells = manual_annotation_cells_with_suffix[manual_annotation_cells_with_suffix['manual_annotation'] == 'Other'].index

    # Ensure the 'Other' cells have the suffix properly formatted before removal
    other_cells_no_suffix = other_cells.str.replace(expected_suffix, '')

    # Remove these 'Other' cells from the adata
    print(f"Removing {len(other_cells)} cells annotated as 'Other' from {adata_key}")
    adata = adata[~adata.obs.index.isin(other_cells_no_suffix)].copy()

    # Add the number of cells in the current adata to the total after removal
    total_cells_after += adata.obs.shape[0]

    # Update the adatas dictionary with the filtered adata
    QC_filtered_adatas[adata_key] = adata

    print(f"Number of cells after removal in {adata_key}: {adata.obs.shape[0]}")
    print("\n")

# Print the total cell counts before and after removal
print(f"Total number of cells before removing 'Other' cells: {total_cells_before}")
print(f"Total number of cells after removing 'Other' cells: {total_cells_after}")

# %% [markdown]
# ## Merge

# %%
adata = ad.concat(list(QC_filtered_adatas.values()), axis=0, label="sample", join="outer", merge="unique")

print(adata)

# %%

# remove and gc the QC_filtered_adatas and adatas
del QC_filtered_adatas
del adatas
gc.collect()

# %% [markdown]
# ## Normalize & HVG

# %%
# before normalizing, save the raw counts
adata.layers["counts"] = adata.X.copy()

# visualize the distribution of total counts in cells
f.plot_total_counts_vs_cells(adata, bins=250, zoom_x_range=(0, 5000), zoom_y_range=(0, 1000))

# %%
# normalize and log the data and store that in a different layer too incase we need to switch between them.
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['log_norm'] = adata.X.copy()
#so we can plot a histogram of the total counts/cells after norm
adata.obs["norm_total_counts"] = adata.X.sum(axis=1)    
f.plot_total_counts_vs_cells(adata, bins=250, zoom_x_range=(0, 5000), zoom_y_range=(0, 1000))

# %%
sc.pp.highly_variable_genes(adata, batch_key="sample")
sc.pl.highly_variable_genes(adata)
#make some variables to print a nice little summary
num_highly_variable_genes_t = adata.var['highly_variable'].sum()
total_genes_t = adata.var.shape[0]

print(f"{num_highly_variable_genes_t} out of {total_genes_t} genes are considered highly variable.")

# %% [markdown]
#  ## Cell cycle

# %%
adata.X = adata.layers["counts"].copy()
adata=f.annotate_cellcycle_mouse(adata)

df = adata.obs[['phase', 'genotype', 'apoe', 'treatment', 'disease']]


# Calculate and print the total percentages for each phase
total_phase_percentages = {}
total_count = adata.obs.shape[0]

for p in adata.obs['phase'].unique():
    phase_count = len(adata.obs[adata.obs['phase'] == p])
    phase_percentage = round(phase_count / total_count * 100, 2)
    total_phase_percentages[p] = phase_percentage
    print(f"The percent of {p} phase in the total dataset is {phase_percentage}%")
print("\n")

# Calculate and print the percentages for each sample and compare to the total dataset
for s in adata.obs['sample'].unique():
    for p in adata.obs['phase'].unique():
        phase_count = adata[(adata.obs['phase'] == p) & (adata.obs['sample'] == s)].shape[0]
        sample_total_count = adata[adata.obs['sample'] == s].shape[0]
        sample_phase_percentage = round(phase_count / sample_total_count * 100, 2)
        total_phase_percentage = total_phase_percentages[p]
        difference = sample_phase_percentage - total_phase_percentage
        difference_statement = f"{difference:.2f}% higher" if difference > 0 else f"{-difference:.2f}% lower"
        print(f"{s} : The percent of {p} phase is {sample_phase_percentage}% ({difference_statement} than the total dataset)")
    print('\n')


for s in adata.obs['genotype'].dropna().unique():
    for p in adata.obs['phase'].unique():
        phase_count = adata[(adata.obs['phase'] == p) & (adata.obs['genotype'] == s)].shape[0]
        sample_total_count = adata[adata.obs['genotype'] == s].shape[0]
        sample_phase_percentage = round(phase_count / sample_total_count * 100, 2)
        total_phase_percentage = total_phase_percentages[p]
        difference = sample_phase_percentage - total_phase_percentage
        difference_statement = f"{difference:.2f}% higher" if difference > 0 else f"{-difference:.2f}% lower"
        print(f"{s} : The percent of {p} phase is {sample_phase_percentage}% ({difference_statement} than the total dataset)")
    print('\n')

# %% [markdown]
# ## Regress MT-genes

# %%
adata.X = adata.layers['log_norm'].copy()

# %%

adata = sc.pp.regress_out(adata, 'pct_counts_mt', copy=True)


# %% [markdown]
# ## Dimensionality Reduction

# %%

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50)

sc.pl.pca(
    adata,
    color=[
        "sample", "sample", "pct_counts_mt", "pct_counts_mt", 
        "genotype", "genotype", "treatment", "treatment", 
        "apoe", "apoe", "disease", "disease"
    ],
    dimensions=[
        (0, 1), (2, 3), (0, 1), (2, 3), 
        (0, 1), (2, 3), (0, 1), (2, 3), 
        (0, 1), (2, 3), (0, 1), (2, 3)
    ],
    ncols=2,
    size=2
)

# %%
sc.pp.neighbors(adata, n_pcs=25)

sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=['sample','genotype', 'apoe', 'treatment', 'pct_counts_mt'], ncols=2)

# %%
sc.pl.umap(adata, color=['sample', 'phase', 'total_counts', 'genotype', 'apoe', 'treatment', 'disease', 'Malat1'], ncols=2)

# %% [markdown]
# ## Cluster

# %%
i = 0.1
while i <= 2.1:
    key = f"leiden_{i}"
    if key not in adata.obs:
        sc.tl.leiden(
            adata,
            resolution=i,
            random_state=0,
            n_iterations=2,
            directed=False,
            key_added=key,
            flavor="igraph"
        )
    i += 0.2
    i = round(i, 2)

sc.pl.umap(adata, color=['leiden_0.1', 'leiden_0.3', 'leiden_0.5', 'leiden_0.7', 'leiden_0.9', 'leiden_1.1', 'leiden_1.3', 'leiden_1.5', 'leiden_1.7', 'leiden_1.9', 'leiden_2.1'], legend_loc="on data", ncols=3)

# %% [markdown]
# ## Pseudobulk w/ Decoupler

# %%


# Get pseudo-bulk profile
dc_adata = dc.get_pseudobulk(
    adata,
    sample_col='Classification',
    groups_col='genotype',
    layer='counts',
    mode='sum',
    min_cells=0,
    min_counts=0
)

dc.plot_psbulk_samples(dc_adata, groupby=['Classification', 'genotype', 'sample'], figsize=(12, 4))

# %%
# Store raw counts in layers
dc_adata.layers['counts'] = dc_adata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(dc_adata, target_sum=1e4)
sc.pp.log1p(dc_adata)
sc.pp.scale(dc_adata, max_value=10)
sc.tl.pca(dc_adata)

dc.swap_layer(dc_adata, 'counts', inplace=True)

# %%
sc.pl.pca(dc_adata, color=['disease', 'genotype', 'apoe', 'treatment', 'genotype_treatment', 'sample'], ncols=3, size=300)
sc.pl.pca(dc_adata, color=['disease', 'genotype', 'apoe', 'treatment'], dimensions=(2,3) , ncols=2, size=300)
sc.pl.pca_variance_ratio(dc_adata)

# %%
dc.get_metadata_associations(
    dc_adata,
    obs_keys = ['genotype', 'disease', 'apoe', 'treatment', 'genotype_treatment', 'psbulk_n_cells', 'psbulk_counts',],  # Metadata columns to associate to PCs
    obsm_key='X_pca',  # Where the PCs are stored
    uns_key='pca_anova',  # Where the results are stored
    inplace=True,
)

# %%
dc.plot_associations(
    dc_adata,
    uns_key='pca_anova',  # Summary statistics from the anova tests
    obsm_key='X_pca',  # where the PCs are stored
    stat_col='p_adj',  # Which summary statistic to plot
    obs_annotation_cols = ['disease', 'treatment', 'genotype_treatment', 'apoe', 'sample'], # which sample annotations to plot
    titles=['Principle component scores', 'Adjusted p-values from ANOVA'],
    figsize=(7, 5),
    n_factors=10,
)

# %%
sc.pl.pca(dc_adata, color='treatment', dimensions=(5,6), size=300)

# %%
print(dc_adata.obs['treatment'].tail())

# %%
# Build DESeq2 object
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=dc_adata,
    design_factors=['disease', 'treatment', 'apoe'],
    refit_cooks=True,
    inference=inference,
)



# %%
dds.deseq2()

# %%
dc_VEH_vs_LPS = DeseqStats(dds, contrast=['treatment', 'LPS', 'VEHICLE'])
dc_E3_vs_E4 = DeseqStats(dds, contrast=['apoe', 'E3', 'E4'])  

# %%
dc_VEH_vs_LPS.summary()


# %%
VEH_vs_LPS_results = dc_VEH_vs_LPS.results_df

# %%
dc.plot_volcano_df(
    VEH_vs_LPS_results,
    sign_thr=0.05, lFCs_thr=0.25,
    x='log2FoldChange',
    y='padj',
    top=50,
    figsize=(8, 4)
)


# %%
dc_E3_vs_E4.summary()

# %%
E3_vs_E4_results = dc_E3_vs_E4.results_df

# %%
dc.plot_volcano_df(
    E3_vs_E4_results,
    x='log2FoldChange',
    y='padj',
    top=20,
    figsize=(8, 4)
)

# %% [markdown]
# ## Transcription factor activity inference

# %%
# mat = VEH_vs_LPS_results[['stat']].T.rename(index={'stat': 'VEH_vs_LPS'})

# %%
# # Retrieve CollecTRI gene regulatory network
# collectri = dc.get_collectri(organism='mouse', split_complexes=False)
# collectri


# %%
# # Infer pathway activities with ulm
# tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=collectri)
# tf_acts

# %%
# dc.plot_barplot(
#     acts=tf_acts,
#     contrast='VEH_vs_LPS',
#     top=25,
#     vertical=True,
#     figsize=(3, 6)
# )

# %%
# Extract logFCs and pvals
logFCs = VEH_vs_LPS_results[['log2FoldChange']].T.rename(index={'log2FoldChange': 'VEH_vs_LPS'})
pvals = VEH_vs_LPS_results[['padj']].T.rename(index={'padj': 'VEH_vs_LPS'})

# # Plot
# dc.plot_volcano(
#     logFCs=logFCs,
#     pvals=pvals,
#     contrast='VEH_vs_LPS',
#     name='E2F4',
#     net=collectri,
#     top=10,
#     sign_thr=0.05,
#     lFCs_thr=0.5
# )

# %%
# dc.plot_network(
#     net=collectri,
#     obs=mat,
#     act=tf_acts,
#     n_sources=[],
#     n_targets=15,
#     node_size=100,
#     figsize=(7, 7),
#     c_pos_w='darkgreen',
#     c_neg_w='darkred',
#     vcenter=True
# )

# %% [markdown]
# ## DEGs/Volcanos\
# 
# Begin by subsettig data to isolate treatment variable. Do wilcoxon rank gene on each(4) subset to get DEGs for each genotype isolating effect of treatment. Plot volcano plots and get lists. find overlap genes between the different genotypes.

# %%
# Set the log normalized data as the main expression data in `adata`
adata.X = adata.layers['log_norm'].toarray().copy()  # Ensure this is dense if needed

# Now create subsets without modifying `.X` in each subset
adata_e3fad = adata[adata.obs['genotype'] == 'E3FAD'].copy()
adata_e4fad = adata[adata.obs['genotype'] == 'E4FAD'].copy()
adata_e3wt = adata[adata.obs['genotype'] == 'E3WT'].copy()
adata_e4wt = adata[adata.obs['genotype'] == 'E4WT'].copy()

adata_subs = [adata_e3fad, adata_e4fad, adata_e3wt, adata_e4wt]
print(adata_subs)
print(adata)

# %%

for adata_subset in adata_subs:
    sc.pp.pca(adata_subset)
    sc.pp.neighbors(adata_subset, n_pcs=25)
    sc.tl.umap(adata_subset)
    sc.pl.umap(adata_subset, color=['sample', 'genotype', 'apoe', 'treatment', 'pct_counts_mt'], ncols=2)

# %%


# %% [markdown]
# ## GO 


