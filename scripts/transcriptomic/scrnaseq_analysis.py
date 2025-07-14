# %% [markdown]
# # Setup

# %%
import re
import numpy as np
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from bbknn import bbknn
import venn
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches
from math import ceil
import seaborn as sns
import decoupler as dc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from gprofiler import GProfiler
import omnipath
import importlib
import gc
import os
import sys
import importlib
from pathlib import Path
import warnings


# add the project root to sys.path for imports
notebook_dir = Path.cwd().resolve()
project_root = notebook_dir.parents[1]
if str(project_root) not in sys.path:
    sys.path.append(str(project_root))

# import and reload helper functions
import func_lib as f
importlib.reload(f)

# ignore pandas performance warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# set random seed for reproducibility
seed = 2

# define where to save figures and tables
FIG_DIR = project_root / "results" / "transcriptomic" / "figures"
TABLE_DIR = project_root / "results" / "transcriptomic" / "tables"

for d in (FIG_DIR, TABLE_DIR):
    d.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Data Import

# %%

# Step 1: Get human MSigDB hallmark gene sets
msigdb_human = dc.get_resource('MSigDB', organism='human')

# Take a quick look
print(msigdb_human.head())

# Step 2: Extract unique human genes to map
unique_human_genes = msigdb_human['genesymbol'].unique().tolist()

# Step 3: Use g:Profiler to map from human to mouse
gp = GProfiler(return_dataframe=True)

# Perform orthology conversion from human to mouse using gene symbols
orthologs = gp.orth(
    organism='hsapiens',
    query=unique_human_genes,
    target='mmusculus'
)

# Check the results
print(orthologs.head())

# Step 4: Merge the orthology map with MSigDB
# Prepare map: original human symbol to mapped mouse symbol
ortholog_map = orthologs[['incoming', 'name']]
ortholog_map.columns = ['genesymbol', 'mouse_genesymbol']

# Merge with MSigDB hallmark gene sets
msigdb_mouse = msigdb_human.merge(ortholog_map, on='genesymbol')

# Drop 'genesymbol' column
msigdb_mouse = msigdb_mouse.drop(columns=['genesymbol'])

# Rename 'mouse_genesymbol' to 'genesymbol', keep all other columns
msigdb_mouse = msigdb_mouse.rename(columns={'mouse_genesymbol': 'genesymbol'})



# Remove duplicates
msigdb_mouse = msigdb_mouse.drop_duplicates()

# Final structure: source (pathway), genesymbol (mouse gene)
print(msigdb_mouse.head())


#subset reactome
msigdb_mouse_reactome = msigdb_mouse[msigdb_mouse['collection'] == 'reactome_pathways']

# ubset hallmark
msigdb_mouse_hallmark = msigdb_mouse[msigdb_mouse['collection'] == 'hallmark']

# subset kegg
msigdb_mouse_kegg = msigdb_mouse[msigdb_mouse['collection'] == 'kegg_pathways']

# %%
adatas = {
    "D1": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/collect_tube_1_batch_3_June_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D2": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/collect_tube_2_batch_3_june_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D3": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/Isol_Microglia_EFAD_TD_august_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D4": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/MICROGLIA_E3E4FAD_TD_23_05_2023_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D5": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/Tube_1_July_20_TD_YYRFC_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D6": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/Tube2_July_20_TD_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D7": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/Tube3_july_20_TD_outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
    "D8": sc.read_10x_mtx(
        "/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/Results_HT/Isol_microglia_EFAD_Sept_outs/outs/filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=False,
        make_unique=True,
        gex_only=False
    ),
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
manual_annotations = pd.read_csv('/Users/loganlarlham/Documents/BINP52/data/transcriptomic/raw/manual_annotation.csv', index_col=0)




# %%
len(treatment_map), len(genotype_map), len(manual_annotations)

# %%
print(f"Number of hashtags in D8: {adatas['D8'].var[adatas['D8'].var['feature_types'] == 'Antibody Capture'].shape[0]}")
for idx, row in adatas["D8"].var.iterrows():
    if row["feature_types"] == "Antibody Capture":
        print(f"gene_ids: {row['gene_ids']}")
        print(f"index: {idx}")
        # Convert the feature name to an integer index
        col_index = adatas['D8'].var_names.get_loc(idx)
        # Use the integer index to access the data
        count = adatas['D8'].X[:, col_index].sum()
        print(f"Number of {idx} in D8: {count}")
        print("\n")

# %% [markdown]
# # Hashtag Correction & Demultiplexing

# %%
# The index of the rows which are feature_types == antibody capture of D8 is incorrect, so we need to fix it



# Print the number of cells in D8
print(f"Number of vars in D8: {adatas['D8'].var.shape[0]}")

# Create a copy of the var DataFrame for modification
var_df = adatas['D8'].var.copy()

# Create a list to store indices of variables to remove
variables_to_remove = []

# Create a dictionary to map old variable names to new ones
rename_dict = {}

# Iterate over the var DataFrame
for index, row in var_df.iterrows():
    if row["feature_types"] == "Antibody Capture":
        if row['gene_ids'] == "Hashtag_1":
            # Rename 'Hashtag_1' to 'K283'
            rename_dict[index] = 'K283'
        elif row['gene_ids'] == "Hashtag_2":
            # Rename 'Hashtag_2' to 'K284'
            rename_dict[index] = 'K284'
        elif row['gene_ids'] == "Hashtag_3":
            # Add the index to the list of variables to remove
            variables_to_remove.append(index)
        elif row['gene_ids'] == "Hashtag_4":
            variables_to_remove.append(index)

# Apply the renaming to var
var_df.rename(index=rename_dict, inplace=True)

# Update the var DataFrame in the AnnData object
adatas['D8'].var = var_df

# Update var_names to match var.index
adatas['D8'].var_names = var_df.index

# Remove the variables from the AnnData object
adatas['D8'] = adatas['D8'][:, ~adatas['D8'].var.index.isin(variables_to_remove)]

print(f"Number of vars in D8 after removing hashtag_3/4: {adatas['D8'].var.shape[0]}")

# Verify the changes
for index, row in adatas["D8"].var.iterrows():
    if row["feature_types"] == "Antibody Capture":
        print(f"gene_ids: {row['gene_ids']}")
        print(f"index: {index}")
        print("\n")

# %% [markdown]
# # Demultiplex by HTO

# %%
for adata in adatas.values():
    hashtag_var_mask = adata.var['feature_types'] == 'Antibody Capture'
    hashtag_vars = adata.var_names[hashtag_var_mask]
    print(hashtag_vars)
    for hashtag_var in hashtag_vars:
        # .X is typically in CSR format, so convert to array for the single column
        # (take .A1 or .toarray()[:, 0] or similar)
        adata.obs[hashtag_var] = adata[:, hashtag_var].X.toarray().ravel()
        print(adata.obs.columns)
        print(adata.obs[hashtag_var])

# %%
#Demultiplex hashtags from datasets
for adata_key, adata in adatas.items():
    print(f"{adata_key}:")
    if adata_key == "D8":
        f.hash_demulitplex(adata, hashtag_vars, number_of_noise_barcodes=1)
    else:
        f.hash_demulitplex(adata, adata.obs.columns)
    print("\n")

# %% [markdown]
# # Cell Filtering & Annotation

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
# # Quality Control

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
# # Remove Non-Microglia (need to add manual annotation step from analysis1 notebook)

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
# # Merge & Subset

# %%
adata = ad.concat(list(QC_filtered_adatas.values()), axis=0, label="sample", join="outer", merge="unique")
adata.obs_names_make_unique()


print(adata)

# %%

# remove and gc the QC_filtered_adatas and adatas
del QC_filtered_adatas
del adatas
gc.collect()

# %%
# Group by 'genotype', 'treatment', and count unique 'classification' values
breakdown = adata.obs.groupby(['genotype', 'treatment'])['Classification'].nunique().reset_index()

# Rename columns for clarity
breakdown.columns = ['Genotype', 'Treatment', 'Number of individuals']

# Print the resulting table
print(breakdown)

# %%
 # Count the number of unique classifications
num_classifications = adata.obs['Classification'].nunique()
num_genotype_treatment = adata.obs['genotype_treatment'].nunique()

# Print the result
print(f"Number of unique classifications: {num_classifications}")
print(f"Number of unique genotype_treatment combinations: {num_genotype_treatment}")

# %% [markdown]
# # Make Subsets for Experimental Conditions

# %%
adata.obs

# %%
adata.uns['name'] = 'Full Dataset (All Samples, only Microglia, MT-regressed, QC-filtered)'

# Now create subsets for each genotype
adata_e3fad = adata[adata.obs['genotype'] == 'E3FAD'].copy()
adata_e3fad.uns['name'] = 'E3FAD'
adata_e4fad = adata[adata.obs['genotype'] == 'E4FAD'].copy()
adata_e4fad.uns['name'] = 'E4FAD'
adata_e3wt = adata[adata.obs['genotype'] == 'E3WT'].copy()
adata_e3wt.uns['name'] = 'E3WT'
adata_e4wt = adata[adata.obs['genotype'] == 'E4WT'].copy()
adata_e4wt.uns['name'] = 'E4WT'
# Now subset for each specific genotype and treatment
adata_e3fad_lps = adata_e3fad[adata_e3fad.obs['treatment'] == 'LPS'].copy()
adata_e3fad_lps.uns['name'] = 'E3FAD_LPS'
adata_e3fad_vehicle = adata_e3fad[adata_e3fad.obs['treatment'] == 'VEHICLE'].copy()
adata_e3fad_vehicle.uns['name'] = 'E3FAD_VEHICLE'
adata_e4fad_lps = adata_e4fad[adata_e4fad.obs['treatment'] == 'LPS'].copy()
adata_e4fad_lps.uns['name'] = 'E4FAD_LPS'
adata_e4fad_vehicle = adata_e4fad[adata_e4fad.obs['treatment'] == 'VEHICLE'].copy()
adata_e4fad_vehicle.uns['name'] = 'E4FAD_VEHICLE'
adata_e3wt_lps = adata_e3wt[adata_e3wt.obs['treatment'] == 'LPS'].copy()
adata_e3wt_lps.uns['name'] = 'E3WT_LPS'
adata_e3wt_vehicle = adata_e3wt[adata_e3wt.obs['treatment'] == 'VEHICLE'].copy()
adata_e3wt_vehicle.uns['name'] = 'E3WT_VEHICLE'
adata_e4wt_lps = adata_e4wt[adata_e4wt.obs['treatment'] == 'LPS'].copy()
adata_e4wt_lps.uns['name'] = 'E4WT_LPS'
adata_e4wt_vehicle = adata_e4wt[adata_e4wt.obs['treatment'] == 'VEHICLE'].copy()
adata_e4wt_vehicle.uns['name'] = 'E4WT_VEHICLE'
# Also creat subset for 'apoe' and 'treatment
adata_e4 = adata[adata.obs['apoe'] == 'E4'].copy()
adata_e4.uns['name'] = 'E4'
adata_e3 = adata[adata.obs['apoe'] == 'E3'].copy()
adata_e3.uns['name'] = 'E3'
adata_e4_lps = adata_e4[adata_e4.obs['treatment'] == 'LPS'].copy()
adata_e4_lps.uns['name'] = 'E4_LPS'
adata_e4_vehicle = adata_e4[adata_e4.obs['treatment'] == 'VEHICLE'].copy()
adata_e4_vehicle.uns['name'] = 'E4_VEHICLE'
adata_e3_lps = adata_e3[adata_e3.obs['treatment'] == 'LPS'].copy()
adata_e3_lps.uns['name'] = 'E3_LPS'
adata_e3_vehicle = adata_e3[adata_e3.obs['treatment'] == 'VEHICLE'].copy()
adata_e3_vehicle.uns['name'] = 'E3_VEHICLE'





adata_subs = [adata_e4, adata_e3, adata_e4_lps, adata_e4_vehicle, adata_e3_lps, adata_e3_vehicle, adata_e3fad, adata_e4fad, adata_e3wt, adata_e4wt, adata_e3fad_lps, adata_e3fad_vehicle, adata_e4fad_lps, adata_e4fad_vehicle, adata_e3wt_lps, adata_e3wt_vehicle, adata_e4wt_lps, adata_e4wt_vehicle]
print(adata_subs)
print(adata)

# %% [markdown]
# # Normalization & Feature Selection

# %%
# before normalizing, save the raw counts
adata.layers["counts"] = adata.X.copy()
f.plot_total_counts_vs_cells(adata, bins=250, zoom_x_range=(0, 5000), zoom_y_range=(0, 1000))
#  same for the subsets
for adata_sub in adata_subs:
    adata_sub.layers["counts"] = adata_sub.X.copy()

# %%
# normalize and log the data and store that in a different layer too incase we need to switch between them.
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['log_norm'] = adata.X.copy()
#so we can plot a histogram of the total counts/cells after norm
adata.obs["norm_total_counts"] = adata.X.sum(axis=1)    
f.plot_total_counts_vs_cells(adata, bins=250, zoom_x_range=(0, 5000), zoom_y_range=(0, 1000))

#  same for the subsets (without plot)      
for adata_sub in adata_subs:
    sc.pp.normalize_total(adata_sub)
    sc.pp.log1p(adata_sub)
    adata_sub.layers['log_norm'] = adata_sub.X.copy()

# %%
sc.pp.highly_variable_genes(adata, batch_key="sample")
sc.pl.highly_variable_genes(adata)
#make some variables to print a nice little summary
num_highly_variable_genes_t = adata.var['highly_variable'].sum()
total_genes_t = adata.var.shape[0]

print(f"{num_highly_variable_genes_t} out of {total_genes_t} genes are considered highly variable.")

#  same for the subsets
for adata_sub in adata_subs:
    sc.pp.highly_variable_genes(adata_sub, batch_key="sample")
    num_highly_variable_genes_t = adata.var['highly_variable'].sum()
    total_genes_t = adata.var.shape[0]

    print(f"{num_highly_variable_genes_t} out of {total_genes_t} genes are considered highly variable.")   

# %% [markdown]
# # Cell Cycle Scoring & Regression

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
# # Regress Mitochondrial genes

# %%
adata.X = adata.layers['log_norm'].copy()

#  same for the subsets
for adata_sub in adata_subs:
    adata_sub.X = adata_sub.layers['log_norm'].copy()

# %%

adata = sc.pp.regress_out(adata, 'pct_counts_mt', copy=True)

#  same for the subsets
for adata_sub in adata_subs:
    adata_sub = sc.pp.regress_out(adata_sub, 'pct_counts_mt', copy=True)
    


# %% [markdown]
# # Dimensionality Reduction

# %%
# Prior to various plots, Lets the levels of the categorical variable 'genotype_treatment' in the dataset so it is ordered correctly in the plots
adata.obs['genotype_treatment'] = pd.Categorical(adata.obs['genotype_treatment'], categories=['E3WT_VEHICLE', 'E3WT_LPS', 'E3FAD_VEHICLE', 'E3FAD_LPS', 'E4WT_VEHICLE', 'E4WT_LPS', 'E4FAD_VEHICLE', 'E4FAD_LPS'], ordered=True)

# %%
adata

# %%

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
plt.savefig(FIG_DIR / "pca_variance_ratio.pdf", dpi=300, bbox_inches='tight', transparent=True)
plt.close()

colors = [
    "sample", "sample", "pct_counts_mt", "pct_counts_mt", 
    "genotype", "genotype", "treatment", "treatment", 
    "apoe", "apoe", "disease", "disease"
]

dimensions = [
    (0, 1), (2, 3), (0, 1), (2, 3), 
    (0, 1), (2, 3), (0, 1), (2, 3), 
    (0, 1), (2, 3), (0, 1), (2, 3)
]

for i, (color, dims) in enumerate(zip(colors, dimensions)):
    fig, ax = plt.subplots()
    sc.pl.pca(
        adata,
        color=color,
        dimensions=[dims],
        size=2,
        show=False,
        ax=ax
    )
    fig.savefig(FIG_DIR / f"pca_{color}_{dims[0]}_{dims[1]}.pdf", dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)

#  same for the subsets (without plot)
for adata_sub in adata_subs:
    sc.tl.pca(adata_sub)

# %%
sc.pp.neighbors(adata, n_pcs=25)

sc.tl.umap(adata, min_dist=0.5)


# same for the subsets
for adata_sub in adata_subs:
    sc.pp.neighbors(adata_sub, n_pcs=25)
    sc.tl.umap(adata_sub)

# %%
sc.pl.umap(adata, color=['sample','genotype', 'apoe', 'treatment', 'pct_counts_mt'], size =10, ncols=2,  frameon=False, add_outline=True)

# %%
sc.tl.umap(adata, min_dist=0.5, n_components=3)

# %%
sc.pl.umap(adata, color=['genotype', 'apoe'], size=15,  add_outline=True, projection='3d')

# %% [markdown]
# # Clustering

# %%
i = 0.1
while i <= 1.9:
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

# %%


sc.pl.umap(adata, color=['leiden_0.3', 'leiden_0.5', 'leiden_0.7', 'leiden_0.9', 'leiden_1.1', 'leiden_1.3', 'leiden_1.5', 'leiden_1.7', 'leiden_1.9'], legend_loc="on data", ncols=2)



# %%
# plot with DEGs by cluster
f.plot_umap_with_gene_list(adata, 'leiden_1.9', n_top_genes=20)

# %% [markdown]
# # Cell-Type Assignment

# %%
def assign_cell_types2(adata):
    # Assign tentative cell type labels in adata.obs called ['anno2'] based on the top genes
    # Initialize the 'anno2' column with empty strings or any default value
    adata.obs['anno2'] = ''

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['0','4']) , 'anno2'] = 'Ribosomal with DAM markers'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['1','6','7']), 'anno2'] = 'DAM'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['5']), 'anno2'] = 'Cycling (G2M/Cdk1+)'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['8', '10','18']), 'anno2'] = 'IEG-enriched Homeostatic'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['2','3','9','15','16','17','22','23']), 'anno2'] = 'Homeostatic'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['11']), 'anno2'] = 'Cytokine Response'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['12']), 'anno2'] = 'MHC-II DAM'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['14', '19']), 'anno2'] = 'BAM-like'

    adata.obs.loc[adata.obs['leiden_1.9'].isin(['13', '20']), 'anno2'] = 'Ribosomal Biogenesis'
    
    adata.obs.loc[adata.obs['leiden_1.9'].isin(['21']), 'anno2'] = 'BAM-like'

    # Assign 'Unknown' to cells where 'anno2' is still empty
    adata.obs.loc[adata.obs['anno2'] == '', 'anno2'] = 'Unknown'

assign_cell_types2(adata)

# %%
adata.obs['anno2'].value_counts()

# %%
# make sure you’re using the log-normalized layer
adata.X = adata.layers['log_norm'].copy()

# get all the distinct labels
all_labels = adata.obs['anno2'].unique().tolist()

# loop: for each label as reference, compare all the others against it
for ref in all_labels:
    others = [lbl for lbl in all_labels if lbl != ref]
    print(f"\n=== Ranking genes: reference = {ref}; testing groups = {others} ===")

    # compute DE
    sc.tl.rank_genes_groups(
        adata,
        groupby='anno2',
        groups=others,
        reference=ref,
        method='wilcoxon',
        use_raw=False,
        rankby_abs=True  # if you want absolute log-fold changes
    )

    # plot the top 30 genes per tested group
    sc.pl.rank_genes_groups(
        adata,
        groups=others,
        n_genes=30,
        sharey=False,            # prevents y-axis autoscaling across panels
        title=f"vs {ref}",       # puts the reference in the figure title
        show=True
    )

# %%
tnfa_nfkb_genes = msigdb_mouse_hallmark[msigdb_mouse_hallmark['geneset'] == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB']['genesymbol'].tolist()


adata.X = adata.layers['log_norm'].copy()

# 1) Pairwise DE vs Cytokine Response for IEG-enriched Homeostatic & DAM
sc.tl.rank_genes_groups(
    adata,
    groupby='anno2',
    groups=['IEG-enriched Homeostatic', 'DAM'],
    reference='Cytokine Response',
    method='wilcoxon',
    use_raw=False,
    rankby_abs=True
)

sc.pl.rank_genes_groups(
    adata,
    groups=['DAM', 'IEG-enriched Homeostatic'],
    n_genes=30,
)

group_ids = ['DAM', 'IEG-enriched Homeostatic']
top_genes = {}

for group_id in group_ids:
    uns = adata.uns['rank_genes_groups']
    scores = uns['scores'][group_id]
    names = uns['names'][group_id]

    df = pd.DataFrame({'gene': names, 'score': scores.astype(float)})
    df_sorted = df.sort_values('score', ascending=True).head(100)
    print(f"Lowest score for {group_id}: {df_sorted['score'].max()}")

    top_genes[group_id] = set(df_sorted['gene'])

unique_to_dam = top_genes['DAM'] - top_genes['IEG-enriched Homeostatic']
unique_to_tim = top_genes['IEG-enriched Homeostatic'] - top_genes['DAM']
shared_genes = top_genes['DAM'] & top_genes['IEG-enriched Homeostatic']
shared_tnfa_nfkb = shared_genes & set(tnfa_nfkb_genes)

print(f"TNFA via NF-κB pathway genes: {sorted(tnfa_nfkb_genes)}\n")
print(f"Genes unique to DAM:\n{sorted(unique_to_dam)}\n")
print(f"Genes unique to IEG-enriched Homeostatic:\n{sorted(unique_to_tim)}\n")
print(f"Genes shared between DAM and IEG-enriched Homeostatic:\n{sorted(shared_genes)}\n")
print(f"Shared in TNF-NF-κB pathway:\n{sorted(shared_tnfa_nfkb)}\n")


# 2) “Vs Rest” DE for all groups
groups_of_interest = ['IEG-enriched Homeostatic', 'DAM', 'Cytokine Response']
sc.tl.rank_genes_groups(
    adata,
    groupby='anno2',
    method='wilcoxon',
    use_raw=False,
    rankby_abs=False,
    n_genes=adata.raw.shape[1] if adata.raw is not None else adata.shape[1]
)

sig_genes_by_group = {}
uns = adata.uns['rank_genes_groups']
for group in groups_of_interest:
    pvals_adj = np.array(uns['pvals_adj'][group])
    gene_names = np.array(uns['names'][group])
    scores = np.array(uns['scores'][group])

    up = set(gene_names[(pvals_adj < 0.05) & (scores > 0)])
    down = set(gene_names[(pvals_adj < 0.05) & (scores < 0)])
    sig_genes_by_group[group] = {'up': up, 'down': down}

tim_up, tim_down = sig_genes_by_group['IEG-enriched Homeostatic'].values()
dam_up, dam_down = sig_genes_by_group['DAM'].values()
cyt_up, cyt_down = sig_genes_by_group['Cytokine Response'].values()

print(f"UP shared by all: {len(tim_up & dam_up & cyt_up)}")
print(f"UP TIM & Cyt: {len(tim_up & cyt_up)}")
print(f"UP DAM & Cyt: {len(dam_up & cyt_up)}")
print(f"UP unique Cyt: {len(cyt_up - (tim_up | dam_up))}\n")

print(f"DOWN shared by all: {len(tim_down & dam_down & cyt_down)}")
print(f"DOWN TIM & Cyt: {len(tim_down & cyt_down)}")
print(f"DOWN DAM & Cyt: {len(dam_down & cyt_down)}")
print(f"DOWN unique Cyt: {len(cyt_down - (tim_down | dam_down))}\n")


# 3) Pairwise vs Cytokine Response (signed) and inline overlap comparison
sc.tl.rank_genes_groups(
    adata,
    groupby='anno2',
    groups=['IEG-enriched Homeostatic', 'DAM'],
    reference='Cytokine Response',
    method='wilcoxon',
    use_raw=False,
    rankby_abs=False
)

# extract pairwise DE
uns = adata.uns['rank_genes_groups']
pairwise_sig = {}
for group in ['IEG-enriched Homeostatic', 'DAM']:
    pvals_adj = np.array(uns['pvals_adj'][group])
    names = np.array(uns['names'][group])
    scores = np.array(uns['scores'][group])
    up = set(names[(pvals_adj < 0.05) & (scores > 0)])
    down = set(names[(pvals_adj < 0.05) & (scores < 0)])
    pairwise_sig[group] = {'up': up, 'down': down}

# inline overlap for IEG-enriched Homeostatic
tim_up_pair = pairwise_sig['IEG-enriched Homeostatic']['up']
tim_down_pair = pairwise_sig['IEG-enriched Homeostatic']['down']
print("--- IEG-enriched Homeostatic vs Cytokine Response ---")
print(f"UP overlap: {len(tim_up & tim_up_pair)} | rest only: {len(tim_up - tim_up_pair)} | pair only: {len(tim_up_pair - tim_up)}")
print(f"DOWN overlap: {len(tim_down & tim_down_pair)} | rest only: {len(tim_down - tim_down_pair)} | pair only: {len(tim_down_pair - tim_down)}\n")

# inline overlap for DAM
dam_up_pair = pairwise_sig['DAM']['up']
dam_down_pair = pairwise_sig['DAM']['down']
print("--- DAM vs Cytokine Response ---")
print(f"UP overlap: {len(dam_up & dam_up_pair)} | rest only: {len(dam_up - dam_up_pair)} | pair only: {len(dam_up_pair - dam_up)}")
print(f"DOWN overlap: {len(dam_down & dam_down_pair)} | rest only: {len(dam_down - dam_down_pair)} | pair only: {len(dam_down_pair - dam_down)}\n")


# 4) Optional TNF/NF-κB pathway hits
if 'tnfa_nfkb_genes' in globals():
    for group, label in zip(groups_of_interest, ['TIM','DAM','Cytokine']):
        up_hits = sig_genes_by_group[group]['up'] & set(tnfa_nfkb_genes)
        down_hits = sig_genes_by_group[group]['down'] & set(tnfa_nfkb_genes)
        print(f"{label}: {len(up_hits)} TNF/NF-κB UP, {len(down_hits)} DOWN")

# %%
genes = [
    'Ccl3', 'Ccl4',
    'Nfkbia', 'Nfkbiz',
    'Il1a', 'Il1b',
    'Tnf', 'Tnfaip3',
    'Atf3', 'Gadd45b',
    'Ccrl2', 'Bcl2a1b',
    'Tlr2', 'Zfp36',
    'Cd83', 'Csf1',
]

# Ensure that adata.X contains the log-normalized data
adata.X = adata.layers['log_norm'].copy()

for gene in genes:
    fig = sc.pl.umap(
        adata,
        color=[gene],
        size=50,
        use_raw=False,
        add_outline=False,
        vmax=5,
        return_fig=True
    )
    fig.savefig(str(FIG_DIR/f'umap_{gene}.pdf'), dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    sc.pl.violin(
        adata,
        [gene],
        groupby='anno2',
        use_raw=False,
        rotation=90,
        stripplot=False,
        ylabel=None,
        ax=ax,
        show=False
    )

    # Adjust x-axis tick label font size
    ax.tick_params(axis='x', labelsize=8)  # Adjust this number to your needs

    fig.savefig(str(FIG_DIR/f'violin_{gene}.pdf'), dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)

# %% [markdown]
# # Wilcoxon DEGs 

# %%

for adata_subset in adata_subs:
    adata_subset.X = adata_subset.layers['log_norm'].copy()
    sc.tl.rank_genes_groups(adata_subset, groupby='treatment', method='wilcoxon')

# %%
e3fad_LPS_sig_genes = f.plot_volcano(adata_e3fad, group='LPS', logfc_threshold=0.25, pval_threshold=0.05, top_n=200 )

print(f"Up-Regulated {len(e3fad_LPS_sig_genes['Up-regulated'])}: {e3fad_LPS_sig_genes['Up-regulated']}")
print(f"Down-Regulated {len(e3fad_LPS_sig_genes['Down-regulated'])} : {e3fad_LPS_sig_genes['Down-regulated']}")


# %%
e4fad_LPS_sig_genes = f.plot_volcano(adata_e4fad, group='LPS', logfc_threshold=0.1, pval_threshold=0.05, top_n=15 )

print(f"Up-Regulated {len(e4fad_LPS_sig_genes['Up-regulated'])} : {e4fad_LPS_sig_genes['Up-regulated']}")
print(f"Down-Regulated {len(e4fad_LPS_sig_genes['Down-regulated'])} : {e4fad_LPS_sig_genes['Down-regulated']}")

# %% [markdown]
# # Cell Type Visualizations

# %%
sc.pl.umap(adata, color='anno2')

# %%

sc.pl.umap(adata, color='anno2')

# %%
# Define the desired order of categories
desired_order = [
    'E3WT_VEHICLE', 'E4WT_VEHICLE', 'E3WT_LPS', 'E4WT_LPS', 
    'E3FAD_VEHICLE', 'E4FAD_VEHICLE', 'E3FAD_LPS', 'E4FAD_LPS'
]

# Create a contingency table (crosstab) of counts
counts = pd.crosstab(
    adata.obs['genotype_treatment'],
    adata.obs['anno2']
)

# Compute proportions by dividing each count by the row sum
proportions = counts.div(counts.sum(axis=1), axis=0)

# Get the tab20 colormap colors
tab20_colors = plt.get_cmap('tab20').colors
# Ensure you only use as many colors as there are cell types
num_cell_types = proportions.shape[1]
colors = tab20_colors[:num_cell_types]


# Enforce the order on the 'genotype_treatment' column
adata.obs['genotype_treatment'] = pd.Categorical(
    adata.obs['genotype_treatment'],
    categories=desired_order,
    ordered=True
)

# Get the categories directly from the Categorical object
categories = adata.obs['genotype_treatment'].cat.categories.tolist()
n_categories = len(categories)

# Determine the number of columns needed for the category UMAPs
n_category_cols = ceil(n_categories / 2)

# Set up the figure and gridspec with adjusted spacing
fig = plt.figure(figsize=(15, 10), constrained_layout=False)
width_ratios = [1]*n_category_cols + [2]  # General UMAP is wider and on the right
gs = fig.add_gridspec(nrows=2, ncols=n_category_cols + 1, width_ratios=width_ratios)

# Adjust the spacing between subplots
fig.subplots_adjust(hspace=0.01, wspace=0.3)  # Reduce hspace and wspace for better layout

# Create the general UMAP axes on the right
ax_general = fig.add_subplot(gs[:, -1])  # Spans both rows, last column

# Plot the general UMAP
assign_cell_types2(adata)
sc.pl.umap(
    adata,
    color='anno2',
    size=10,
    palette=colors,
    title='Clustered UMAP of\nMicroglia',  # Added newline
    ax=ax_general,
    show=False,  # Prevents the plot from displaying immediately
)
# Set aspect ratio to 'equal' to make the plot square
ax_general.set_aspect('equal')

# Optionally, set x and y limits to be the same for the general UMAP
xlim = ax_general.get_xlim()
ylim = ax_general.get_ylim()
min_limit = min(xlim[0], ylim[0])
max_limit = max(xlim[1], ylim[1])
ax_general.set_xlim(min_limit, max_limit)
ax_general.set_ylim(min_limit, max_limit)

# Adjust the title font size for the general UMAP
ax_general.set_title(ax_general.get_title(), fontsize=14)

# Now plot the category UMAPs without legends
for idx, cat in enumerate(categories):
    # Compute row and column positions
    row = idx % 2
    col = idx // 2  # Category plots occupy the first n_category_cols columns
    ax = fig.add_subplot(gs[row, col])

    # Plot UMAP for the current category
    assign_cell_types2(adata)
    sc.pl.umap(
        adata,
        color='anno2',
        palette=colors,
        groups=None,
        mask_obs=(adata.obs['genotype_treatment'] == cat),
        size=20,
        title=f'{cat} Microglia',  # Added newline
        ax=ax,
        show=False,   # Prevents individual plots from showing
        legend_loc='none'  # Removes the legend
    )
    # Set aspect ratio to 'equal' to make the plot square
    ax.set_aspect('equal')

    # Optionally, set x and y limits to be the same
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    min_limit = min(xlim[0], ylim[0])
    max_limit = max(xlim[1], ylim[1])
    ax.set_xlim(min_limit, max_limit)
    ax.set_ylim(min_limit, max_limit)

    # Adjust the title font size to prevent overlap
    ax.set_title(ax.get_title(), fontsize=12)

# # Save and show the figure
plt.savefig(str(FIG_DIR / "clustered_umap_all.pdf"),transparent=True,bbox_inches="tight")
plt.show()
plt.close()

# %%
# Density plots for each genotype_treatment group
sc.tl.embedding_density(adata, basis='umap', groupby='genotype_treatment', key_added='umap_density_genotype_treatment')
sc.pl.embedding_density(adata, basis='umap', key='umap_density_genotype_treatment', ncols=4, show=True, return_fig=True)

plt.savefig(str(FIG_DIR/"density_umap_genotype.pdf"), transparent=True, bbox_inches='tight')

# %%


# create and save stacked counts bar chart
counts = pd.crosstab(adata.obs['anno2'], adata.obs['phase'])
colors = plt.get_cmap('tab20').colors[:counts.shape[1]]
ax = counts.plot(
    kind='bar',
    stacked=True,
    color=colors,
    figsize=(12, 8)
)
plt.xlabel('Cell Group (anno2)')
plt.ylabel('Number of Cells')
plt.title('Cell Cycle Phase Counts per anno2 Group')
plt.legend(title='Cell Cycle Phase', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tight_layout()
plt.savefig(str(FIG_DIR / 'phase_counts_per_anno2.pdf'), transparent=True, bbox_inches='tight')
plt.show()
plt.close()
del counts

# create and save stacked proportions bar chart
counts = pd.crosstab(adata.obs['anno2'], adata.obs['phase'])
proportions = counts.div(counts.sum(axis=1), axis=0)
colors = plt.get_cmap('tab20').colors[:proportions.shape[1]]
ax = proportions.plot(
    kind='bar',
    stacked=True,
    color=colors,
    figsize=(12, 8)
)
plt.xlabel('Cell Group (anno2)')
plt.ylabel('Proportion of Cells')
plt.title('Cell Cycle Phase Proportions per anno2 Group')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
plt.legend(title='Cell Cycle Phase', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tight_layout()
plt.savefig(str(FIG_DIR / 'phase_proportions_per_anno2.pdf'), transparent=True, bbox_inches='tight')
plt.show()
plt.close()
del proportions, counts

# %%

assign_cell_types2(adata)

# build contingency tables
counts = pd.crosstab(adata.obs['genotype_treatment'], adata.obs['anno2'])
print(counts)
proportions = counts.div(counts.sum(axis=1), axis=0)
print(proportions)

# choose colors
tab20 = plt.get_cmap('tab20').colors
colors = tab20[:counts.shape[1]]

# plot and save raw counts
ax = counts.plot(kind='bar', stacked=True, color=colors, figsize=(12, 8))
plt.xlabel('Genotype Treatment')
plt.ylabel('Number of Cells')
plt.title('Cell Counts per Genotype Treatment and Cell Type')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tight_layout()
plt.savefig(str(FIG_DIR / 'stacked_bar_chart_counts.pdf'), transparent=True, bbox_inches='tight')
plt.show()
plt.close()

# plot and save proportions
ax = proportions.plot(kind='bar', stacked=True, color=colors, figsize=(12, 8))
plt.xlabel('Genotype Treatment')
plt.ylabel('Proportion of Cells')
plt.title('Cell Type Proportions per Genotype Treatment')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
plt.legend(title='Expression Profile', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tight_layout()
plt.savefig(str(FIG_DIR / 'stacked_bar_chart_proportions.pdf'), transparent=True, bbox_inches='tight')
plt.show()
plt.close()

del counts, proportions

# %%

# Create a mapping from cell types to safe column names
cell_types = adata.obs['anno2'].unique()

# Set Order of 'treatment' so that it appears in the correct order in the plot.
adata.obs['treatment'] = adata.obs['treatment'].cat.reorder_categories(['VEHICLE', 'LPS'])

# Set Order of 'apoe' so that it appears in the correct order in the plot.
adata.obs['apoe'] = adata.obs['apoe'].cat.reorder_categories(['E3', 'E4'])

# Function to create safe column names
def make_safe_column_name(name):
    return re.sub(r'\W+', '_', name)

# Create a mapping from original cell type names to safe column names
cell_type_map = {ct: make_safe_column_name(ct) for ct in cell_types}

# Make a dataframe with each sample as a row
anv_df = pd.DataFrame(index=adata.obs['Classification'].unique())

# Add the 'apoe', 'treatment', and 'disease' columns
anv_df['apoe'] = adata.obs.groupby('Classification')['apoe'].first().values
anv_df['treatment'] = adata.obs.groupby('Classification')['treatment'].first().values
anv_df['disease'] = adata.obs.groupby('Classification')['disease'].first().values

# Convert columns to category type
anv_df['apoe'] = anv_df['apoe'].astype('category')
anv_df['treatment'] = anv_df['treatment'].astype('category')
anv_df['disease'] = anv_df['disease'].astype('category')

# Add the percentage of each cell type in each sample
for ct in cell_types:
    ct_safe = cell_type_map[ct]
    anv_df[ct_safe] = adata.obs.groupby('Classification')['anno2'].apply(
        lambda x: (x == ct).sum() / len(x) * 100
    ).values

# Loop through each cell type to perform analysis and plotting
for ct in cell_types:
    ct_safe = cell_type_map[ct]
    print(f"Processing cell type: {ct}")

    # Prepare data for plotting and analysis by melting the dataframe
    df_long = anv_df.melt(
        id_vars=['apoe', 'treatment', 'disease'],
        var_name='cell_type',
        value_name='percent'
    )
    
    # Filter for the current cell type
    df_long = df_long[df_long['cell_type'] == ct_safe]
    
    # Check if there is data for this cell type
    if df_long.empty:
        print(f"No data for cell type {ct}.")
        continue  # Skip to the next cell type

    # Create a combined group that includes apoe, disease, and treatment.
    df_long['group'] = (df_long['apoe'].astype(str) + "_" +
                        df_long['disease'].astype(str) + "_" +
                        df_long['treatment'].astype(str))
        # Define the desired order for groups.
    # (Adjust these names to match your actual category names if needed)
    desired_order = [
        "E3_WT_VEHICLE", "E3_WT_LPS",
        "E3_FAD_VEHICLE", "E3_FAD_LPS",
        "E4_WT_VEHICLE", "E4_WT_LPS",
        "E4_FAD_VEHICLE", "E4_FAD_LPS"
    ]

    # Convert 'group' to a categorical variable with this order
    df_long['group'] = pd.Categorical(df_long['group'], categories=desired_order, ordered=True)
    df_long['apoe'] = df_long['apoe'].cat.reorder_categories(['E3', 'E4'], ordered=True)
    df_long['disease'] = df_long['disease'].cat.reorder_categories(['WT', 'FAD'], ordered=True)
    df_long['treatment'] = df_long['treatment'].cat.reorder_categories(['VEHICLE', 'LPS'], ordered=True)
    # Define custom x positions for each group.
    # Here, the first two (E3 WT) are placed close together, then a gap before E3 FAD, then a larger gap between E3 and E4.
    positions = {
        "E3_WT_VEHICLE": 0.0,
        "E3_WT_LPS": 0.4,
        "E3_FAD_VEHICLE": 1.2,
        "E3_FAD_LPS": 1.6,
        "E4_WT_VEHICLE": 3.0,
        "E4_WT_LPS": 3.4,
        "E4_FAD_VEHICLE": 4.2,
        "E4_FAD_LPS": 4.6,
    }
    # Map these positions into a new column
    df_long['x'] = df_long['group'].map(positions)
    # Define color palettes for E3 and E4.
    e3_wt_palette = sns.color_palette("Blues", 2)   # two shades for WT vs FAD in E3
    e3_fad_palette = sns.color_palette("Greens", 2)  # two shades for WT vs FAD in E3
    e4_wt_palette = sns.color_palette("Oranges", 2)  # two shades for WT vs FAD in E4
    e4_fad_palette = sns.color_palette("Purples", 2)    # two shades for WT vs FAD in E4


    # Map groups to colors (using the same color for the two treatments within each disease group)
    color_dict = {
        "E3_WT_VEHICLE": e3_wt_palette[0],
        "E3_WT_LPS": e3_wt_palette[1],
        "E3_FAD_VEHICLE": e3_fad_palette[0],
        "E3_FAD_LPS": e3_fad_palette[1],
        "E4_WT_VEHICLE": e4_wt_palette[0],
        "E4_WT_LPS": e4_wt_palette[1],
        "E4_FAD_VEHICLE": e4_fad_palette[0],
        "E4_FAD_LPS": e4_fad_palette[1],
    }

    # Compute group statistics (mean and standard deviation)
    group_stats = df_long.groupby('group')['percent'].agg(['mean', 'std']).reindex(desired_order)
    # Add x positions from our mapping
    group_stats['x'] = [positions[g] for g in desired_order]
            
    
    # Example: Perform a Three-Way ANOVA (apoe, disease, treatment and interactions)
    formula = ('percent ~ C(apoe) + C(disease) + C(treatment) + '
               'C(apoe):C(disease) + C(apoe):C(treatment) + C(disease):C(treatment) + '
               'C(apoe):C(disease):C(treatment)')
    model = ols(formula, data=df_long).fit()
    # ANOVA table
    anova_table = sm.stats.anova_lm(model, typ=2)
    print(f"Three-Way ANOVA Results for cell type {ct}:")
    print(anova_table)

    # Print full model summary
    print("\nOLS Model Summary:")
    print(model.summary())

    # Perform Tukey's HSD test using the new group definitions
    tukey = pairwise_tukeyhsd(
        endog=df_long['percent'],
        groups=df_long['group'],
        alpha=0.05
    )
    print("All group levels in order:", tukey.groupsunique)
    print("Baseline (first level):", tukey.groupsunique[0])
    print(f"\nTukey HSD Post-Hoc Test Results for cell type {ct}:")
    print(tukey.summary())

    # Optional: Visualize the Tukey HSD test results
    tukey.plot_simultaneous(comparison_name=df_long['group'].unique()[0])
    plt.title(f'Tukey HSD Post-Hoc Test for {ct}')
    plt.xlabel('Mean Difference')
    plt.show()


    # ---------------------------------------------------------------------
    # 1) grab the stats out of the summary table
    res = pd.DataFrame(
        data=tukey._results_table.data[1:], 
        columns=tukey._results_table.data[0]
    )
    # make sure numeric columns are floats
    for col in ['meandiff','lower','upper']:
        res[col] = res[col].astype(float)

    # 2) build a “pair” label
    res['pair'] = res['group1'] + ' – ' + res['group2']

    # 3) sort however you like (here by mean difference)
    res = res.sort_values('meandiff')

    # 4) manual error‐bar plot
    fig, ax = plt.subplots(figsize=(6, len(res) * 0.3))
    y = np.arange(len(res))
    x = res['meandiff']
    err_low  = x - res['lower']
    err_high = res['upper'] - x

    ax.errorbar(x, y, 
                xerr=[err_low, err_high], 
                fmt='o', capsize=4)
    ax.axvline(0, color='gray', linestyle='--')

    ax.set_yticks(y)
    ax.set_yticklabels(res['pair'])
    ax.set_xlabel('Mean difference')
    ax.set_title(f'All pairwise Tukey HSD comparisons for {ct}')

    plt.tight_layout()
    plt.show()


    # Plot using matplotlib's bar
    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    # Set a smaller bar width for “skinnier” bars
    bar_width = 0.35

    # Plot the bars with error bars
    ax.bar(
        group_stats['x'],
        group_stats['mean'],
        yerr=group_stats['std'],
        width=bar_width,
        color=[color_dict[g] for g in desired_order],
        capsize=5,
        edgecolor='black'
    )

    # Overlay individual sample points (as diamonds) with slight x jitter
    for g in desired_order:
        subset = df_long[df_long['group'] == g]
        x_pos = positions[g]
        # Create a small random jitter around the x position
        jitter = np.random.uniform(-0.05, 0.05, size=len(subset))
        ax.plot(
            x_pos + jitter,
            subset['percent'],
            'D',
            markersize=7,
            markeredgecolor='black',
            markerfacecolor=color_dict[g],  # Use same color as the corresponding bar
            linestyle='none'
        )


    # Set custom x-axis ticks and labels
    ax.set_xticks([positions[g] for g in desired_order])
    # Replace underscores with dashes for prettier labels (e.g. "E3-WT-VEHICLE")
    ax.set_xticklabels([g.replace("_", "-") for g in desired_order], rotation=45, fontsize=16)

    # Add y-axis padding
    ymin, ymax = ax.get_ylim()
    padding = (ymax - ymin) * 0.1  # 10% padding
    ax.set_ylim(ymin, ymax + padding)

    plt.title(f'Percentage of {ct} Cells Across Groups', fontsize=24, pad=40)
    plt.ylabel(f'Percentage of Cells in \n {ct}', fontsize=16)
    plt.xlabel('Group (apoe-disease-treatment)', fontsize=14)


    # Transform {ct} to a safe name for saving, removing spaces\
    safe_ct = re.sub(r'\W+', '_', ct)
    # Remove brackets
    safe_ct = re.sub(r'[()\+/\\]', '', safe_ct)

    plt.savefig(str(FIG_DIR / f"{safe_ct}_cells_across_groups.pdf"),transparent=True, bbox_inches='tight')
    plt.show()
    plt.close()

# %%

# Create a list of patches for the legend
patches = [
    mpatches.Patch(color=color_dict[g], label=g.replace("_", "-"))
    for g in desired_order
]

# Create an empty figure with no axes visible
fig, ax = plt.subplots(figsize=(4, 3))

# Create the legend in the center. You can adjust loc and number of columns as needed.
legend = ax.legend(handles=patches,
                    loc='center', 
                    frameon=False,   # Remove box if desired
                    ncol=1,          # Number of columns in the legend
                    title='Condition', 
                    title_fontsize=16,
                    fontsize=14)

# Remove the axes since you just want the legend
ax.axis('off')

# Use tight layout for a neat bounding box
plt.tight_layout()

# Save the legend as a separate file
plt.savefig(str(FIG_DIR/"bar_dot_legend.pdf"), transparent=True, bbox_inches='tight')
plt.show()

# %%



# Step 1: Subset to groups of interest
groups_of_interest = ['IEG-enriched Homeostatic', 'Cytokine Response', 'DAM']
adata_subset = adata[adata.obs['anno2'].isin(groups_of_interest)].copy()
adata_subset.obs['anno2'] = adata_subset.obs['anno2'].astype('category')
adata_subset.obs['anno2'] = adata_subset.obs['anno2'].cat.reorder_categories(groups_of_interest, ordered=True)

# Step 2: Define marker groups
homeostatic_genes = ['Tmem119', 'Cx3cr1', 'P2ry12', 'Olfml3', 'Sall1']
DAM_genes = ['Apoe', 'Trem2', 'Tyrobp', 'Cst7', 'Ctsd', 'Clec7a', 'Lpl', 'Cd63', 'Cd9', 'B2m', 'Ctsb', 'Fabp5', 'Fth1']
IEG_genes = ['Fos', 'Jun', 'Egr1', 'Dusp1', 'Jund', 'Junb', 'Fosb', 'Cebpb', 'Ier2', 'Ier3']
NFKB_canon_genes = ['Nfkbia', 'Nfkbiz', 'Ccl3', 'Ccl4', 'Il1a', 'Il1b', 'Tnf', 'Tnfaip3', 'Atf3', 'Gadd45b', 'Ccrl2', 'Bcl2a1b', 'Tlr2']

marker_dict = {
    "Homeostatic": homeostatic_genes,
    "IEG": IEG_genes,
    "DAM": DAM_genes,
    "NF-κB Pathway": NFKB_canon_genes
}

# Flatten list of all genes
all_genes = [gene for group in marker_dict.values() for gene in group]
all_genes = [gene for gene in all_genes if gene in adata_subset.var_names]  # Ensure all genes exist

# Step 3: Compute scaled expression (z-score) per gene per cell
adata_scaled = sc.pp.scale(adata_subset, copy=True)

# Step 4: Create base dotplot using original data (for size)
dotplot = sc.pl.DotPlot(
    adata_subset,
    var_names=marker_dict,
    groupby='anno2',
    use_raw=False,
    vmin=-2,
    vmax=2,
    var_group_rotation=0,
)

# Step 5: Inject scaled expression for color (adata_scaled)
dot_color_df = pd.DataFrame(index=dotplot.dot_color_df.index, columns=dotplot.dot_color_df.columns)

for gene in all_genes:
    gene_idx = adata_scaled.var_names.get_loc(gene)
    scaled_values = adata_scaled.X[:, gene_idx]

    for group in groups_of_interest:
        cells = adata_scaled.obs['anno2'] == group
        mean_scaled = scaled_values[cells].mean()
        dot_color_df.loc[group, gene] = mean_scaled

# Cast to float
dot_color_df = dot_color_df.astype(float)

# Step 6: Apply custom color data to the plot
dotplot.dot_color_df = dot_color_df

dotplot = dotplot.style(
    dot_edge_color='black',
    cmap='bwr',
).legend(
    colorbar_title="Z-score Expression").savefig(str(FIG_DIR/
    "dotplot_3_groups.pdf"),
    dpi=300,
    bbox_inches='tight',
    transparent=True,)

del adata_subset, adata_scaled

# %%
adata.X = adata.layers['log_norm'].copy()
sc.tl.rank_genes_groups(adata, groupby='anno2', method='wilcoxon', rankby_abs=False, pts=True)

adata.X = adata.layers['log_norm'].copy()
adata_scale = sc.pp.scale(adata, copy=True)

fig = sc.pl.rank_genes_groups_matrixplot(
    adata_scale,
    n_genes=10,
    groupby='anno2',
    key='rank_genes_groups',
    min_logfoldchange=0.9,
    vmax=3,
    vmin=-3,
    cmap='bwr',
    show=False,
    return_fig=True,
    colorbar_title='Z-score',
)

# First, trigger figure creation
fig.make_figure()

# Now safely adjust fonts
axes = fig.get_axes()
axes['mainplot_ax'].tick_params(axis='x', labelsize=18)
axes['mainplot_ax'].tick_params(axis='y', labelsize=18)

for text in axes['group_extra_ax'].texts:
    text.set_fontsize(24)

if 'color_legend_ax' in axes:
    axes['color_legend_ax'].tick_params(labelsize=14)

    plt.savefig(str(FIG_DIR / "rank_genes_groups_matrixplot.pdf"), transparent=True, bbox_inches="tight")

# %%
# ## Updated pseudobulk code, no subset individual contrasts 

# %%

# results containers
vehlps_stats = {}
dc_vehlps_results = {}
dc_vehlps_results_sig = {}
dc_vehlps_upreg = {}
dc_vehlps_downreg = {}



# 1) build one pseudobulk across all cells
adata.X = adata.layers['counts'].toarray()
psb = dc.get_pseudobulk(
    adata,
    sample_col='Classification',
    groups_col='genotype',
    mode='sum',
    min_cells=0,
    min_counts=0
)
psb.layers['counts'] = psb.X.copy()

# 2) normalize, log, scale, PCA
sc.pp.normalize_total(psb, target_sum=1e4)
sc.pp.log1p(psb)
sc.pp.scale(psb, max_value=10)
sc.tl.pca(psb)
dc.swap_layer(psb, 'counts', inplace=True)

# 3) ensure combined factor
psb.obs['genotype-treatment'] = (
    psb.obs['genotype'] + '-' + psb.obs['treatment']
).astype('category')


# Check PCA
sc.pl.pca_variance_ratio(psb, show = False)
plt.savefig(str(FIG_DIR / f"Pseudobulk_pca_scree.pdf"), dpi=300, bbox_inches='tight', transparent=True)
sc.pl.pca(psb, color=['disease', 'genotype', 'apoe', 'treatment', 'genotype-treatment', 'sample'], ncols=3, size=300, show = False)
plt.savefig(str(FIG_DIR / f"Pseudobulk_pca.pdf"), dpi=300, bbox_inches='tight', transparent=True)

# Check metadata associations
dc.get_metadata_associations(
    psb,
    obs_keys = ['apoe','disease','treatment', 'psbulk_n_cells', 'psbulk_counts',],  # Metadata columns to associate to PCs
    obsm_key='X_pca',  # Where the PCs are stored
    uns_key='pca_anova',  # Where the results are stored
    inplace=True,
)

dc.plot_associations(
    psb,
    uns_key='pca_anova',  # Summary statistics from the anova tests
    obsm_key='X_pca',  # where the PCs are stored
    stat_col='p_adj',  # Which summary statistic to plot
    obs_annotation_cols = ['apoe','disease','treatment'], # which sample annotations to plot,
    figsize=(7, 5),
)
plt.savefig(str(FIG_DIR / f"Pseudobulk_metadata"), dpi=300, bbox_inches='tight', transparent=True)


# 4) loop genotypes
for genotype in ('E3WT','E3FAD','E4WT','E4FAD'):
    key = f'{genotype}_vehlps'
    contrast = [ 'genotype-treatment',
                 f'{genotype}-LPS',
                 f'{genotype}-VEHICLE' ]

    # DESeq2
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        adata=psb,
        design_factors=['genotype-treatment'],
        refit_cooks=True,
        inference=inference
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=contrast)
    stats.summary()

    # store
    vehlps_stats[key]         = stats
    df                        = stats.results_df
    dc_vehlps_results[key]    = df
    sig                       = df.loc[
        (df['padj'] < 0.05) &
        (df['log2FoldChange'].abs() > 0.25)
    ].sort_values(by='log2FoldChange', key=abs, ascending=False)
    dc_vehlps_results_sig[key] = sig
    dc_vehlps_upreg[key]      = set(sig.loc[sig['log2FoldChange'] > 0].index)
    dc_vehlps_downreg[key]    = set(sig.loc[sig['log2FoldChange'] < 0].index)

    # print lists
    print(f"{genotype} up ({len(dc_vehlps_upreg[key])}):",
          sorted(dc_vehlps_upreg[key]))
    print(f"{genotype} down ({len(dc_vehlps_downreg[key])}):",
          sorted(dc_vehlps_downreg[key]))


    #My volcano
    fig = f.custom_volcano_plot(
        df=dc_vehlps_results[key],
        logfc_col='log2FoldChange',
        pval_col='padj',
        lfc_thresh=0.25,
        pval_thresh=0.05,
        top_n=15,
        x_lim=(-6,6),
        y_lim=(0,5),
        title=f"{genotype}: LPS vs VEHICLE",
        use_adjust_text=True
    )

    fig.savefig(str(FIG_DIR / f"{genotype}_volcano.pdf"), dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)



from matplotlib.colors import to_rgba

# 5) Venn diagrams
up_sets   = [dc_vehlps_upreg[f'{g}_vehlps']   for g in ('E3WT','E3FAD','E4WT','E4FAD')]
down_sets = [dc_vehlps_downreg[f'{g}_vehlps'] for g in ('E3WT','E3FAD','E4WT','E4FAD')]

# Mapping from base labels to full condition keys
group_keys = {
    'E3WT': 'E3_WT',
    'E3FAD': 'E3_FAD',
    'E4WT': 'E4_WT',
    'E4FAD': 'E4_FAD'
}

# Save original generate_colors function
original_generate_colors = venn.generate_colors

for name, sets, suffix in (
    ('Up', up_sets, '_LPS'),
    ('Down', down_sets, '_VEHICLE')
):
    # Extract RGBA colors with alpha=0.4
    colors = [to_rgba(color_dict[group_keys[g] + suffix], alpha=0.4) for g in ('E3WT','E3FAD','E4WT','E4FAD')]

    # Monkey-patch generate_colors to return these
    venn.generate_colors = lambda *args, **kwargs: colors

    # Create dataset dictionary
    label_dict = dict(zip(['E3WT','E3FAD','E4WT','E4FAD'], sets))

    # Generate and show Venn diagram
    fig = venn.venn(label_dict)
    plt.title(f"{name}-regulated Genes in LPS vs Vehicle by Genotype")
    plt.savefig(str(FIG_DIR / f"{name}_regulated_genes_venn.pdf"), dpi=300, bbox_inches='tight', transparent=True)
    plt.show()

# Restore original generate_colors function
venn.generate_colors = original_generate_colors

# 6) pairwise intersections printout
pairs = [
    ('E3FAD','E4FAD'),
    ('E3WT','E4WT'),
    ('E3FAD','E3WT'),
    ('E4FAD','E4WT'),
]
for name, col in (('Up','up_sets'), ('Down','down_sets')):
    print(f"\n{name}-regulated intersections:")
    for a,b in pairs:
        s = locals()[col][['E3WT','E3FAD','E4WT','E4FAD'].index(a)].intersection(
            locals()[col][['E3WT','E3FAD','E4WT','E4FAD'].index(b)]
        )
        print(f" {a} ∩ {b}: {s}")

# %% [markdown]
# ## Bivariate DEG plot

# %%
fig = f.bivariate_quadrant_plot(
    dc_vehlps_results['E3FAD_vehlps'],
    dc_vehlps_results['E4FAD_vehlps'],
    x_label='E3FAD LPS vs Vehicle',
    y_label='E4FAD LPS vs Vehicle',   
    use_adjust_text=True,
    label_genes=NFKB_canon_genes,
    x_lim=(-2, 2),
    y_lim=(-4, 4),
    dpi=300,
    title="E3FAD vs E4FAD LPS Responses"
)

fig.savefig(str(FIG_DIR / "E3FAD_vs_E4FAD_LPS_Responses.pdf"),dpi=300,bbox_inches="tight",transparent=True)
plt.close(fig)

# %% [markdown]
# # Comparison Of Wilcoxon vs Pseudobulk DEGs

# %%
# Format full gene list for display
def format_gene_list(gene_set):
    gene_list = sorted(gene_set)
    return 'None' if len(gene_list) == 0 else ', '.join(gene_list)

# Helper to plot Venn diagrams + full gene lists
def plot_venn_with_genes(sets_dict, title):
    sets = list(sets_dict.values())
    names = list(sets_dict.keys())
    
    labels = venn.get_labels(sets, fill=['number'])
    fig, ax = venn.venn2(labels, names=names)
    
    set1, set2 = sets
    unique_set1 = set1 - set2
    unique_set2 = set2 - set1
    intersection = set1 & set2

    text = (
        f"Unique to {names[0]} ({len(unique_set1)} genes):\n{format_gene_list(unique_set1)}\n\n"
        f"Unique to {names[1]} ({len(unique_set2)} genes):\n{format_gene_list(unique_set2)}\n\n"
        f"Intersection ({len(intersection)} genes):\n{format_gene_list(intersection)}"
    )
    
    plt.subplots_adjust(bottom=0.4)  # More room for big lists
    plt.figtext(0.5, 0.02, text, wrap=True, ha='center', fontsize=9)
    plt.title(title)
    plt.show()

# Upregulated E3FAD
e3fad_up = {
    "Pseudobulk": dc_vehlps_upreg["E3FAD_vehlps"],
    "Wilcoxon": set(e3fad_LPS_sig_genes['Up-regulated'])
}
plot_venn_with_genes(e3fad_up, title="E3FAD Upregulated Genes Venn Diagram")

# Upregulated E4FAD
e4fad_up = {
    "Pseudobulk": dc_vehlps_upreg["E4FAD_vehlps"],
    "Wilcoxon": set(e4fad_LPS_sig_genes['Up-regulated'])
}
plot_venn_with_genes(e4fad_up, title="E4FAD Upregulated Genes Venn Diagram")

# Downregulated E3FAD
e3fad_down = {
    "Pseudobulk": dc_vehlps_downreg["E3FAD_vehlps"],
    "Wilcoxon": set(e3fad_LPS_sig_genes['Down-regulated'])
}
plot_venn_with_genes(e3fad_down, title="E3FAD Downregulated Genes Venn Diagram")

# Downregulated E4FAD
e4fad_down = {
    "Pseudobulk": dc_vehlps_downreg["E4FAD_vehlps"],
    "Wilcoxon": set(e4fad_LPS_sig_genes['Down-regulated'])
}
plot_venn_with_genes(e4fad_down, title="E4FAD Downregulated Genes Venn Diagram")

# %% [markdown]
# # Functional & TF Analysis

# %% [markdown]
# ## Decoupler GSEA on Hallmark
# 

# %%
gsea_results = {}

for comp in dc_vehlps_results:
    print(f"Running GSEA for {comp}...")
    gsea_results[comp] = dc.get_gsea_df(
        dc_vehlps_results[comp],
        stat='stat',
        net=msigdb_mouse_hallmark,
        source='geneset',
        target='genesymbol',
        times=1000,
        min_n=5,
        seed=42,
        verbose=True
    )

# %%
for comp, gsea_df in gsea_results.items():
    print(f"\nGSEA results for {comp} (FDR p-value < 0.05):")
    sig_gsea_results = gsea_df[gsea_df["FDR p-value"] < 0.05]
    sig_gsea_results = sig_gsea_results.sort_values(by="NES", ascending=True)
    print(sig_gsea_results.head(25))
    # Plot
    genotype = comp.split('_')[0]            # e.g. 'E3WT'
    prefix = genotype[:2] + '_' + genotype[2:]  # e.g. 'E3_WT'

    vehicle_key = f"{prefix}_VEHICLE"  # e.g. 'E3_WT_VEHICLE'
    lps_key     = f"{prefix}_LPS"      # e.g. 'E3_WT_LPS'

    # 2. build a list of colors for each NES value
    nes = sig_gsea_results["NES"]
    bar_colors = [
        color_dict[vehicle_key] if v < 0 else color_dict[lps_key]
        for v in nes
    ]

    # 3. plot with those colors
    plt.figure(figsize=(10, 6))
    plt.barh(
        sig_gsea_results["Term"],
        nes,
        color=bar_colors
    )
    plt.axvline(x=0, color='black', linewidth=1)
    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.title(f"Significant Hallmark Pathways for {genotype} LPS vs Vehicle (FDR < 0.05)")
    plt.tight_layout()
    plt.savefig(str(FIG_DIR / f"{genotype}_gsea_hallmark.pdf"), dpi=300, bbox_inches="tight", transparent=True)
    plt.show()

# %% [markdown]
# ## Transcription Factor Inference with Decoupler

# %%

# Retrieve the CollecTRI gene regulatory network once
collectri = dc.get_collectri(organism='mouse', split_complexes=False)

# Dictionary to store results from each comparison
ulm_results = {}

# Loop over each dataset in your dictionary
for comparison_name, df_res in dc_vehlps_results.items():
    print(f"Processing comparison: {comparison_name}")

    # Prepare `mat` for decoupler
    mat = df_res[['stat']].T.rename(index={'stat': comparison_name})

    # Run ULM
    tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=collectri)
    
    ulm_results[comparison_name] = {
        'tf_acts': tf_acts,
        'tf_pvals': tf_pvals
    }

    # Extract the activity and p-value rows as Series
    tf_acts_series = tf_acts.loc[comparison_name]
    tf_pvals_series = tf_pvals.loc[comparison_name]

    # Get the top 12 TFs by activity
    top12_tfs = tf_acts_series.sort_values(ascending=False).head(12)

    # Get the corresponding p-values
    top12_pvals = tf_pvals_series[top12_tfs.index]

    # Combine into a single DataFrame for clean printing
    top12_df = pd.DataFrame({
        'activity': top12_tfs,
        'p-value': top12_pvals
    })

    print(f"\nTop 12 TFs by activity for {comparison_name}:")
    print(top12_df)


    # Parse comparison name
    parts = comparison_name.split('_')
    genotype = parts[0]
    contrast = parts[1]
    apoe_val = genotype[0:2]
    status_val = genotype[2:]
    
    barplot_title = f"Top differentially activated TFs in LPS vs Vehicle-treated {apoe_val} {status_val}"
    network_title = f"Network plot of top differentially activated TF in LPS vs Vehicle-treated {apoe_val} {status_val}"
    suptitle = f"TF inference of {genotype} LPS vs Vehicle"


    # === BARPLOT ===
    fig_bar, ax_bar = plt.subplots(figsize=(7, 6))
    dc.plot_barplot(
        acts=tf_acts,
        contrast=comparison_name,
        top=12,
        vertical=True,
        ax=ax_bar,
        cmap='bwr'
    )
    ax_bar.set_xlabel("Transcription Factor")
    ax_bar.set_ylabel("Inferred Activity")
    ax_bar.set_title(barplot_title)
    fig_bar.suptitle(suptitle, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save and show barplot
    fig_bar.savefig(str(FIG_DIR / f"ULM_{comparison_name}_barplot.pdf"), dpi=300, transparent=True)
    plt.show()
    plt.close(fig_bar)

    # === NETWORK PLOT ===
    fig_net = dc.plot_network(
        net=collectri,
        obs=mat,
        act=tf_acts,
        n_sources=6,
        n_targets=15,
        node_size=50,
        figsize=(7, 7),
        c_pos_w='darkgreen',
        c_neg_w='darkred',
        t_cmap='PiYG',
        vcenter=True,
        return_fig=True
    )
    fig_net.suptitle(network_title, fontsize=14)

    # Save and show network plot
    fig_net.savefig(str(FIG_DIR / f"ULM_{comparison_name}_network.pdf"), dpi=300, transparent=True)
    plt.show()
    plt.close(fig_net)

# %% [markdown]
# # Save Results

# %%
# differential expression tables
for name, df in dc_vehlps_results.items():
    df.to_csv(str(TABLE_DIR / f"dc_vehlps_{name}.csv"), index=True)

for name, df in dc_vehlps_results_sig.items():
    df.to_csv(str(TABLE_DIR / f"dc_vehlps_sig_{name}.csv"), index=True)

# GSEA result tables
for name, df in gsea_results.items():
    df.to_csv(str(TABLE_DIR / f"gsea_{name}.csv"), index=True)

# build a minimal AnnData with filtered obs/uns
keys_to_remove = (
    'Antigen-presenting response (HLA)',
    'Cytokines response 1 (CRM-1)',
    'Cytokines response 2 (CRM-2)',
    'Disease associated (DAM)',
    'Homeostatic (HM)',
    'Interferon response (IRM)',
    'Ribosomal response (RM)',
    'Transitioning CRM'
)

obs_keys = [
    k for k in adata.obs.keys()
    if not re.match(r'leiden_|K\d+|G\d+', k)
    and k not in keys_to_remove
]

uns_keys = [
    k for k in adata.uns.keys()
    if not re.match(r'leiden_', k)
]

adata_minimal = adata.copy()
adata_minimal.obs = adata.obs[obs_keys].copy()
adata_minimal.uns = {k: adata.uns[k] for k in uns_keys}

# save the minimal object into your tables directory (ignored by git)
adata_minimal.write_h5ad(str(TABLE_DIR / "Logan_adata.h5ad"))

# clean up
del adata_minimal, obs_keys, uns_keys, keys_to_remove


