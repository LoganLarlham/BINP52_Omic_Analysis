# ==============================================================================  
#                   __                                                          
#                 /\ \__                                                        
# _____   _ __  ___\ \ ,_\    __    ___    ___ ___      __                      
#/\ '__`\/\`'_\/ __`\ \ \/  /'__`\ / __`\/' __` __`\  /'__`\                    
#\ \ \L\ \ \ \/\ \L\ \ \ \_/\  __//\ \L\ /\ \/\ \/\ \/\  __/                    
# \ \ ,__/\ \_\ \____/\ \__\ \____\ \____\ \_\ \_\ \_\ \____\                   
#  \ \ \/  \/_/\/___/  \/__/\/____/\/___/ \/_/\/_/\/_/\/____/                   
#   \ \_\ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   
#    \/_/                                                                        
# ==============================================================================
#        PROTEOMIC DIFFERENTIAL EXPRESSION ANALYSIS        
# ==============================================================================
# Author   : Logan Larlham
# Date     : April 2nd, 2025
# Purpose  : In this script we will be analyzing the LS-MS/MS
#            data from our APOE-LPS project, focusing on FAD female samples.
# ==============================================================================
# Data Description: ----
# Sheet 1: Sample annotation (AnimalID, raw data filename, APOE genotype,
#          disease status, sex, treatment).
# Sheet 2: ProteinGroup table with gene names and raw log₂ LFQ intensities
#          per sample.
# Sheet 3: ProteinGroup table with normalized log₂ LFQ intensities
#          (used for all downstream analyses).
# ==============================================================================
# Analysis Steps:
# 1. Load & preprocess data: filter out “N/A” genes, pivot to long format,
#    build a SummarizedExperiment with expression matrix and metadata.
# 2. save & review normalization plot.
# 3. Assess missingness: generate frequency, count, missing-value heatmap,
#    and detection plots.
# 4. Impute missing values with MinProb (q = 0.5); produce imputation diagnostic.
# 5. Compute sample-to-sample Pearson correlations; plot full and top-500
#    variance heatmaps.
# 6. Subset to FAD female samples; relevel factors for Disease, Treatment,
#    Genotype, Sex.
# 7. Fit limma model (Genotype × Treatment); extract coefficients, t-stats,
#    p-values, adjusted p-values, B-statistics.
# 8. Visualize differential expression: volcano plots, PCA plots, and heatmaps
#    of the top 50 proteins for each contrast.
# 9. Save full topTable CSVs and significant protein lists per contrast.
# 10. Plot pathway-specific heatmaps for NF-κB and TGF-β pathway proteins.
# 11. Integrate RNA & proteomics DE: intersect leading-edge genes, compute
#    Spearman correlations, plot cross-omic scatter plots.
# 12. Perform gene-set analysis via limma::camera on Hallmark, KEGG, and
#    Reactome gene sets; collate results.
# 13. Export session info (package versions) for reproducibility.
# ==============================================================================

# ------------------------------------------------------------------------------
# Libraries ----
# ------------------------------------------------------------------------------
library(DEP)
library(limma)
library(edgeR)
library(fgsea)
library(msigdbr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(tidyverse)
library(readxl)
library(data.table)
library(pheatmap)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(svglite)

# ------------------------------------------------------------------------------
# Plot saving function ----
# ------------------------------------------------------------------------------
save_plot <- function(plot_obj, filename,
                      width = 8, height = 6, dpi = 300) {
  filepath <- file.path("DE_results", paste0(filename, ".svg"))
  if (inherits(plot_obj, c("Heatmap", "HeatmapList"))) {
    grDevices::svg(filepath, width = width, height = height, bg = "transparent")
    ComplexHeatmap::draw(plot_obj)
    dev.off()
  } else {
    ggplot2::ggsave(
      filename = filepath,
      plot     = plot_obj,
      device   = "svg",
      bg       = "transparent",
      dpi      = dpi,
      width    = width,
      height   = height,
      units    = "in"
    )
  }
}

# ------------------------------------------------------------------------------
# Data loading ----
# ------------------------------------------------------------------------------
anno <- read_xlsx(
  "~/Documents/Summer_proj_2024/Proteomics/BioMS_H2503_9121_dDIA_norm_250319_PGReport_R_OUT.xlsx", 
  sheet = 1
)
d_log <- read_xlsx(
  "~/Documents/Summer_proj_2024/Proteomics/BioMS_H2503_9121_dDIA_norm_250319_PGReport_R_OUT.xlsx", 
  sheet = 2
)
d_lognorm <- read_xlsx(
  "~/Documents/Summer_proj_2024/Proteomics/BioMS_H2503_9121_dDIA_norm_250319_PGReport_R_OUT.xlsx",
  sheet = 3
)

# ------------------------------------------------------------------------------
# Data wrangling ----
# ------------------------------------------------------------------------------
d_lognorm <- d_lognorm %>% filter(PG.Genes != "N/A")

duplicated_genes <- c("Macf1","Mocs2","Nrxn2","Tmpo","Tor1aip2")
d_long <- d_lognorm %>%
  pivot_longer(cols = 7:85, names_to = "sample", values_to = "intensity") %>%
  select(PG.Genes, PG.ProteinGroups, sample, intensity) %>%
  mutate(gene_label = if_else(PG.Genes %in% duplicated_genes,
                              paste0(PG.Genes, "_", PG.ProteinGroups),
                              PG.Genes))

gene_id_map <- d_long %>%
  distinct(gene_label, PG.Genes, PG.ProteinGroups) %>%
  rename(name = gene_label, symbol = PG.Genes, ID = PG.ProteinGroups)

dt_long <- as.data.table(d_long)
dt_wide <- dcast(dt_long, sample ~ gene_label, value.var = "intensity", fill = NA)

anno_subset <- anno %>%
  select(AnimalID, Condition, Sex, Genotype, Treatment) %>%
  mutate(
    Disease   = substr(Genotype, 3, nchar(Genotype)),
    Genotype  = substr(Genotype, 1, 2)
  )

dt_wide_annotated <- dt_wide %>%
  left_join(anno_subset, by = c("sample" = "AnimalID")) %>%
  relocate(Genotype, Disease, Sex, Treatment, .after = sample)

metadata_cols <- c("Genotype","Disease","Sex","Treatment")
expr_cols <- setdiff(colnames(dt_wide_annotated), c("sample", metadata_cols))

expr_matrix <- dt_wide_annotated %>%
  select(sample, all_of(expr_cols)) %>%
  column_to_rownames("sample") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("name") %>%
  left_join(gene_id_map, by = "name")

experimental_design <- dt_wide_annotated %>%
  as_tibble() %>%
  select(sample, Disease, Sex, Genotype, Treatment) %>%
  rename(label = sample) %>%
  mutate(condition = interaction(Disease, Sex, Genotype, Treatment, sep = "_")) %>%
  group_by(condition) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

protein_table <- expr_matrix %>%
  select(name, ID, all_of(experimental_design$label))
protein_table[,3:ncol(protein_table)] <-
  lapply(protein_table[,3:ncol(protein_table)], as.numeric)

LFQ_columns <- 3:ncol(protein_table)
se <- make_se(protein_table, LFQ_columns, experimental_design)

# ------------------------------------------------------------------------------
# Normalization ----
# ------------------------------------------------------------------------------
norm_plot <- plot_normalization(se)
print(norm_plot)
save_plot(norm_plot, "normalization")

# ------------------------------------------------------------------------------
# Missingness and expression ----
# ------------------------------------------------------------------------------
p_freq   <- plot_frequency(se)
p_num    <- plot_numbers(se)
p_miss   <- plot_missval(se)
p_detect <- plot_detect(se)

save_plot(p_freq,   "frequency")
save_plot(p_num,    "numbers")
save_plot(p_miss,   "missval")
save_plot(p_detect, "detect")

# ------------------------------------------------------------------------------
# Factor levels ----
# ------------------------------------------------------------------------------
colData(se)$Disease   <- relevel(factor(colData(se)$Disease), ref = "WT")
colData(se)$Treatment <- relevel(factor(colData(se)$Treatment), ref = "vehicle")
colData(se)$Genotype  <- relevel(factor(colData(se)$Genotype), ref = "E3")
colData(se)$Sex       <- relevel(factor(colData(se)$Sex), ref = "female")

# ------------------------------------------------------------------------------
# Imputation ----
# ------------------------------------------------------------------------------
se_minprob_05 <- impute(se, fun = "MinProb", q = 0.5)
save_plot(plot_imputation(se, se_minprob_05), "imputation_plot")

# ------------------------------------------------------------------------------
# Sample correlation ----
# ------------------------------------------------------------------------------
exprs  <- assay(se_minprob_05)
cor_mat <- cor(exprs, use = "pairwise.complete.obs", method = "pearson")

svg("DE_results/sample_correlation_heatmap.svg", width = 8, height = 8)
pheatmap(cor_mat,
         main = "Sample-to-sample Pearson correlation",
         annotation_col = as.data.frame(colData(se_minprob_05)[, c("Disease","Sex","Genotype","Treatment")]))
dev.off()

vars   <- apply(exprs, 1, var)
top500 <- names(sort(vars, decreasing = TRUE))[1:500]
cor500 <- cor(exprs[top500, ], use = "pairwise.complete.obs")

svg("DE_results/sample_correlation_heatmap_top500.svg", width = 8, height = 8)
pheatmap(cor500,
         main = "Sample-to-sample Pearson correlation (Top 500 variance proteins)",
         annotation_col = as.data.frame(colData(se_minprob_05)[, c("Disease","Sex","Genotype","Treatment")]))
dev.off()

# ------------------------------------------------------------------------------
# Full heatmap ----
# ------------------------------------------------------------------------------
pheatmap(
  exprs,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  scale           = "row",
  show_rownames   = FALSE,
  show_colnames   = TRUE,
  annotation_col  = as.data.frame(colData(se_minprob_05)[, c("Disease","Sex","Genotype","Treatment")]),
  fontsize_col    = 8,
)
# ------------------------------------------------------------------------------
# Differential Expression Analysis ----
# ------------------------------------------------------------------------------
se_fad_female <- se_minprob_05[, 
                               colData(se_minprob_05)$Disease == "FAD" & 
                                 colData(se_minprob_05)$Sex == "female"
]

# Function to plot contrast results
plotContrastResults <- function(contrast_name, fit_obj, se_obj, group_var, desc_title = NULL,
                                output_dir = "DE_results/", heatmap_width = 16, heatmap_height = 12) {
  
  # Section 2: Differential testing with limma
  topRes <- topTable(fit_obj, coef = contrast_name, number = Inf, sort.by = "P")
  
  # Section 3: Volcano plot
  volcano_plot <- topRes %>%
    rownames_to_column("Gene") %>%
    mutate(Significant = case_when(
      logFC > 0.5 & P.Value < 0.05  ~ "Up",
      logFC < -0.5 & P.Value < 0.05 ~ "Down",
      TRUE                           ~ "Not"
    )) %>%
    ggplot(aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.3) +
    labs(title = desc_title %||% paste("Volcano:", contrast_name),
         x = "log2FC", y = "-log10(P.Value)") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_label_repel(
      data = \(d) d |> filter(Significant != "Not"),
      aes(label = Gene),
      size = 5, colour = "black", fill = "white",
      label.size = 0.25, label.r = unit(0.15, "lines"),
      box.padding = 0.3, segment.color = "black", max.overlaps = Inf
    ) +
    coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0, 5))
  ggsave(
    filename = file.path(output_dir, paste0("volcano_", contrast_name, ".svg")),
    plot = volcano_plot,
    device = svglite::svglite,
    bg = "transparent",
    width = heatmap_width,
    height = heatmap_height / 2,
    units = "in",
    dpi = 300
  )
  
  # Section 4: PCA plot
  expr_data <- assay(se_obj)
  pca <- prcomp(t(expr_data), scale. = TRUE)
  pca_df <- as.data.frame(pca$x) %>%
    rownames_to_column("sample") %>%
    left_join(as.data.frame(colData(se_obj)) %>% rownames_to_column("sample"), by = "sample") %>%
    mutate(Genotype = factor(Genotype, levels = c("E3", "E4")))
  pca_plot <- ggplot(
    pca_df,
    aes(x = PC1, y = PC2, colour = .data[[group_var]], shape = Genotype)
  ) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(E3 = 16, E4 = 17)) +
    labs(title = "PCA") +
    theme_minimal()
  ggsave(
    filename = file.path(output_dir, paste0("pca_", contrast_name, ".svg")),
    plot = pca_plot,
    device = "svg",
    bg = "transparent",
    width = heatmap_width,
    height = heatmap_height / 2,
    units = "in",
    dpi = 300
  )
  
  # Section 5: Heatmap of top 50 proteins
  topRes50 <- head(topRes, 50)
  expr_top <- expr_data[rownames(topRes50), , drop = FALSE]
  annot_df <- as.data.frame(colData(se_obj))
  annot_df <- annot_df[, !(names(annot_df) %in% c("replicate", "condition", "ID", "label")), drop = FALSE]
  pheatmap(
    mat             = t(expr_top),
    scale           = "column",
    annotation_row  = annot_df,
    main            = desc_title %||% contrast_name,
    cellwidth       = 10,
    cellheight      = 8,
    filename        = file.path(output_dir, paste0("heatmap_", contrast_name, ".pdf")),
    width           = heatmap_width,
    height          = heatmap_height,
    silent          = TRUE
  )
  
  invisible(list(volcano = volcano_plot, pca = pca_plot))
}

# Fit limma model and perform contrasts
design <- model.matrix(~ Genotype * Treatment, data = colData(se_fad_female))
colnames(design) <- make.names(colnames(design))
fit <- lmFit(assay(se_fad_female), design) |> eBayes(trend = TRUE, robust = TRUE)
e4_contrast <- makeContrasts(LPS_E4 = TreatmentLPS + GenotypeE4.TreatmentLPS, levels = design)
fit_e4 <- contrasts.fit(fit, e4_contrast) |> eBayes(trend = TRUE, robust = TRUE)

# Plot results for model contrasts
plots_e3 <- plotContrastResults("TreatmentLPS", fit, se_fad_female, "Genotype", "LPS effect in E3FAD")
plots_e4 <- plotContrastResults("LPS_E4", fit_e4, se_fad_female, "Genotype", "LPS effect in E4FAD")
plots_int <- plotContrastResults("GenotypeE4.TreatmentLPS", fit, se_fad_female, "Treatment", "Interaction: LPS response E4 vs E3")

print(plots_e3$pca); print(plots_e4$pca); print(plots_int$pca)

#### Save plots ####
save_plot(plots_e3$volcano, "volcano_LPS_E3")
save_plot(plots_e4$volcano, "volcano_LPS_E4")
save_plot(plots_int$volcano, "volcano_interaction")
save_plot(plots_int$pca, "pca_interaction")


# ------------------------------------------------------------------------------
# Save topTables and significant proteins for each contrast ----
# ------------------------------------------------------------------------------
# Initialize lists
top_tables_list <- list()
sig_protein_list <- list()

# Define contrast names and corresponding coefficient IDs
contrast_info <- list(
  LPS_E3      = "TreatmentLPS",
  LPS_E4      = "LPS_E4",
  Interaction = "GenotypeE4.TreatmentLPS"
)

# Loop through contrasts
for (nm in names(contrast_info)) {
  coef_id <- contrast_info[[nm]]
  # Select fit object
  fit_obj <- switch(nm,
                    LPS_E3      = fit,
                    LPS_E4      = fit_e4,
                    Interaction = fit
  )
  
  # Extract full topTable
  tt <- topTable(fit_obj, coef = coef_id, number = Inf, sort.by = "P")
  top_tables_list[[nm]] <- tt
  
  # Filter significant proteins (P.Value < 0.05)
  sig_prot <- tt %>% filter(adj.P.Val < 0.05)
  sig_protein_list[[nm]] <- sig_prot
  
  # Save to CSV
  write.csv(tt,       file = paste0("DE_results/topTable_", nm, ".csv"),   row.names = TRUE)
  write.csv(sig_prot, file = paste0("DE_results/sigProteins_", nm, ".csv"), row.names = TRUE)
}

# ------------------------------------------------------------------------------
# Nf-kB/TNF Pathway Proteins Categorized by Function ----
# ------------------------------------------------------------------------------
receptors_adaptors <- c(
  "Tnfrsf1a",  # TNF receptor 1 (canonical NF-κB)
  "Tnfrsf1b",  # TNF receptor 2
  "Cd40",      # CD40 receptor (immune activation)
  "Tnfrsf13c", # BAFFR (non-canonical pathway)
  "Ltbr",      # Lymphotoxin-beta receptor (non-canonical NF-κB)
  "Il1r1",     # IL-1 receptor
  "Il1rap",    # IL-1 receptor accessory protein
  "Tradd",     # TNF receptor adaptor protein
  "Fadd",      # Adaptor in apoptotic and NF-κB signaling
  "Ripk1",     # Kinase bridging TNFR signaling to NF-κB or cell death
  "Ripk3",     # Necroptosis and inflammation
  "Traf2",     # Adaptor relaying TNFR signals to IKK
  "Tnfrsf8",   # CD30, involved in TNF-related signaling
  "Tnfrsf21"   # Death receptor 6 (DR6), possible NF-κB activation
)

ikk_complex <- c(
  "Chuk",     # IKKα (activates NF-κB via IκB degradation)
  "Ikbkb",    # IKKβ (key kinase in canonical pathway)
  "Ikbkg",    # NEMO (regulatory scaffold of IKK complex)
  "Map3k7",   # TAK1 (activates IKK complex downstream of TNFR/IL1R)
  "Tab1",     # TAK1-binding protein
  "Tab2",
  "Tab3",
  "Ikbke",    # IKKε (non-canonical pathway, innate immunity)
  "Map3k14"   # NIK (central in non-canonical NF-κB activation)
)

nfkb_subunits <- c(
  "Rela",     # NF-κB p65 (transcriptional activator)
  "Relb",     # Non-canonical NF-κB subunit
  "Rel",      # c-Rel (immune-specific NF-κB subunit)
  "Nfkb1",    # p105/p50 (canonical NF-κB DNA-binding)
  "Nfkb2"     # p100/p52 (non-canonical NF-κB)
)

nfkb_modulators <- c(
  "Nfkbia",   # IκBα (inhibitor of NF-κB)
  "Nfkbib",   # IκBβ (inhibitor of NF-κB)
  "Nfkbie",   # IκBε
  "Nfkbiz",   # IκBζ (NF-κB transcriptional co-factor)
  "Bcl3",     # Transcriptional co-activator (p50/p52)
  "Cebpb"     # C/EBPβ (cooperates with NF-κB in inflammation)
)

negative_regulators <- c(
  "Tnfaip3",  # A20 (negative regulator of NF-κB)
  "Cyld",     # Deubiquitinase, suppresses NF-κB
  "Otulin",   # Deubiquitinase, inhibits LUBAC-mediated activation
  "Sharpin",  # Component of LUBAC complex
  "Rbck1",    # HOIP, E3 ligase in LUBAC
  "Raif1",    # HOIL-1L, part of LUBAC
  "Fbwx11",   # Targets IκBα for degradation
  "Pias1",    # SUMO ligase, modulates NF-κB activity
  "Sumo3",    # SUMO protein modifying NF-κB components
  "Pdlim2"    # Terminates NF-κB signaling via nuclear export
)

misc_proteins <- c(
  "Commd1",   # Suppresses NF-κB transcriptional activity
  "Cul2",     # Ubiquitin ligase scaffold, interacts with NF-κB repressors
  "Il1b",     # IL-1β (NF-κB target and feedback activator)
  "Tlr2",     # Toll-like receptor 2 (activates NF-κB)
  "Tnf"       # TNFα (canonical activator of TNFR/NF-κB axis)
)

# Combine all groups
nfkb_genes <- c(
  receptors_adaptors,
  ikk_complex,
  nfkb_subunits,
  nfkb_modulators,
  negative_regulators,
  misc_proteins
)

# ------------------------------------------------------------------------------
# Check Nf-kB-related proteins in each contrast and print
# ------------------------------------------------------------------------------
for (nm in names(top_tables_list)) {
  tt <- top_tables_list[[nm]]
  subset <- tt[rownames(tt) %in% nfkb_genes, ] %>%
    arrange(abs(P.Value))
  print(paste("TopTable for contrast:", nm))
  print(subset)
}
# ------------------------------------------------------------------------------
# Heatmap for Nf-kB pathway proteins (scaled by row)
# ------------------------------------------------------------------------------
# Prepare expression matrix for selected NF-kB pathway genes
expr_subset <- assay(se_fad_female)[
  intersect(nfkb_genes, rownames(assay(se_fad_female))),
  , drop = FALSE
]

# Create and sort sample annotation
annot_df <- as.data.frame(colData(se_fad_female))[, c("Disease", "Sex", "Genotype", "Treatment")]
annot_df$sample <- rownames(annot_df)
annot_df <- annot_df %>%
  arrange(Genotype, Treatment)

# Reorder expression matrix columns to match sorted annotation
expr_ordered <- expr_subset[, annot_df$sample]

# Set sample names as rownames in annotation_col

# Define gene categories
pathway_lists <- list(
  ReceptorsAdaptors   = receptors_adaptors,
  IKKComplex           = ikk_complex,
  NfkbSubunits         = nfkb_subunits,
  NfkbModulators       = nfkb_modulators,
  NegativeRegulators   = negative_regulators,
  Miscellaneous        = misc_proteins
)

# Build row annotation (protein category labels)
row_annot_df <- bind_rows(
  lapply(names(pathway_lists), function(cat) {
    data.frame(Gene = pathway_lists[[cat]], Category = cat, stringsAsFactors = FALSE)
  })
)

# Keep only genes present in the ordered expression matrix
row_annot_df <- row_annot_df[row_annot_df$Gene %in% rownames(expr_ordered), ]
row_annot_df <- as.data.frame(row_annot_df)       
rownames(row_annot_df) <- NULL                
row_annot_df <- column_to_rownames(row_annot_df, var = "Gene")

# Save heatmap as SVG
svg("DE_results/heatmap_NFkB_TNF.svg", width = 8, height = 10)
pheatmap(
  expr_ordered,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annot_df[, !(colnames(annot_df) %in% "sample")],
  annotation_row = row_annot_df,
  main = "NF-kB / TNF Pathway Protein Expression",
  show_rownames = TRUE
)
dev.off()
# ------------------------------------------------------------------------------
# TGF-β Pathway Proteins Categorized by Function (Mouse, camelCase)
# ------------------------------------------------------------------------------

ligands <- c(
  "Tgfb1",  
  "Tgfb2",  
  "Tgfb3"   
)

receptors_coreceptors <- c(
  "Tgfbr1", 
  "Tgfbr2", 
  "Tgfbr3", 
  "Eng"    
)

core_smads <- c(
  "Smad2",  
  "Smad3",
  "Smad4"  
)

inhibitory_smads <- c(
  "Smad7",  
  "Smad6"   
)

extracellular_regulators <- c(
  "Ltbp1", "Ltbp2", "Ltbp3", "Ltbp4",  
  "Fbn1",                               
  "Thbs1",                             
  "Itgav", "Itgb6", "Itgb8",          
  "Furin"                              
)

smad_modulators <- c(
  "Zfyve9",  
  "Trim33",  
  "Ski", "Skil", 
  "Smurf1", "Smurf2",  
  "Usp15"    
)

non_smad_branch <- c(
  "Traf6",           
  "Map3k7",          
  "Tab1", "Tab2", "Tab3",
  "Mapk14",         
  "Mapk8", "Mapk9",  
  "Rhoa", "Rock1", "Rock2" 
)

transcriptional_cofactors <- c(
  "Crebbp", "Ep300", 
  "Sp1",    
  "Jun",            
  "Foxp3",         
  "E2f4"          
)

negative_regulators <- c(
  "Nedd4l", 
  "Ppp1ca", "Ppp2ca", 
  "Dcaf7"           
)

# Combine all groups into one master vector
# ------------------------------------------------------------------------------
tgfb_genes <- c(
  ligands,
  receptors_coreceptors,
  core_smads,
  inhibitory_smads,
  extracellular_regulators,
  smad_modulators,
  non_smad_branch,
  transcriptional_cofactors,
  negative_regulators
)

# Check TGF-β–related proteins in each contrast and print
# ------------------------------------------------------------------------------
for (nm in names(top_tables_list)) {
  tt <- top_tables_list[[nm]]
  subset <- tt[rownames(tt) %in% tgfb_genes, ] %>%
    arrange(abs(P.Value))
  print(paste("TopTable for contrast:", nm))
  print(subset)
}

# Heatmap for TGF-β pathway proteins (scaled by row)
# ------------------------------------------------------------------------------
# Prepare expression matrix for selected TGF-β pathway genes
expr_subset <- assay(se_fad_female)[
  intersect(tgfb_genes, rownames(assay(se_fad_female))),
  , drop = FALSE
]

# Create and sort sample annotation
annot_df <- as.data.frame(colData(se_fad_female))[, c("Disease", "Sex", "Genotype", "Treatment")]
annot_df$sample <- rownames(annot_df)
annot_df <- annot_df %>%
  arrange(Genotype, Treatment)

# Reorder expression matrix columns to match sorted annotation
expr_ordered <- expr_subset[, annot_df$sample]

# Define TGF-β pathway gene categories
pathway_lists <- list(
  Ligands                  = ligands,
  ReceptorsCoreceptors     = receptors_coreceptors,
  CoreSmads                = core_smads,
  InhibitorySmads          = inhibitory_smads,
  ExtracellularRegulators  = extracellular_regulators,
  SmadModulators           = smad_modulators,
  NonSmadBranch            = non_smad_branch,
  TranscriptionalCofactors = transcriptional_cofactors,
  NegativeRegulators       = negative_regulators
)

# Build row annotation (gene category labels)
row_annot_df <- bind_rows(
  lapply(names(pathway_lists), function(cat) {
    data.frame(Gene = pathway_lists[[cat]], Category = cat, stringsAsFactors = FALSE)
  })
)

# Filter to genes actually present
row_annot_df <- row_annot_df[row_annot_df$Gene %in% rownames(expr_ordered), ]
row_annot_df <- as.data.frame(row_annot_df)
rownames(row_annot_df) <- NULL
row_annot_df <- column_to_rownames(row_annot_df, var = "Gene")

# Generate heatmap
svg("DE_results/heatmap_TGFb.svg", width = 8, height = 10)
pheatmap(
  expr_ordered,
  scale            = "row",
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  annotation_col   = annot_df[, setdiff(colnames(annot_df), "sample")],
  annotation_row   = row_annot_df,
  main             = "TGF-β Pathway Protein Expression",
  show_rownames    = TRUE
)
dev.off()

# ------------------------------------------------------------------------------
# Plotting of individual proteins ----
# ------------------------------------------------------------------------------
# Pull expression values for Tradd
tradd_expr <- assay(se_minprob_05)["Tradd", ]
meta <- as.data.frame(colData(se_minprob_05))
meta$expr <- tradd_expr
meta$group <- interaction(meta$Genotype, meta$Treatment, sep = "_")

# Plot
ggplot(meta, aes(x = group, y = expr, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  labs(title = "Tradd expression by Genotype and Treatment",
       x = "Genotype_Treatment",
       y = "log2 LFQ intensity") +
  theme_minimal()

# ------------------------------------------------------------------------------
# Plotting individual proteins with DEP plot_single-like plot ----
# ------------------------------------------------------------------------------
proteins <- c(    "Ccl3", "Ccl4",
                  "Nfkbia", "Nfkbiz", "Nfkb1",
                  "Il1a", "Il1b",
                  "Tnf", "Tnfaip3",
                  "Atf3", "Gadd45b",
                  "Ccrl2", "Bcl2a1b",
                  "Tlr2", "Zfp36",
                  "Cd83", "Csf1")

chunks <- split(proteins, ceiling(seq_along(proteins) / 6))


plot_group <- function(prot_group, se, new_order) {
  prot_group <- intersect(prot_group, rownames(se))
  if (length(prot_group) == 0) {
    warning("No proteins from this group found in the dataset – skipping.")
    return(NULL)
  }
  mat_raw <- assay(se)
  
  matc <- mat_raw - rowMeans(mat_raw, na.rm=TRUE)
  df_pts <- matc %>%
    as.data.frame() %>%
    rownames_to_column("protein") %>%
    gather(sample, value, -protein) %>%
    filter(protein %in% prot_group) %>%
    left_join(
      colData(se) %>% 
        as.data.frame() %>% 
        rownames_to_column("sample"),
      by = "sample"
    ) %>%
    mutate(
      replicate = factor(replicate),
      condition = factor(condition, levels = new_order)
    )
  df_sum <- df_pts %>%
    group_by(protein, condition) %>%
    summarize(
      mean   = mean(value, na.rm=TRUE),
      sd     = sd(value,   na.rm=TRUE),
      n      = n(),
      error  = qnorm(0.975) * sd / sqrt(n),
      CI.L   = mean - error,
      CI.R   = mean + error
    ) %>%
    ungroup()
  
  ggplot(df_sum, aes(condition, mean)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_col(colour = "black", fill = "grey80") +
    geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.2) +
    geom_point(
      data = df_pts,
      aes(condition, value, colour = replicate),
      position = position_jitter(width = 0.1, height = 0),
      size = 2
    ) +
    facet_wrap(~ protein, scales = "free_x", ncol = 2) +
    labs(
      x      = "Condition",
      y      = expression(log[2]~"Intensity"~"(±95% CI)"),
      colour = "Replicate"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      panel.grid       = element_line(colour = "grey90")
      # axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

new_order <- c("FAD_female_E3_vehicle", "FAD_female_E3_LPS",
               "FAD_female_E4_vehicle",   "FAD_female_E4_LPS")

plots <- lapply(chunks, function(prot_group) {
  plot_group(prot_group, se = se_fad_female, new_order = new_order)
})

for(p in plots) print(p)

# ------------------------------------------------------------------------------
# Transcriptomic comparison ----
# ------------------------------------------------------------------------------
export_dir <- "/Users/loganlarlham/Documents/Summer_proj_2024/Notebooks/Analysis/export"

dc_files <- list.files(
  path       = export_dir,
  pattern    = "^dc_vehlps_.*\\.csv$",
  full.names = TRUE
)
dc_results <- lapply(dc_files, read.csv, row.names = 1)
names(dc_results) <- basename(dc_files) %>%
  sub("^dc_vehlps_(.*)\\.csv$", "\\1", .)
print(names(dc_results))  # e.g. "E3FAD_vehlps", etc.

gsea_files <- list.files(
  path       = export_dir,
  pattern    = "^gsea_.*\\.csv$",
  full.names = TRUE
)
gsea_results <- lapply(gsea_files, read.csv, row.names = 1)
names(gsea_results) <- basename(gsea_files) %>%
  sub("^gsea_(.*)\\.csv$", "\\1", .)
print(names(gsea_results))

gprot_files <- list.files(
  path       = export_dir,
  pattern    = "^gsea_prot_.*\\.csv$",
  full.names = TRUE
)
gsea_prot_results <- lapply(gprot_files, read.csv, row.names = 1)
names(gsea_prot_results) <- basename(gprot_files) %>%
  sub("^gsea_prot_(.*)\\.csv$", "\\1", .)
print(names(gsea_prot_results))


rna_contrasts  <- c(E3 = "E3FAD_vehlps", E4 = "E4FAD_vehlps")
prot_contrasts <- c(E3 = "LPS_E3",       E4 = "LPS_E4")


for (geno in names(rna_contrasts)) {
  rna_nm  <- rna_contrasts[geno]
  prot_nm <- prot_contrasts[geno]
  
  grna  <- gsea_results[[rna_nm]]
  gprot <- gsea_prot_results[[prot_nm]]
  
  le_rna_all   <- unique(unlist(grna$leadingEdge))
  common_genes <- intersect(le_rna_all,
                            rownames(top_tables_list[[prot_nm]]))
  
  rna_df <- dc_results[[rna_nm]] %>%
    as_tibble(rownames = "gene") %>%
    filter(gene %in% common_genes) %>%
    select(gene, log2FoldChange) %>%
    rename(rna_log2FC = log2FoldChange)
  
  prot_df <- top_tables_list[[prot_nm]] %>%
    as_tibble(rownames = "gene") %>%
    filter(gene %in% common_genes) %>%
    select(gene, logFC) %>%
    rename(prot_log2FC = logFC)
  
  crossomic_df <- inner_join(rna_df, prot_df, by = "gene")
  
  nfkb_genes <- unlist(
    grna$leadingEdge[ grna$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ]
  )
  tgfb_genes <- unlist(
    grna$leadingEdge[ grna$pathway == "HALLMARK_TGF_BETA_SIGNALING" ]
  )
  
  crossomic_df <- crossomic_df %>%
    mutate(
      pathway = case_when(
        gene %in% nfkb_genes & gene %in% tgfb_genes ~ "Both",
        gene %in% nfkb_genes                  ~ "TNFA",
        gene %in% tgfb_genes                  ~ "TGF"
      ),
      pathway = factor(pathway, levels = c("TNFA","TGF","Both"))
    )
  
  label_df <- crossomic_df %>%
    filter(pathway %in% c("TNFA", "TGF")) %>%
    mutate(pathway = droplevels(pathway)) %>%
    group_by(pathway) %>%
    summarize(
      n_obs = sum(!is.na(rna_log2FC) & !is.na(prot_log2FC)),
      rho   = if (n_obs >= 3) 
        cor(rna_log2FC, prot_log2FC, method = "spearman", use = "complete.obs")
      else 
        NA_real_,
      pval  = if (n_obs >= 3) 
        cor.test(rna_log2FC, prot_log2FC, method = "spearman", exact = FALSE)$p.value
      else 
        NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      rho_str = format(rho,  digits = 2),
      pval_str = format(pval, digits = 2),
      x_pos = -0.5,
      y_pos = if_else(pathway == "TNFA", 1.00, 0.75),
      label = paste0("rho == ", rho_str, "~','~p == ", pval_str)
    )
  
  p <- ggplot(crossomic_df, aes(
    x     = rna_log2FC,
    y     = prot_log2FC,
    color = pathway
  )) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c(TNFA = "#D55E00", TGF = "#0072B2", Both = "#009E73")
    ) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    geom_text_repel(aes(label = gene), size = 3) +
    geom_text(
      data    = label_df,
      aes(x = x_pos, y = y_pos, label = label, color = pathway),
      parse   = TRUE,
      hjust   = 1,
      size    = 4
    ) +
    coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
    labs(
      title = paste("Spearman correlation of RNA vs Protein log2FC", geno, "Leading edge genes"),
      x     = "RNA log2 fold-change",
      y     = "Protein log2 fold-change",
      color = "Pathway\nmembership"
    ) +
    theme_minimal()
  
  save_plot(p, paste0("spearman_leading_edge_", geno))
}

# ------------------------------------------------------------------------------
# Transcriptomic comparison for full pathway genes ------
# ------------------------------------------------------------------------------
hs_hallmark <- msigdbr(
  species  = "Mus musculus",
  collection = "H"
) %>% 
  select(gs_name, gene_symbol)

full_nfkb <- hs_hallmark %>% 
  filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>% 
  pull(gene_symbol) %>% unique()

full_tgfb <- hs_hallmark %>% 
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>% 
  pull(gene_symbol) %>% unique()

for(geno in names(rna_contrasts)) {
  rna_nm  <- rna_contrasts[geno]
  prot_nm <- prot_contrasts[geno]

  grna  <- gsea_results[[rna_nm]]
  gprot <- gsea_results[[prot_nm]]
  
  genes_union <- union(full_nfkb, full_tgfb)

  rna_df <- as_tibble(dc_results[[rna_nm]], rownames = "gene") %>%
    filter(gene %in% genes_union) %>%
    select(gene, log2FoldChange) %>%
    rename(rna_log2FC = log2FoldChange)
  
  prot_df <- as_tibble(top_tables_list[[prot_nm]], rownames = "gene") %>%
    filter(gene %in% genes_union) %>%
    select(gene, logFC) %>%
    rename(prot_log2FC = logFC)

  crossomic_df <- inner_join(rna_df, prot_df, by = "gene")
  
  crossomic_df <- crossomic_df %>%
    mutate(pathway = case_when(
      gene %in% full_nfkb & gene %in% full_tgfb ~ "Both",
      gene %in% full_nfkb                  ~ "TNFA",
      gene %in% full_tgfb                  ~ "TGF"
    )) %>%
    mutate(pathway = factor(pathway, levels = c("TNFA", "TGF", "Both")))
  
  label_df <- crossomic_df %>%
    filter(pathway %in% c("TNFA", "TGF")) %>%
    droplevels() %>%
    split(.$pathway) %>%
    map_df(~ {
      ct <- cor.test(
        x      = .x$rna_log2FC,
        y      = .x$prot_log2FC,
        method = "spearman"
      )
      tibble(
        pathway = .x$pathway[1],
        rho     = format(ct$estimate, digits = 2),
        pvalue  = format(ct$p.value, digits   = 2),
        label  = paste0("rho == ", rho, "~','~p == ", pvalue),
        x_pos   = -0.5,
        y_pos  = if_else(pathway == "TNFA",  1.00,  
                         0.75)  
      )
    })
  print("full pathway corr")
  print(label_df)

  p <- ggplot(crossomic_df, aes(
    x     = rna_log2FC,
    y     = prot_log2FC,
    color = pathway
  )) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c(TNFA = "#D55E00", TGF = "#0072B2", Both = "#009E73")
    ) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    geom_text_repel(aes(label = gene), size = 3) +
    geom_text(
      data    = label_df,
      aes(x = x_pos, y = y_pos, label = label, color = pathway),
      parse   = TRUE,
      hjust   = 1,
      size    = 4
    ) +
    coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
    labs(
      title = paste("Spearman correlation of RNA vs Protein log2FC in", geno, "all genes"),
      x     = "RNA log2 fold-change",
      y     = "Protein log2 fold-change",
      color = "Pathway\nmembership"
    ) +
    theme_minimal()
  
  print(p)
  ggplot2::ggsave(
    filename = paste0("DE_results/spearman_allpathway_", geno, ".pdf"),
    plot     = p,
    device   = "pdf",           
    bg       = "transparent",  
    width    = 12,              
    height   = 10,            
    units    = "in",
    dpi      = 300
  )
}
# ------------------------------------------------------------------------------
# Transcriptomic comparison for All logfc changes ---------
# ------------------------------------------------------------------------------
for (geno in names(rna_contrasts)) {
  rna_nm  <- rna_contrasts[geno]
  prot_nm <- prot_contrasts[geno]
  

  rna_df <- as_tibble(dc_results[[rna_nm]], rownames = "gene") %>%
    select(gene, log2FoldChange) %>%
    rename(rna_log2FC = log2FoldChange)
  

  prot_df <- as_tibble(top_tables_list[[prot_nm]], rownames = "gene") %>%
    select(gene, logFC) %>%
    rename(prot_log2FC = logFC)
  

  crossomic_df <- inner_join(rna_df, prot_df, by = "gene")
  

  ct <- cor.test(crossomic_df$rna_log2FC,
                 crossomic_df$prot_log2FC,
                 method = "spearman")
  rho    <- format(ct$estimate, digits = 2)
  pval   <- format(ct$p.value, digits = 2)
  label  <- paste0("rho == ", rho, "~','~p == ", pval)
  

  p <- ggplot(crossomic_df, aes(x = rna_log2FC, y = prot_log2FC)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    annotate("text", x = -0.5, y = max(crossomic_df$prot_log2FC, na.rm = TRUE),
             label = label, parse = TRUE, hjust = 0, size = 4) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-2, 2)) +
    labs(
      title = paste("Spearman correlation of RNA vs Protein log2FC values after LPS in", geno ),
      x     = "RNA log2 fold-change",
      y     = "Protein log2 fold-change"
    ) +
    theme_minimal()
  
  ggplot2::ggsave(
    filename = paste0("DE_results/spearman_allDiffExp_", geno, ".pdf"),
    plot     = p,
    device   = "pdf",           
    bg       = "transparent",   
    width    = 12,             
    height   = 6,             
    units    = "in",
    dpi      = 300
  )
  print(p)
}
# ------------------------------------------------------------------------------
# Transcriptomic comparison for significant LogFC -----------
# ------------------------------------------------------------------------------
for (geno in names(rna_contrasts)) {
  print(geno)
  rna_nm  <- rna_contrasts[geno]
  prot_nm <- prot_contrasts[geno]
  
  # Read in significant gene list exported from Python
  sig_df <- read.csv(
    file.path(export_dir, paste0("dc_vehlps_sig_", geno, "FAD_vehlps.csv")),
    row.names = 1,
    stringsAsFactors = FALSE
  )
  sig_genes <- rownames(sig_df)
  
  # assemble RNA DE table, filter to significant genes
  rna_df <- as_tibble(dc_results[[rna_nm]], rownames = "gene") %>%
    filter(gene %in% sig_genes) %>%
    select(gene, log2FoldChange) %>%
    rename(rna_log2FC = log2FoldChange)
  
  # assemble Prot DE table, filter to significant genes
  prot_df <- as_tibble(top_tables_list[[prot_nm]], rownames = "gene") %>%
    filter(gene %in% sig_genes) %>%
    select(gene, logFC) %>%
    rename(prot_log2FC = logFC)
  
  # merge
  crossomic_df <- inner_join(rna_df, prot_df, by = "gene")
  
  # build Spearman label
  ct <- cor.test(crossomic_df$rna_log2FC,
                 crossomic_df$prot_log2FC,
                 method = "spearman")
  rho    <- format(ct$estimate, digits = 2)
  pval   <- format(ct$p.value, digits = 2)
  label  <- paste0("rho == ", rho, "~','~p == ", pval)
  
  # plot
  p <- ggplot(crossomic_df, aes(x = rna_log2FC, y = prot_log2FC)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    geom_text_repel(aes(label = gene), size = 3) +
    annotate("text", x = -2.5, y = 1.5,
             label = label, parse = TRUE, hjust = 0, size = 4) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
    labs(
      title = paste("Spearman correlation of significant RNA vs Protein log2FC in", geno, "Sig DiffExp"),
      x     = "RNA log2 fold-change",
      y     = "Protein log2 fold-change"
    ) +
    theme_minimal()
  
  print(p)
  ggplot2::ggsave(
    filename = paste0("DE_results/spearman_siglogfc_", geno, ".pdf"),
    plot     = p,
    device   = "pdf",           # base‐R pdf() device
    bg       = "transparent",   # make the background transparent
    width    = 12,              # inches (tweak to suit)
    height   = 10,              # inches
    units    = "in",
    dpi      = 300              # ignored for vector output, but harmless
  )
}
# ------------------------------------------------------------------------------
# Limma::romer instead of GSEA on Hallmark ----
# ------------------------------------------------------------------------------
exprs <- assay(se_fad_female)     


msig     <- msigdbr(species = "Mus musculus", collection = "H")
hallmark <- split(msig$gene_symbol, msig$gs_name)
index    <- ids2indices(hallmark, rownames(exprs), remove.empty = TRUE)

contrast_list <- list(
  LPS_E3      = 3,         
  LPS_E4      = e4_contrast, 
  Interaction = 4        
)


camera_results_H <- list()

for (nm in names(contrast_list)) {
  con <- contrast_list[[nm]]
  message("=== camera: ", nm, " ===")
  
  cam <- camera(
    y        = exprs,
    index    = index,
    design   = design,
    contrast = con,

    sort     = TRUE               
  )
  
  print(head(cam, 10))
  
  camera_results_H[[nm]] <- cam
}

# ---------------------------------------------------------------------
# Limma::camera on KEGG MEDICUS gene sets ----
# ---------------------------------------------------------------------


exprs <- assay(se_fad_female)   


msig       <- msigdbr(species = "Mus musculus",
                      collection = "C2",
                      subcollection = "CP:KEGG_MEDICUS")
kegg_med   <- split(msig$gene_symbol, msig$gs_name)


is_ref <- sapply(names(kegg_med), function(pathway_name) {
  parts <- strsplit(pathway_name, "_")[[1]]
  length(parts) >= 3 && parts[3] == "REFERENCE"
})


kegg_med_ref <- kegg_med[is_ref]


length(kegg_med)    
length(kegg_med_ref)  
names(kegg_med_ref)[1:10]

index      <- ids2indices(kegg_med_ref, rownames(exprs), remove.empty = TRUE)


contrast_list <- list(
  LPS_E3      = 3,
  LPS_E4      = e4_contrast,
  Interaction = 4
)

camera_results_kegg <- list()

for (nm in names(contrast_list)) {
  con <- contrast_list[[nm]]
  message("=== camera (KEGG): ", nm, " ===")
  
  cam <- camera(
    y        = exprs,
    index    = index,
    design   = design,
    contrast = con,
    sort     = TRUE      

  )
  
  print(head(cam, 10))    
  camera_results_kegg[[nm]] <- cam
}

camera_df <- bind_rows(
  camera_results_kegg[["LPS_E3"]] %>%
    rownames_to_column("Pathway") %>%
    mutate(Contrast = "LPS effect in E3"),
  camera_results_kegg[["LPS_E4"]] %>%
    rownames_to_column("Pathway") %>%
    mutate(Contrast = "LPS effect in E4"),
  camera_results_kegg[["Interaction"]] %>%
    rownames_to_column("Pathway") %>%
    mutate(Contrast = "Interaction: LPS response in E4 vs E3")
) %>%
  group_by(Pathway) %>%
  filter(any(FDR < 0.25)) %>%
  ungroup() %>%
  mutate(
    ShapeCat = case_when(
      FDR < 0.05  ~ "FDR < 0.05",
      FDR < 0.25  ~ "FDR < 0.25",
      TRUE        ~ "Not significant"
    ),
    ShapeCat = factor(
      ShapeCat,
      levels = c("FDR < 0.05", "FDR < 0.25", "Not significant")
    )
  )

p <- ggplot() +
  geom_point(
    data = camera_df,
    aes(
      x      = Contrast, 
      y      = Pathway,
      colour = Direction,
      shape  = ShapeCat
    ),
    size = 3
  ) +
  ggstar::geom_star(
    data      = filter(camera_df, ShapeCat == "FDR < 0.05"),
    aes(x = Contrast, y = Pathway, fill = Direction),
    starshape   = 1,
    size        = 5,
    colour      = "black",
    show.legend = FALSE 
  ) +
  scale_shape_manual(
    name   = "FDR threshold",
    drop   = FALSE,
    values = c(
      "FDR < 0.05"      = 8,  
      "FDR < 0.25"      = 17,  
      "Not significant" = 16   
    ),
    guide = guide_legend(
      override.aes = list(
        shape  = c(8,17,16),  
        colour = "black",     
        size   = 5    
      )
    )
  ) +
  scale_colour_manual(
    name   = "Direction",
    values = c(Up = "red", Down = "blue")
  ) +
  scale_fill_manual(
    name  = "Direction",
    values = c(Up = "red", Down = "blue"),
    guide = FALSE
  ) +
  labs(
    title = "Pathway Enrichment Analysis",
    x     = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y     = element_text(size = 10),
    axis.text.x     = element_text(angle = 25, hjust = 1, size = 14),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p)


ggsave(
  filename = "DE_results/kegg_camera_dotplot_with_stars_and_legend.svg",
  plot     = p,
  width    = 10,
  height   = 12,
  device   = "svg"
)

# ---------------------------------------------------------------------
# Limma::camera on Reactome gene sets ----
# ---------------------------------------------------------------------


exprs <- assay(se_fad_female)


msig      <- msigdbr(species = "Mus musculus",
                     collection = "C2",
                     subcollection = "CP:REACTOME")
reactome  <- split(msig$gene_symbol, msig$gs_name)
index     <- ids2indices(reactome, rownames(exprs), remove.empty = TRUE)


contrast_list <- list(
  LPS_E3      = 3,
  LPS_E4      = e4_contrast,
  Interaction = 4
)

camera_results_reac <- list()

for (nm in names(contrast_list)) {
  con <- contrast_list[[nm]]
  message("=== camera (Reactome): ", nm, " ===")
  
  cam <- camera(
    y        = exprs,
    index    = index,
    design   = design,
    contrast = con,
    sort     = TRUE
  )
  
  print(head(cam, 20))
  camera_results_reac[[nm]] <- cam
}
# ------------------------------------------------------------------------------
# Save session info ----
# ------------------------------------------------------------------------------
writeLines(capture.output(session_info()), 'DE_results/session_info.txt')
