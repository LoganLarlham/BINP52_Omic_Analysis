```text
   _    ____    __  __ _                      _ _        
    / \  |  _ \  |  \/  (_) ___ _ __ ___   __ _| (_) __ _  
   / _ \ | | | | | |\/| | |/ __| '__/ _ \ / _` | | |/ _` | 
  / ___ \| |_| | | |  | | | (__| | | (_) | (_| | | | (_| | 
 /_/   \_\____/  |_|  |_|_|\___|_|  \___/ \__, |_|_|\__,_| 
  __  __       _ _   _                    |___/            
 |  \/  |_   _| | |_(_)       ___  _ __ ___ (_) ___ ___    
 | |\/| | | | | | __| |_____ / _ \| '_ ` _ \| |/ __/ __|   
 | |  | | |_| | | |_| |_____| (_) | | | | | | | (__\__ \   
 |_|__|_|\__,_|_|\__|_| _    \___/|_| |_| |_|_|\___|___/   
 |_   _| |__   ___  ___(_)___                              
   | | | '_ \ / _ \/ __| / __|                             
   | | | | | |  __/\__ \ \__ \                             
   |_| |_| |_|\___||___/_|___/                             
                                                                                
```                     
https://github.com/LoganLarlham/BINP52_Omic_Analysis
                                                                    
#  Overview  
  transcriptomic and proteomic analysis of LPS response in EFAD mice.

  ##  Environments  
   The environment.yml contains all dependencies for Python-based transcriptomic analysis.  
   The renv lockfiles include all R packages required for proteomic analysis.  

  ##  Project Structure  
   data/  
     transcriptomic/  
       raw/         
       processed/   
     proteomic/  
       raw/          
       processed/   
   scripts/  
     transcriptomic/ # Jupyter notebook and helper Python modules  
     proteomic/      # R analysis script  
   results/  
     transcriptomic/  
       figures/  tables/  
     proteomic/  
       figures/  tables/  

   Usage  
   scripts/transcriptomic/ contains one large  notebook that performs all  quality control, clustering, pseudobulk differential expression, and gene set enrichment analysis of microglial single-cell RNA-seq data  
   scripts/proteomic/analysis.R runs the mass spectrometry–based proteomic analysis pipeline  

   Raw Data  
   Note  data/*/raw/ contains large files that are not tracked due to size.  Request access if needed.  

   Analysis Overview  
   Two complementary pipelines are implemented to characterize microglial responses to lipopolysaccharide priming in FAD mice with human APOE alleles  
   Transcriptomic Analysis (scripts/transcriptomic)  
     • single-cell RNA sequencing of isolated microglia  
     • quality filtering, normalization, and mitochondrial/ribosomal regression using Scanpy and Anndata  
     • dimensionality reduction (PCA, UMAP) and Leiden clustering for identification of distinct microglial states  
     • differential abundance testing by three-way ANOVA using statsmodels  
     • pseudobulk differential expression with Decoupler and PyDESeq2  
     • gene set enrichment analysis on Wald statistics with Decoupler and MSigDB hallmark sets  
   Proteomic Analysis (scripts/proteomic/analysis.R)  
     • data acquired by DIA–PASEF on timsTOF HT and processed in Spectronaut v19  
     • log2-transformed label-free quantification (LFQ) values imported into R and structured with SummarizedExperiment  
     • quality control, missing value imputation (MinProb) and normalization with DEP  
     • differential abundance analysis using limma with empirical Bayes shrinkage  
     • competitive gene set testing against KEGG pathways using Camera  

   Tools and Dependencies  
   Python environment (environment.yml)  
     • scanpy, anndata, numpy, pandas, scipy, statsmodels, decoupler, pydeseq2, matplotlib  
   R environment (renv)  
     • tidyverse (tidyr, dplyr), SummarizedExperiment, DEP, limma, camera, ggplot2  

   Code Structure  
   transcriptomic notebooks  
     scrnaseq_analysis.ipynb
     func_lib.py  

   proteomic R script  
     analysis.R  

   Results and Outputs  
   Processed matrices, figures, and tables are available under results/transcriptomic and results/proteomic directories  
   Generated figures include UMAP embeddings, cluster proportion plots, volcano plots, GSEA heatmaps, PCA of proteomic data, volcano and pathway enrichment plots
