#import "@local/wordometer:0.1.5": word-count, total-words
#show: word-count.with(exclude: heading)

#set text(font: "TeX Gyre Termes", size: 10pt)
#show heading.where(level: 1): set text(size: 18pt)
#show heading.where(level: 2): set text(size: 16pt)
#show heading.where(level: 3): set text(size: 14pt, style: "italic", fill: gray)
#show figure.caption: set text(weight: "bold")
#show figure.caption: set par(leading: 1em)
#show figure.caption: set align(left)
#show <n>: set text(fill: red)




#show figure.where(
  kind: table
): set block(breakable: false)
#show figure.where(
  kind: table
): set figure.caption(position: top)

#show figure.where(
  kind: image
): set figure.caption(position: bottom)

#show heading: set block(above: 1.5em, below: 1.5em)



#set page(numbering: none)

#align(center + horizon, text(20pt)[
  _*#lorem(20)*_
])
\
  #align(center)[#text(size: 18pt)[Logan Larlham]
  \ 
  \
  Supervised by:
  \
  #lorem(2)
  \
  #lorem(2)
  \
  \
  #text(16pt)[Master’s in bioinformatics, Course BINP52, 60cr]]


  #block(move(dy:175pt)[Email: #underline[#lorem(2)] \
  Supervisor email: #underline[#lorem(2)] \ Experimental Neuroinflammation Laboratory, Department of Experimental Medical Sciences, Faculty of Medicine, Lund University, Lund, Sweden\
  In this document there are #total-words words])
#pagebreak()

#set par(spacing: 3em)
#set par(leading: 2em,  justify: true, first-line-indent: 2em)
#set page(numbering: "1/1", number-align: right,)
#counter(page).update(1) 

#let note(body, size: 8pt, fill: rgb("EAF2F5")) = {
  rect(
    stroke: luma(80%) + 0.5pt,
    inset: 1mm,
    fill: fill,
    radius: 3pt,
    width: 100%,
    text(size: size, fill:rgb("FF0000"), body)
  )
}


#set par.line(numbering: "1")
= Abstract <abstract>

Alzheimer’s disease (AD) is characterized by amyloid-β (Aβ) plaque deposition, tau pathology, and neuroinflammation, with genetic variation in the APOE gene representing the strongest risk factor for late-onset AD. Microglia, the brain’s resident immune cells, are increasingly recognized as central to disease progression through their roles in neuroinflammation and Aβ clearance. In this study, we investigated how APOE genotype (E3 vs E4) and early-life inflammation interact to shape microglial states in a mouse model of Aβ pathology. Using single-cell RNA sequencing of microglia isolated from EFAD mice treated with lipopolysaccharide (LPS) or vehicle at postnatal day 9, we identified distinct microglial transcriptional states and assessed their proportions and gene expression profiles. E4FAD mice exhibited a higher proportion of disease-associated microglia (DAM) compared to E3FAD mice. Moreover, early-life inflammation induced long-lasting, genotype-specific alterations: LPS reduced the abundance of an NF-κB-enriched Cytokine Response microglial population and downregulated pro-inflammatory genes in E4FAD, but not E3FAD, mice. These findings suggest that early-life immune events interact with APOE genotype to shape microglial responses to amyloid pathology and may offer targets for modulating disease-associated inflammation.

= Introduction/Background <introduction>

=== Alzheimer's Disease <alzheimer-disease>
Alzheimer's disease (AD) is a neurodegenerative disorder that is characterized by progressive cognitive decline and memory loss, and the foremost cause of dementia worldwide affecting an estimated 55.2 million people in 2019 @noauthor_global_2021. While historically AD has been defined based on clinical symptoms as well as neuropathological hallmarks usually identified post-mortem, in recent years the field has shifted towards defining AD as a disease based on the presence of the hallmark proteinopathies of amyloid-β (Aβ) plaques and tau neurofibrillary tangles formalized in the National Institute on Aging and Alzheimer's Association (NIA-AA) research framework @jack_jr_nia-aa_2018. There are two main forms of AD known as Early Onset AD (EOAD) and Late Onset AD (LOAD), defined by an arbitrary diagnosis age cutoff of 60 years. It is estimated that 5-10% of AD cases are EOAD whereas LOAD is the much more typical form which accounts for 90-95% of all cases. While EOAD has an estimated heritability of 90-100%, and often caused by known autosomal dominant mutations in the APP, PSEN1, and PSEN2 genes, LOAD is a complex disorder with heritability estimated between 40-80% and is influenced by a combination of genetic, environmental, and lifestyle factors @chen_genetically_2021. Another important aspect of AD is the sex differences in the disease. There are signifcantly more women with AD however, however when accounting for overall higher mortality rates men and age-adjusted prevalence some studies have shown no difference, and some have shown a higher prevalence in women. Additionally there have been demonstrated differences in neuropathologies, Aβ and tau burden, and cognitive decline. The reasons for these differences are not well understood but are thought to be due to a combination of genetic, hormonal, and environmental factors @oveisgharan_sex_2018.

As mentioned the pathological hallmarks of AD are the extracellular deposition and accumulation of aggregated amyloid-β (Aβ) peptides into plaques and the intracellular accumulation of hyperphosphorylated tau protein into neurofibrillary tangles. The amyloid cascade hypothesis has been the dominant theory of AD pathogenesis since it was first formally proposed in 1992 @hardy_alzheimers_1992 which posits that Aβ deposition is the primary caustive event of AD. The hypothesis has been supported by a large body of evidence over decades and in light of the success of recent Aβ targeting antibody drugs @knopman_lecanemab_2023 the hypothesis has gained renewed strength. John Hardy, author of the 1992 paper, stated "...the trial, which achieved its primary and secondary end points, is proof that reducing brain amyloid has a clear clinical benefit, consistent with the amyloid hypothesis." @hardy_anti-amyloid_2023. Still the amyloid hypothesis has been very controversial and despite the success of anti-amyloid drugs remains unproven and faces critiques that it does not explain the full range of AD pathologies and observations. While the question of whether Aβ is the cause of AD is still debated, it is clear it plays a central role.

#note()[Maybe I need a quick sentence mentioning Ab42 vs Ab40?]

=== Microglia in Alzheimer's Disease <microglia-alzheimers-disease>
In addition to the two classic proteinopathies the disease is also charactrized by synaptic and neuronal loss, neuroinflammation, microgliosis, and astrocytosis. Microglia are the primary innate immune cells of the central nervous system (CNS) and serve many functions in both the developing and adult brain, including phagocytosis, synaptic pruning, neuronal surveillance, antigen presentation, and cytokine release. They are found in all regions of the brain, comprising approximately 5-12% of the total cell population of the adult mouse brain depending on the region(lawson 1990). In his 1907 characterization of the disease Alois Alzheimer described alterations to the glia cells, however, they were not considered to have a central role in the disease until the 1980s when evidence showed that microglia were embedded within Aβ plaques in an "activated" state(Selkoe 1989) and later that they could internalize Aβ @paresce_microglial_1996. These and other observations made microglia a focus of research in AD. 
Historically microglia were often labelled as either 'resting' or 'activated' based on their morphology and expression of surface markers. In the 2000s evidence challenging this simplistic classification emerged, showing that microglia are dynamic and constantly changing in response to their evironment . Classifying microglia as either 'M1' or 'M2' also became popular around this time, but in light of newer evidence has also been shown to be inadequate at describing the heterogeneity of microglia and use is now discouraged @paolicelli_microglia_2022.

In the past 15 years microglia have become the subject of intense research thanks to the advent of new genetic technologies which have suggested a driving role for microglia in AD pathogenesis. In the early 2010s, a series of genome wide association studies (GWAS) identified several genetic risk factors for AD and found many of the high risk genetic variants were in genes expressed in microglia and associated with immune function including TREM2, CD33, ABCA7, MS4A, CR1, CLU, and BIN1 (@harold_genome-wide_2009 @lambert_genome-wide_2009, @seshadri_genome-wide_2010, @hollingworth_common_2011, @naj_common_2011, @guerreiro_using_2013, @jonsson_variant_2013). As new methods were developed for the investigation of the transcriptome of microglia (@butovsky_identification_2014, @bennett_new_2016) insights were gained into mechanisms of the previously identified genetic risk factors (Wang2015). It also allowed for the exploration of patterns of gene expression by age, brain region, and perturbations @grabert_microglial_2016. These analyses began to suggest that microglia were less homogenous than previously thought. The heterogeneity of microlgia, in particular in 'activated' states, became apparent with the advancement of single cell RNA seqeuencing techniques and the publication of the landmark papers by Keren-Shaul et al. @keren-shaul_unique_2017 and Krassmann et al. @krasemann_trem2-apoe_2017. These papers leveraged the transcriptional profile of individual microglia to identify a distinct activation state in mouse models of AD that they termed Disease Associated Microglia (DAM) and microglial neurodegenerative phenotype (MgnD) respectively. These microglia exhibited downregulation of expression of genes associated with homeostatic microglia (Cx3cr1, Tmem119, P2ry12) and upregulation of genes associated with inflammation and phagocytosis (Trem2, Apoe, Cst7, Tyrobp, B2m, Fth1, Lyz2, Axl, Csf1, Lpl, Cd9, Itgax, Clec7a). Additionally Keren-Shaul et al 2017 used Trem2 knockout mice to show that the development of the DAM phenotype was dependent on Trem2. Numerous other studies have now also attempted to identify DAM/MgnD microglia(@hammond_single-cell_2019, @mathys_single-cell_2019, @li_developmental_2019, @masuda_microglia_2020, @sala_frigerio_major_2019, @friedman_diverse_2018, @ellwanger_prior_2021, @millet_exhausted-like_2024, @mancuso_xenografted_2024, @lee_apoe_2023, @grubman_transcriptional_2021). These studies were generally able to identify cells with expression profiles with striking similarity to DAM/MgnD, however it is important to keep in mind that this expression profile has not been formally and rigorously defined leading different groups to use slightly different gene marker lists, different names, and different criteria for defining subpopulations. There are also differences due to use of different mouse models of AD (some which model Aβ pathology, tau pathology, or both), alterations to the genetic background with knockout or knockin of genes, differences in age of the mice, and technical differences in isolation, dissociation, and sequencing. Despite these difficulties, the DAM/MgnD phenotype is now widely accepted as a microglial activation state associated with AD pathology in mouse models. These studies also identified other microglial subtypes with distinct expression profiles again however, it is difficult to compare the profiles. Some patterns of microglia expression have been identified across multiple studies, such as those associated with interferon response genes (Ifit2, irf7, Ifitm3, Cxcl10, Oasl1, etc), MHCII genes (H2-Aa, H2-Ab1, H2-Eb1, H2-DMa, etc), and G2/M phase genes (Top2a, Mki67, Cenpe, Mcm5, Birc5, etc) (@chen_genetically_2021). An issue with the DAM/MgnD expression profile has been the lack of clear equivalent in Human single nucleus RNA sequencing (snRNAseq) studies. Several studies using post mortem human brain tissue of AD patients have now been published, some have identified populations with expression profiles bearing similarities to the DAM/MgnD phenotype or interferon response microglia, however whether these constitute the same populations as those in mouse models is still debated(@olah_single_2020, @n_human_2023, @prater_human_2023, @nguyen_apoe_2020, @tuddenham_cross-disease_2024, @safaiyan_white_2021, @rosenzweig_sex-dependent_2024). The explosion of research in microglial hetergeneity using single cell approaches has led to the identification of many new microglial subtypes, but there is still a lack of consensus on nomenclature and the biological significance of these subtypes. Given the huge amount of data generated by these studies and the differences in methodology considerable effort needs to be made to create a coherent understanding of the microglial heterogeneity in AD.(@reid_microglial_2025, @pettas_profiling_2022). 

#text(fill: red)[(maybe this is all too review-y? perhaps should narrow focus on only the most important recent findings (Keren-Shaul, Krassmann, Olah, and maybe one or two others)?)]


=== Apolipoprotein E <apolipoprotein-e>
Prior to GWAS studies only one gene had been been identified as a risk factor for AD, the apolipoprotein E (APOE) gene. Since its indentification as a risk factor in 1993 by Corder et al @corder_gene_1993 it has remained the strongest genetic risk factor despite dozens of other genes being identified. The APOE gene is located on chromosome 19 and has three common alleles, APOE2, APOE3, and APOE4 which result from the two single nucleotide polymorphisms rs429358 and rs7412 @seripa_genetics_2011. APOE3 is the most common allele accounting for an estimated 78% in the global population, while APOE4 is the second most common at 15% and APOE2 is the rarest at 6%@eisenberg_worldwide_2010. There is a fourth allele, APOE1 (sometimes called APOE3r), but has only been found in a small number of families@seripa_genetics_2011.
The three alleles result in 6 genotypes, 3 heterozygous and 3 homozygous, which produce three different proteoforms: apoE2, apoE3, and apoE4. The proteoforms differ at two amino acid positions, 112 and 158, with the E2 proteoform having cysteine at both positions, E3 having cysteine at position 112 and arginine at position 158, and E4 having arginine at both positions. The protein is 34kDa and composed of 299 amino acids, and is the primary lipid and cholesterol transporter of the CNS. It is primarily produced by astrocytes, but is also produced by microglia, oligodendrocytes, pericytes and neurons, and its expression by microglia increases signifcantly in the presence of AB.
APOE3 is considered to be the "neutral" allele, while APOE2 is protective and APOE4 is deleterious. APOE2 carriers have a 50% reduced risk of developing AD compared to homozygous APOE3 carriers @li_developmental_2019 ,while APOE4 carriers have a 3 fold increased risk for a single allele and 9-15 fold increased risk for homozygous carriers compared to APOE3 carriers. In addition to increased risk of AD, the apoE4 proteoform also induces increased Aβ plaque load and density @schmechel_increased_1993, increased tau pathology @baek_effect_2020, impairs autophagy @parcon_apolipoprotein_2018 @simonovitch_effects_2019 , and reduces brain Aβ clearance @castellano_human_2011. These observations in addition to the upregulation of APOE in microglia associated with Aβ plaques has led to the hypothesis that apoE has a proteoform-specific effect on microglial function, which affets the pathogenesis of AD. A full review of APOE and its interaction with AD can't be provided here but is available elsewhere @fernandez-calle_apoe_2022 .
  
=== 5xFAD and EFAD <5xfad-efad>
Due to the difficulties posed by obtaining and studying human Alzheimer’s disease (AD) brain tissue, transgenic mouse models have become common for investigating disease pathogenesis. One widely used model is the 5xFAD mouse, which was engineered to express five familial AD (FAD)-linked mutations in the amyloid precursor protein (APP) and presenilin-1 (PSEN1) and presenilin-2 (PSEN2) genes @oakley_intraneuronal_2006. These mutations drive an increased production of amyloid-β (Aβ), particularly Aβ42, leading to rapid cerebral Aβ accumulation. By two months of age, 5xFAD mice begin to show amyloid plaques and gliosis, and by six months, they display neuronal loss and impaired memory function, making them a valuable model for studying Aβ-related pathology.

However, 5xFAD mice have several important limitations. Most notably, they do not develop tau-related pathology, the second hallmark of AD, which must be considered when interpreting findings. Additionally, these mice express the murine Apoe gene, which differs significantly from human APOE: mouse Apoe shares only 199 out of 299 amino acids with its human counterpart and exists in a single proteoform @tai_efad_2017. 

To address these issues, 5xFAD mice have been crossed with human APOE knock-in mice to generate the EFAD model, which expresses human APOE2, APOE3, or APOE4 alleles. EFAD mice better replicate the isoform-specific effects observed in humans: for example, E4FAD mice show increased Aβ deposition and greater neuroinflammation compared to E2FAD or E3FAD mice. The introduction of human APOE in place of murine APOE delays the accumulation of amyloid plaques from ~2 months in 5xFAD to ~6 months, which further demonstrates their functional differences @tai_efad_2017. Although the EFAD model improves the study of APOE genotype effects on Aβ pathology, it remains an incomplete model of human AD, as it still lacks tau pathology and does not fully recapitulate the multifactorial aspects of sporadic late-onset AD. Nevertheless, EFAD mice are a critical tool for exploring the interaction between APOE genotype and Aβ pathology under controlled experimental conditions.

#note()[Maybe I need a little paragraph discussing the sex-related differences in AD, and interaction with APOE?]

=== Microglial Priming <microglial-priming>
Microglia have been shown to possess an innate immune memory (IIM), a phenomenon in which exposure to an inflammatory stimulus has long-lasting/sustained effects on the microglial response to subsequent stimuli(NeherCunningham2019). This has been demonstrated in human AD patients and in mouse models of AD. In an amyloid model mouse of AD, microglia injected with Lipopolysaccharide (LPS) at 3 months of age show distinct changes in microgial expression profiles and alter neuropathology @wendeln_innate_2018. Whether the IIM effect is protective or detrimental is unclear, as studies have shown contradictory effects although it appears to be dependent on the timing and dose of the LPS injection. Previous published research by members of the lab that LPS injection at 6 weeks of age in the 5xFAD model lead to improved long term memory function, reduced microglia soma size, and increased microglial internalization of Aβ @yang_lps_2023. (NEED a couple more sentences here... importance of being able to change microglial response given its recently appreciated role in AD pathogenesis, and sense how strong or how this effect is changed based on timing? something bout lecanamb AB clearing and how ppriming could be used to improve this approach? )

The current understanding of the importance of neuroinflammation, its' regulation by microglia, the modulatory effects of APOE genotype, and the priming effect of microglia have informed the design of this study. The aim of this study is to investigate using single cell RNA seq analysis and proteomic analysis of microglia from E3/E4FAD mice how APOE genotype and early life inflammation interact to alter the microglial phenotype. It is part of a larger research project which combines these technologies with mouse behavioral analysis and immunohistochemstry to characterize the effect of these two variables on the amyloid pathology. 

= Materials and methods <materials-and-methods>
== Animals  <animals>
Humanized APOE3 and APOE4 knock-in mice were crossed with the 5xFAD mouse model to generate E3FAD and E4FAD mice, as previously described (Youmans et al., 2012). Wild-type (WT) littermates were used as controls. Mice were housed in a controlled environment (12-hour light/dark cycle, constant temperature and humidity) with free access to food and water. All procedures were conducted in compliance with the European Union laboratory regulations for animal experiments and were approved by the Animal Research Committee of Lund University.

#figure(
  caption:[Summary of Number of Animals by Genotype and Treatment \ #text(weight:"regular")[]],
table(
  stroke: none,
  align: center,
  columns: 3,
  fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },

  table.header(
    [Genotype],
    [Treatment],
    [\# of Animals],
  ),
  table.hline(),
  table.cell(rowspan: 2, align: horizon)[E3WT], [Vehicle],
  [4], [LPS], [3],
  table.hline(stroke: 0.5pt),
  table.cell(rowspan: 2, align: horizon, fill:none)[E3FAD], [Vehicle],
  [4], [LPS], [3],
  table.hline(stroke: 0.5pt),
  table.cell(rowspan: 2, align: horizon)[E4WT], [Vehicle],
  [4], [LPS], [3],
  table.hline(stroke: 0.5pt),
  table.cell(rowspan: 2, align: horizon, fill: none)[E4FAD], [Vehicle],
  [4], [LPS], [5],
  table.hline(stroke: 1pt),
  table.cell(rowspan: 1, colspan: 2, align: horizon)[Total Animals],
  table.cell(rowspan: 1)[30], 
))


== Experimental procedure and sample processing <experimental-procedure>
E3FAD, E4FAD, and their wild-type (WT) counterparts were administered a single intraperitoneal (i.p.) injection of either lipopolysaccharide (LPS, 1 mg/kg) or vehicle (0.9% saline) at postnatal day 9 (P9) (see Table 1 for group assignments). At six months of age, mice were sacrificed by transcardial perfusion. Brains were immediately dissected, and the right hemisphere was processed for microglial cell isolation.

Microglial enrichment was performed using magnetic-activated cell sorting (MACS) with CD11b MicroBeads and LS columns. Enriched microglia were further purified by fluorescence-activated cell sorting (FACS) using an ARIA III flow cytometer. Cells positive for CD45, CX3CR1, and CD11b were collected for downstream analyses.

== Cell Sequencing <cell-sequencing>
From each experimental sample, approximately 2,500 sorted microglial cells were loaded for single-cell capture and cDNA library preparation using the 10x Genomics Chromium Single Cell 3ʹ v3 reagent kit and workflow. Individual samples were barcoded using cell-hashing antibodies and pooled into groups of four samples per sequencing run. Libraries were sequenced on an Illumina platform.

Raw sequencing data were processed with CellRanger (10x Genomics, v6.1.2) using the mouse reference genome GRCm39 (mm39). This pipeline included alignment, filtering, barcode counting, and UMI (unique molecular identifier) counting to generate feature-barcode matrices for each library.

== Quality Control, Normalization, and Data Correction <qc-normalization-dc>
Individual libraries were first filtered to remove cells identified as empty droplets or doublets using the HashSolo algorithm. Additional quality control (QC) filters were applied to each library: cells with fewer than 2,000 total counts, fewer than 300 expressed genes, greater than 5% mitochondrial gene expression, or fewer than 5% ribosomal gene expression were removed. Genes expressed in fewer than three cells were also excluded. QC metrics and visualizations for each sample, both pre- and post-filtering, are included in Supplementary Figure 1 and Supplementary Table 1. 

After filtering, individual libraries were merged into a single dataset using Anndata's concat function with an outer join along the cell axis and a unique merge for the gene axis. Raw counts were saved as a separate layer before normalization. Normalization was performed using Scanpy’s normalize_total function followed by log-transformation with log1p.

Following initial normalization, the merged dataset was examined for the presence of contaminating non-microglial cells. Based on clustering, expression of canonical marker genes, and manual annotation, cells identified as non-microglia were removed. The dataset was re-normalized after removal of non-microglial cells using the same normalization strategy.

Highly variable genes (HVGs) were then identified using Scanpy’s highly_variable_genes function with sample ID as the batch key. These HVGs were used for subsequent dimensionality reduction steps, including principal component analysis (PCA).

The merged dataset was assessed for batch effects using PCA. Several batch correction methods were tested, including Harmony (integration via latent space correction), BBKNN (graph-based nearest neighbor correction), and ComBat (empirical Bayes adjustment). Ultimately, batch effects were not observed to be significant, and no batch correction was applied.

A high contribution of mitochondrial genes to variance across cells was observed during HVG and PCA analysis. To reduce this effect, mitochondrial gene expression was regressed out using Scanpy’s regress_out function prior to downstream analyses.

#figure(
  caption:[Summary of Datasets Before and After Quality Control\ #text(weight:"regular")[]],
table(
  stroke: none,
  align: center,
  columns: 5,
  fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },

  table.header(
    table.cell()[Library],
    table.vline(),
    table.cell(colspan: 2, align: horizon)[Before],
    table.vline(),
    table.cell(colspan: 2, align: horizon)[After],
  ),
    [],
    [\# of Cells],
    [\# of Genes],
    [\# of Cells],
    [\# of Genes],
  table.hline(),
  [D1], [3403], [33989], [2903], [14987],
  [D2], [5041], [33989], [4698], [15514],
  [D3], [5029], [33989], [4081], [15463],
  [D4], [3892], [33989], [3558], [15052],
  [D5], [3825], [33989], [3606], [15260],
  [D6], [4408], [33989], [3603], [15788],
  [D7], [3643], [33989], [3103], [14935],
  [D8], [3903], [33989], [3601], [15926],
  table.hline(stroke: 0.5pt),
  [Merged], [-], [-], [29153], [17341]
))


== Clustering and Annotation <clustering-annotation>

Principal component analysis (PCA) was computed on normalized, log-transformed data. Neighborhood graphs were constructed using the top 25 principal components as input after analyzing the scree plot, and UMAP embeddings were calculated with a min_dist parameter of 0.5 to preserve both local and global structures.

Leiden clustering @traag_louvain_2019 was performed iteratively across a range of resolutions (0.1–1.9 in steps of 0.1). UMAPs colored by cluster identities were visually inspected at each resolution to identify regions of high hetergeneity and monitor for under- and overclustering. Annotation started with highest resolution clustering to capture fine-grained heterogeneity and then merged clusters that showed high similarity in their gene expression profiles. Overclustering was identified when newly created clusters did not significantly differ from neighboring clusters in their top differentially expressed genes (DEGs).

DEGs were identifiedusing the Wilcoxon rank-sum test (rank_genes_groups function). DEGs were calculated for each cluster against the total dataset and neighbouring clusters to combine clusters into with sufficiently different expression profiles. Expression profile annotations were assigned based on the combination of differential expression profiles, cluster relationships, and comparison to canonical gene signatures described in previous studies, including Keren-Shaul et al. @keren-shaul_unique_2017 , Olah et al. @olah_single_2020 , Prater et al. @prater_human_2023 , Mancuso et al. @mancuso_xenografted_2024 , Sala-Frigerio et al. @sala_frigerio_major_2019 , Mathys et al. @mathys_single-cell_2019, Millet et al. @millet_exhausted-like_2024, @hammond_single-cell_2019, @n_human_2023.

Clusters were annotated into 9 expression profiles; Homeostatic, DAM, MHC-II/Antigen Presentation DAM, Cytokine Response, Cycling(G2M/Cdk1+),  BAM-like, Neuronal surveillance/Ribosomal biogenesis, DAM with Ribosomal upregulation, and Immediate Early Gene (IEG) Enriched Homeostatic. Naming of profiles was based on similarity to previously described profiles in the literature, though should not be considered a formal definition. 

== Cluster Proportions <cluster-proportions>
The proportion of each cluster was calculated as the number of cells in the cluster divided by the total number of cells in the sample. A three way ANOVA was performed to assess the effect of genotype (E3 vs E4), treatment (LPS vs Vehicle), and disease status (FAD vs WT) on the proportion of each cluster. An ordinary least squares (OLS) linear regression model was fit to the data with the proportion of each cluster as the dependent variable and the variable main effects and interactions as predictors.  Tukey’s Honest Significant Difference (HSD) test was used for post-hoc pairwise comparisons between conditions. All statistical analyses were performed using the statsmodels package in Python.

== Pseudobulk Analysis <pseudobulk-analysis>

Pseudobulk analysis was performed to more robustly identify DEGs between condition groups as it has been shown to significantly reduce the number of false positives and increase the accuracy of subsequent enrichment analyses as compared to the Wilcoxon rank-sum test @squair_confronting_2021. The Python package 'Decoupler' was used to generate pseudobulk profiles for each individual mouse using the 'sum' method. The pseudobulk profiles were then used to identify DEGs between groups using 'PyDeseq2', a python implementation of the DESeq2 algorithm @muzellec_pydeseq2_2023 @love_moderated_2014. The Wald test was used generate p-values with cook's filtering and correction for multiple testing using Benjamini-Hochberg method. The model was fit with disease staus (FAD vs WT), APOE genotype (E3 vs E4), and treatment (LPS vs Vehicle) as factors with their interactions, and specific comparisons were extracted using the 'results' function. Genes were considered differentially expressed if they had a adjusted p-value < 0.05 and a log2 fold change > |0.25|. The pseudobulk DEGs were then used for downstream analysis including Gene Ontology (GO) enrichment analysis and transcription factor inference.

== Gene Set Enrichment Analysis (GSEA)
Gene Set Enrichment Analysis (GSEA) is a statistical method used to determine whether selected sets of genes show statistically significant differences between two biological states @subramanian_gene_2005. Rather than focusing only on statistically significant DEGs, as other methods do, GSEA evaluates the entire set of genes ranked by a metric which reflects their different expression between conditions, thereby better capturing the pathway-level enrichment differences.

For each contrast of interest, genes were ranked by their Wald test statistic derived from the pseudobulk differential expression analysis. GSEA was then performed using the implementation provided in the Decoupler package @muzellec_pydeseq2_2023, with the MSigDB Hallmark gene sets @a_molecular_2015 as the reference collection. This approach enabled the identification of pathways significantly enriched in up- or down-regulated genes across our conditional contrasts of interest.

= Results <results>

=== Cluster Expression Profiles 

Following dimensionality reduction and Leiden clustering, we identified nine distinct microglial expression profiles based on differentially expressed genes (DEGs) between clusters and reference to canonical expression signatures described in previous studies. Cluster annotation was guided by comparisons to known microglial states in both mouse and human brain studies, including those by Keren-Shaul et al. @keren-shaul_unique_2017, Sala-Frigerio et al. @sala_frigerio_major_2019, Olah et al. @olah_single_2020, Mathys et al. @mathys_single-cell_2019, and others. The final annotated expression profiles were: Homeostatic, DAM, MHC-II/Antigen Presentation DAM, Cytokine Response, Cycling (G2M/Cdk1+), BAM-like, Neuronal Surveillance/Ribosomal Biogenesis, Ribosomal DAM, and Immediate Early Gene (IEG) Enriched Homeostatic.

These expression profiles are visualized in Figure 1. The UMAP plots (Figure 1a) show cluster organization across the entire dataset and by condition. Figure 1b displays the top 10 DEGs for each cluster when compared to the rest, and Figure 1c shows density plots of cluster distributions per condition.

#figure(caption:[ Distinct Expression Clusters and Their Proportions by Condition \ #text(weight:"regular", size: 11pt)[*a.* UMAP visualizations of the dataset colored by manually annotated expression profiles, entire dataset, and by condition. *b.* Heatmap of the top 10 differentially expressed genes for each cluster expression profile showing scaled log-normalized expression values (z-scores) *c.* UMAP density plots for each condition.]],
image("figs/fig1.svg")
)

The Homeostatic cluster was characterized by high expression of canonical microglial maintenance genes, including P2ry12, Tmem119, and Cx3cr1, and served as a reference for identifying more activated or altered states. The IEG cluster shared elevated expression of homeostatic markers but also showed upregulation of immediate early genes such as Fos, Jun, and Egr1, a pattern that has been interpreted as both a stimulus-responsive microglial state and, alternatively, as a dissociation-induced artifact @marsh_dissection_2022  @millet_exhausted-like_2024.

Two small clusters were also identified. The Cycling cluster was defined by elevated expression of cell cycle and proliferation-associated genes such as Stmn1, Top2a, and H2az1, consistent with previously reported populations of proliferative microglia in AD models and aged brain tissue @sala_frigerio_major_2019 @n_human_2023 @ellwanger_prior_2021. The BAM-like (Border-Associated Macrophage-like) cluster was characterized by high expression of perivascular macrophage markers including Pf4, Mrc1, and Ms4a7. Although microglial enrichment was performed during FACS, this population likely reflects a small number of CNS-associated macrophages retained during sorting.

The DAM (Disease-Associated Microglia) cluster showed a transcriptional profile closely matching that originally described by Keren-Shaul et al. @keren-shaul_unique_2017, including high expression of Trem2, Apoe, Cst7, Tyrobp, and Clec7a, accompanied by suppression of homeostatic genes. Three additional clusters were identified with overlapping DAM features but distinct gene expression signatures. The Cytokine Response cluster exhibited elevated expression of inflammatory chemokines and NF-κB target genes such as Ccl3, Ccl4, and Il1b. The MHC-II/Antigen Presentation DAM cluster was enriched for antigen presentation genes including Cd74, H2-Ab1, H2-Aa, and H2-Eb1. The Ribosomal DAM cluster retained DAM-like features with additional upregulation of ribosomal genes.

=== Expression Profile Propotions Across Conditions

After identifying different expression profiles in the dataset, the relative abundance of each cell type was assessed across experimental conditions. The proportions observed in the eight experimental groups are shown in Figure 2a, with absolute cell counts presented in Supplementary Figure 2 and corresponding values listed in Supplementary Table 2.
Disease versus wild-type is the main cause of variation in cell type composition, with the most pronounced differences observed in the homeostatic and DAM clusters. Homeostatic cells accounted for 39–56% of the population in WT mice, compared to only 3–17% in FAD mice, representing a statistically significant reduction (ANOVA p < 1e-8; linear model coefficient = +46.2%, p < 0.001). The APOE genotype was not associated with a statistically significant change in homeostatic cell proportions, nor was there a significant difference with LPS treatment. 

#figure(caption:[ Changes in Cluster Proportions by Condition \ #text(weight:"regular", size: 11pt)[*a.* Stacked Bar plot of the proportion of each cluster in each condition. *b.* Bar plots of the proportion of specific clusters in each condition, significant differences determine by three way anova and Tukey's post-hoc test and indicated by asterisks. ]],
image("figs/fig2.svg")
)

The DAM cluster was absent in all WT groups except for a small proportion (1.9%) in the E3WT LPS group and ranged from 23-52% in FAD mice. Additionally, E4FAD mice have a larger porportion of DAM microglia than E3FAD mice (interaction coefficient = +12.3%, p = 0.012), which is consistent with previous single cell studies. The three clusters representing more "activated" states similar to DAM (Cytokine Response, MHC-II/Antigen Presentation DAM, and Ribosomal DAM), also show zero or very low proportions in WT mice, and are significantly increased in FAD mice.

The Cytokine Response cluster shows an interesting pattern in which the E4FAD LPS group has a significantly lower proportion of cells compared to the E4FAD vehicle group (Tukey p-adj = 0.028), while the E3FAD LPS group is not significantly different from the E3FAD vehicle group (Tukey p-adj = 0.807). This suggests an APOE isoform-specific effect of LPS priming on this microglial population. This interpretation is supported by a statistically significant three-way interaction term in the ANOVA (p = 0.038) and a notable effect size in the linear model (coefficient = −15.1%, p = 0.038), indicating that LPS treatment reduces the proportion of Cytokine Response cells specifically in E4FAD mice.

=== E4FAD but Not E3FAD or WT Mice Exhibit Reduced NfkB-related Expression after LPS Priming

To more robustly identify differentially expressed genes (DEGs) between conditions, pseudobulk analysis was performed on the dataset. This involves summing the expression each gene in all cells from an individual mouse to create a single expression profile, aiding in overcoming the sampling variability in single cell RNA sequencing. The pseudobulk profiles can then be used to identify DEGs with methods typically used for bulk RNA sequencing data, such as DESeq2 which produces fewer false positives than the more common Wilcoxon rank-sum test used in single cell analysis. #note()[Do I need to include justification for pseudobulk analysis? Is this justification more appropriate in the methods section? or discussion? should I include a figure with PCAs/PC metadata association/diagnostic plots?] 
DEGs were identified between the LPS and vehicle treatment groups in the four genetic backgrounds, the results of which are visible in the volcano plots in Figure 3a. There are no more than 1 significant DEG which is shared between any of the comparisons. However, as a general pattern it can be seen that the FAD mice have a much larger number of downregulate genes in the LPS treatment group compared to the vehicle group, while the WT mice show a much smaller number of downregulated genes.
#figure(caption:[Differentially expressed Genes and Pathways in LPS versus Vehicle Treated Mice\ #text(weight:"regular", size: 11pt)[*a.* Volcano plots of pseudobulk DEGs in LPS versus vehicle comparisons in the four genotypes with top genes labelled. *b.* Correlation plot of LogFC values of genes in LPS versus vehicle treatment in E3FAD compared to E4FAD. Pearson's correlation computed on genes above the logFC threshold of 0.25 *c.* GSEA of significant pathways in the LPS versus vehicle comparisons computed on Wald test statistics.]],
image("figs/fig3.svg")
)
When the E3FAD and E4FAD LPS treatment group DEGs are compared, a statistically significant correlation is observed (Pearson's r = -0.44, p = 1.5e-5), indicating that genes that are upregulated in E3FAD mice are likely to be downregulated in E4FAD mice. A number of these genes are those associated with the NfkB pathway, including CCl4, Tlr2, and Nfkbia. 
Based on the pseudobulk DEGs, Gene Set Enrichment Analysis (GSEA) was performed to identify pathways that were significantly enriched in the LPS treatment groups. The results (Fig 3c) show that the top upregulated pathways in both E3WT and E3FAD groups are TNF-a signaling via NfkB, while it is the top downregulated pathway in the E4FAD group. It does not appear in the significant pathways in the E4WT group, however the top upregulated pathway in this group is Interferon gamma response. Additionally, the E3WT and E3FAD groups show signifcant upregulation of the Inflammatory response pathway and TGF beta signaling, while the E4FAD group shows downregulation of the Inflammatory response pathway and no upregulation of inflammatory pathways. 

We also used the pseudobulk DEGs to infer transcription factor activity using the CollecTRI regulon collection. 

#note()[- Was going to add trasncription factor inference here but decided not to include it.
#strike[Add a figure and short section showing the top 6-12 trasncription factors that are upregulated on the basis of the pseudobulk analysis. Can include network plots. Can show that these top transcription factors are different between E3FAD and E4FAD, in being upstream of the NfkB pathway.]]

=== Cytokine Response Microglia

Given that the dampening of NF-κB-related pathways in the E4FAD LPS group coincided with a reduced proportion of Cytokine Response microglia, the transcriptional profile of this cluster was examined to assess its potential contribution to the observed pseudobulk differences.

The Cytokine Response cluster was marked by elevated expression of genes associated with pro-inflammatory signaling and canonical NF-κB pathway activation. These included the chemokines Ccl3 and Ccl4, the interleukins Il1a and Il1b, and the TNF family member Tnf. Several key regulators and targets of NF-κB signaling, including Nfkbia, Nfkbiz, and Tnfaip3, were also upregulated. Additional genes enriched in this cluster were involved in stress response (Atf3, Gadd45b), immune regulation (Cd83, Zfp36, Bcl2a1b), pattern recognition (Tlr2), cytokine signaling (Csf1), and chemokine receptor activity (Ccrl2).

These expression patterns are shown in Figure 4. UMAP projections reveal that these transcripts are highly expressed within the Cytokine Response cluster, with limited expression in other microglial states. Violin plots further illustrate the cluster-specific enrichment of these genes.

#figure(caption:[Genes in the Cytokine Response Cluster \ #text(weight:"regular", size: 11pt)[UMAP plots of the dataset colored by the expression of specific genes which are differentially upregulated in the Cytokine Response cluster and which are part of the NF-κB pathway paired with the violin plots in all clusters.]],
image("figs/fig4.svg")
)

The selective expression of inflammatory genes within this cluster, combined with its reduced abundance in E4FAD mice following LPS treatment, suggests that changes in this microglial state contribute to the largest observed pathway-level differences in the more conservative pseudobulk analysis. 

In addition to the upregulation of pro-inflammatory signaling genes, the Cytokine Response cluster was found to show elevated expression of several immediate early genes (IEGs), including Fos, Jun, Egr1, and Dusp1. This overlap prompted further comparison with the IEG-enriched cluster, which also showed increased expression of these same IEGs relative to the rest of the dataset. Despite this similarity, clear differences were observed between the two clusters. The IEG-enriched cluster retained high expression of canonical homeostatic markers such as Tmem119, Cx3cr1, and P2ry12, which were largely absent in the Cytokine Response cluster. In contrast, the Cytokine Response cluster was enriched for DAM-associated genes including Apoe, Clec7a, Trem2, and Ctsd, both in comparison to the total population and when directly compared to the IEG-enriched cluster.

When compared to the DAM cluster, the Cytokine Response cluster remained distinct. The DAM cluster did not show enrichment for IEGs or NF-κB pathway-associated transcripts, suggesting that the Cytokine Response cluster represents a distinct activation state with overlapping features of both DAM and IEG-enriched microglia. These relationships are visualized in the dot plot in figure 5.

#figure(
  caption: [Comparison of Overlapping and Distinct Features Between Microglial States \ 
  #text(weight: "regular", size: 11pt)[Dot plot showing the expression of immediate early genes (IEGs), homeostatic markers, DAM-associated genes, and NF-κB pathway-related transcripts across the Cytokine Response, IEG-enriched, and DAM clusters. Dot size represents the fraction of cells expressing each gene, calculated using a threshold of log-normalized expression > 0 (equivalent to raw counts ≥ 1). Color indicates average expression level scaled across clusters (Z-score normalized).] ],
  image("figs/fig5.svg")
)

=== REeeallyllly could use a section on the proteomics data

#note()[Hope i have time to add this]

= Discussion <discussion>

=== Microglial Heterogeneity in Alzheimer's Disease

Neuroinflammation has become recognized as a central phenomenon in the pathogenesis of Alzheimer’s disease (AD), both regulating and responding to amyloid-β production, plaque formation, tau hyperphosphorylation, and synaptic and neuronal loss, within a complex and dynamic system that remains incompletely understood @heneka_neuroinflammation_2024. While all central nervous system (CNS) cell types have been implicated, accumulating evidence has positioned microglia—the resident innate immune cells of the brain—as central players in these processes. As research into microglial involvement in AD has expanded, so too has the understanding of their diverse functional states. However, the molecular mechanisms governing transitions between these states remain poorly defined.

Advances in laboratory techniques and in particular single-cell transcriptomics have given us an entirely new view of microglial heterogeneity.  The outdated binary classification of microglia into “resting” or “activated” states has been replaced by an understanding of multiple distinct transcriptional states influenced by disease progression, brain region, age, microenvironment, and genetic background. This shift has been driven by studies describing specific microglial activation profiles, among which the Disease-Associated Microglia (DAM) program, described by Keren-Shaul et al. @keren-shaul_unique_2017, remains the most influential. This transcriptional program was proposed to occur in two stages, dependent on TREM2 signaling, and in response to amyloid deposition. Subsequent studies have identified a broader array of microglial states, including cycling, antigen-presenting, phagocytic, and various pro-inflammatory states. The current challenge in this area of research is integrating a rapidly expanding body of evidence produced with varied model systems, methodologies, nomenclatures and sometimes conflicting interpretations.

APOE genotype remains the strongest genetic risk factor for late-onset Alzheimer’s disease, and in the context of microglial heterogeneity, its importance has become increasingly apparent. APOE is highly expressed by microglia in response to amyloid pathology and has been shown to modulate a number of microglial functions including phagocytosis, lipid metabolism, and inflammatory signaling. The data presented here provide additional evidence for APOE-dependent modulation of microglial state. As expected, homeostatic microglia represented the majority population in wild-type mice and were nearly absent in the presence of amyloid pathology. FAD mice showed a marked shift toward DAM and related activated states, with a significantly higher proportion of DAM microglia observed in E4FAD compared to E3FAD mice. This finding is consistent with previous studies using humanized APOE models.

=== Early-Life Immune Priming Alters Microglia Trajectory

LPS priming of microglia has been shown to influence their function both acutely and over extended periods following the initial insult @perry_microglial_2014. Unpublished data from members of our lab have demonstrated that LPS administration at postnatal day 9 alters the microglial response to amyloid pathology later in life. Changes were observed in the transcriptional levels of specific cytokines, and a microglial population not present in vehicle-treated mice was identified (@yang_lps_2023, preprint). These findings, combined with the well-established role of APOE in modulating microglial state and inflammatory signaling, motivated the present analysis to investigate how human APOE isoforms may interact with early-life immune priming to shape microglial phenotypes in the context of AD-related pathology.(should I add a sentence here saying why post-natal day 9?)

The data presented here provide additional evidence for the impact of APOE genotype on microglial phenotypes and the modulatory effect of LPS priming. As expected, a majority of microglia in non-disease (WT) mice were assigned to the homeostatic cluster, with near-zero representation in the disease-associated transcriptomic profiles. In contrast, FAD mice exhibited a substantial shift, with large proportions of their microglial populations classified into disease-associated clusters and a corresponding reduction in homeostatic cells. An increased proportion of microglia in the DAM cluster was observed in E4FAD mice compared to their E3FAD counterparts, consistent with previous studies using mice expressing human APOE alleles.

=== Genotype-Specific Modulation of Inflammatory Signal

The analysis shows a genotype-specific priming effect on the expression of particular cytokines and inflammatory genes, primarily those in the canonical NF-κB signaling pathway. In mice treated with LPS, a downregulation of these genes was observed specifically in E4FAD animals, while an upregulation was seen in the LPS-treated E3FAD group (Fig. 3). The genes associated with these pathways were not expressed uniformly across the microglial population, indicating that the observed transcriptional changes were not global but instead reflect the behavior of a specific microglial subset. This subset, identified as the Cytokine Response cluster, was not present in non-disease mice, indicating that its emergence is dependent on amyloid pathology. Cluster-level analysis suggests that the observed pathway-level effect in E4FAD mice treated with LPS is driven by a reduction in the proportion of microglia occupying this state.

=== Contextualizing the Cytokine Response Cluster

A Cytokine Response cluster, defined by increased expression of cytokines, chemokines, and NF-κB–related genes, has been observed in other studies. Lee et al. @lee_apoe_2023 used E3FAD and E4FAD mice challenged with LPS and performed single-cell sequencing 24 hours later; a cluster enriched for NF-κB–related transcripts was identified, with higher abundance in E3 compared to E4. Mancuso et al. @mancuso_xenografted_2024 identified a transcriptionally similar population in their xenograft model, significantly enriched in the presence of amyloid pathology, although no APOE genotype–specific differences were reported. Hammond et al. @hammond_single-cell_2019 analyzed healthy mice across a range of ages and identified microglial clusters with similar chemokine and cytokine expression profiles, which increased in abundance with age. This microglial phenotype has now been identified across multiple models and conditions; however, its functional role in Alzheimer’s disease remains unclear.

Beyond the identification of this cluster, cytokine signaling more broadly remains poorly characterized in the context of AD. Cytokines have been proposed to contribute to increased amyloid-β production and accumulation @venegas_microglia-derived_2017, and microglial phagocytosis of amyloid-β has been shown to induce cytokine release @halle_nalp3_2008. Other studies have reported that interleukin-1β reduces microglial clearance of amyloid-β @heneka_nlrp3_2013 and that cytokine production is associated with neuronal loss. However, there is also evidence that LPS-induced inflammation can promote amyloid-β clearance @herber_microglial_2007. These findings highlight the complexity of cytokine signaling in AD, which is governed by context-specific feedback loops that depend critically on the specific cytokines, chemokines, and signaling pathways involved, as well as their timing and intensity

Continued investigation of these microglial populations will be important for understanding their function, how they develop and progress, and how they transition across brain regions and in response to disease states. The analysis presented here suggests that this population can be modulated by pre-symptomatic intervention, which may represent a viable therapeutic strategy. This transcriptomic analysis forms part of a larger study examining the interaction between early-life LPS priming and APOE genotype. Unpublished data from this study indicate improved short-term memory performance in LPS-treated E4FAD mice, an effect not observed in their E3 counterparts. Taken together, these findings suggest that modulating microglial phenotype at an early timepoint may yield beneficial effects on Alzheimer’s disease–related outcomes.

=== Biological vs. Artifactual: Evaluating Immediate Early Gene Expression

One of the subclusters identified in our dataset shows upregulation of immediate early genes (IEGs), most prominently from the Fos and Jun families, while also retaining expression of canonical homeostatic microglial markers. This pattern, referred to here as IEG-enriched, has been previously reported by Marsh et al. @marsh_dissection_2022 who showed that it can be induced by enzymatic dissociation during single-cell processing. They found that the use of mechanical dissociation or the addition of transcriptional inhibitors, such as actinomycin D, during enzymatic dissociation could prevent the emergence of this expression pattern. Based on these findings, they concluded that this IEG-rich profile represents an ex vivo stress response, and advocated for the use of transcriptional inhibitors as a standard practice for microglial single-cell experiments.

However, more recent studies suggest that this signature may not be purely artifactual. Mancuso et al. @mancuso_xenografted_2024 generated a single-cell microglial transcriptomic dataset using the App^NL-G-F amyloid mouse model with transplanted human microglia expressing APOE2, APOE3, or APOE4, and included transcriptional inhibitors specifically to address the ex vivo activation signature reported by Marsh. Despite this, they still identified a cluster they termed “Transitioning Cytokine Response Microglia” (tCRM), which they describe as “…show high levels of homeostatic genes but also express CRM markers.” The top 15 DEGs for this cluster (FOS, JUN, DUSP1, KLF2, HSPA1A1A, IER2, IER3, RHOB, JUNB, Ch25H, HSPA1AB, JUND, FOSB, CEBPD, RGS1) are primarily immediate early genes and closely resemble the ex vivo stress signature identified by Marsh. While Mancuso et al. did not report any APOE isoform-specific differences in this population, they found that it was significantly reduced in APOE knockout mice, suggesting that despite the use of inhibitors, this expression profile may reflect a real APOE-dependent microglial state.

Further evidence for the biological relevance of this type of expression pattern comes from Millet et al. @millet_exhausted-like_2024, who performed scRNA-seq on microglia from E3 and E4 × 5xFAD mice aged to 96 weeks. They identified a population they termed “Terminally Inflammatory Microglia” (TIM), which was enriched in aged animals and which they describe as “marked by concomitant expression of inflammatory genes such as S100a8 and S100a9 and immediate early response genes such as Fos, Jun, and Egr1.” Although their protocol did not include transcriptional inhibitors and they acknowledge resemblance to the ex vivo stress profile reported by Marsh, they argue that TIMs represent a bona fide microglial state. They support this by showing that similar populations are seen in human snRNAseq datasets prepared using both mechanical and enzymatic dissociation, and that the TIM cluster is more prevalent in aged brains and modulated by APOE genotype.

A more recent study by Mulenge et al. @mulenge_transcriptomic_2024 extended the work of Marsh by evaluating the effects of both dissociation and sorting. They confirmed that transcriptional inhibitors can suppress IEG expression but also showed that FACS sorting itself induces an ex vivo activation signature. Specifically, they found that Zfp36, Dusp1, Jun, Ccl4, Ccl3, Tnf, Socs3, and Hspa1a were consistently upregulated across sorted datasets. However, they noted that this expression was generally elevated across the entire dataset and not confined to specific clusters.

In our data, by contrast, these IEG and cytokine-related genes are not uniformly upregulated across the full microglial population but are instead restricted to the Cytokine Response and IEG-enriched clusters. This spatial restriction, along with parallels to expression profiles seen in Mancuso and Millet even in the presence of transcriptional inhibitors, supports the interpretation that these clusters represent true biological microglial states, not merely artifacts of tissue processing.

=== Conclusions and Future Directions

Continued investigation of these microglial populations will be important for understanding their function, how they develop and progress, and how they transition across brain regions and in response to disease states. The analysis presented here suggests that this population can be modulated by pre-symptomatic intervention, which may represent a viable therapeutic strategy. This transcriptomic analysis forms part of a larger study examining the interaction between early-life LPS priming and APOE genotype. As part of this effort, additional data are being generated to validate these transcriptional findings through complementary proteomic analysis using mass spectrometry from the same experimental model. Immunohistochemical staining of brain sections will be used to validate microglial responses and assess amyloid plaque burden, providing a spatial context for the observed molecular changes. These datasets will be integrated with behavioral testing results, which have already indicated improved short-term memory performance in LPS-treated E4FAD mice, an effect not observed in their E3 counterparts. Together, these analyses aim to establish a more comprehensive understanding of how early-life immune events intersect with genetic risk to influence microglial function and disease progression in Alzheimer’s disease.


= Supplemental Information <supplemental-information>
#show figure.where(
  kind: table
): set figure(supplement:[Supporting Table])

#note()[TO PUT IN SUPPLEMENTARY DATA:
- QA plots before & after filtering
  - mito %
  - ribo %
  - \# genes / cell 
  - PCA coloured by library, condition
  - PCA scree plot
- Wilcoxon DEG top genes for each cluster against rest and against closest clusters? probably as a table or attached .xlxs? with logfc/pvals_adj
- Pseudobulk plots
  - PCA plots 
  - metadata pca heatmap
  - 


]

#show table.cell.where(y: 0, ): set text(size: 8pt)
#figure(
  caption:[Detailed Summary of Quality Control Filtering for Each Dataset\ #text(weight:"regular")[]],
table(
  stroke: none,
  align: bottom,
  columns: 11,
  fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },

  table.header(
    table.cell()[Dataset],
    table.cell()[Initial Cells],
    table.cell()[ (Mito)],
    table.cell()[ (Ribo)],
    table.cell()[ (Low Counts)],
    table.cell()[ (Low Genes)],
    table.cell()[Unique Cells Removed],
    table.cell()[Initial Genes],
    table.cell()[Genes Removed],
    table.cell()[Final Cells],
    table.cell()[Final Genes],

  ),
  table.hline(),
  [D1], [5500], [97], [197], [465], [2], [500], [33989], [19002], [5000], [14987],
  [D2], [3683], [85], [150], [400], [3], [450], [33989], [18200], [3233], [15789],
  [D3], [4717], [110], [180], [420], [1], [480], [33989], [17800], [4237], [16189],
  [D4], [3968], [90], [160], [380], [2], [420], [33989], [18700], [3548], [15289],
  [D5], [4220], [95], [170], [390], [1], [440], [33989], [18000], [3780], [15989],
  [D6], [3890], [88], [155], [370], [2], [410], [33989], [18900], [3480], [15089],
  [D7], [4405], [102], [165], [450], [1], [500], [33989], [17600], [3905], [16389],
  [D8], [3600], [80], [145], [360], [2], [400], [33989], [19100], [3200], [14889],
  table.hline(stroke: 0.5pt),
  [Merged], [-], [-], [-], [-], [-], [-], [-], [-], [28643], [17429]
))

#pagebreak()
#show bibliography: set par(leading: 1em, first-line-indent: -0.5em)
#show bibliography: set block(spacing:2em)
#bibliography("./ref.bib", style: "vancouver")




// [Cells Removed (Low Counts)],
// [Cells Removed (Low Genes)],
// [Unique Cells Removed],
// [Initial Genes],
// [Genes Removed],
// [Final Cells],
// [Final Genes],