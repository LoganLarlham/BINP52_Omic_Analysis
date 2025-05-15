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

#note()[
- Two sentence mentioning AD and APOE
- One/two sentence of importance of microglia in AD
- One sentence emphasizing the characterization of microglia (with scRNAseq) in AD
- One sentence describing aim of study
- One sentence mentioining the exprimental design
- One sentence summarizing the results]

Alzheimer's disease (AD) is the most common form of dementia, and the most dominant genetic risk factors are the APOE variants. In the past 10 years it has become clear that microglia, brain-resident immune cells, play a central role in AD pathogenesis through their involvement in neuroinflammation and amyloid clearance. Using single-cell RNA sequencing, this study examines how human APOE variants and early-life inflammation modulate microglial responses in a mouse model of amyloid pathology. The results show that the APOE genotype has a significant effect on the microglial response to amyloid pathology, with E4FAD mice showing a higher proportion of DAM microglia than E3FAD mice. The results also show that early life inflammation can lead to long-term changes in microglia that affect their response to amyloid pathology in a genotype-specific manner.

= Introduction/Background <introduction>

#note()[PLOS introduction instructions:
The introduction should:

Provide background that puts the manuscript into context and allows readers outside the field to understand the purpose and significance of the study
Define the problem addressed and why it is important
Include a brief review of the key literature
Note any relevant controversies or disagreements in the field
Conclude with a brief statement of the overall aim of the work and a comment about whether that aim was achieved
]

#note()[
    - One Paragraph low level information of AD, two types Early Onset AD (EOAD) and Late onset (LOAD). LOAD most common form of dementia, higher incidence in females, more severe pathology. 
    - One paragraph 
      - Talk about main hallmarks of AD, amyloid plaques, neurofibrillary tangles, neuroinflammation
      - Sentence describing the long preclinical phase of AD
    - Mention of previous understanding of AD 
    - Mention of newer understanding of AD in which microglia play a central role (De strooper?)]


=== Alzheimer's Disease <alzheimer-disease>
Alzheimer's disease (AD) is a neurodegenerative disorder that is characterized by progressive cognitive decline and memory loss, and the foremost cause of dementia worldwide affecting an estimated 55.2 million people in 2019 @noauthor_global_2021. While historically AD has been defined based on clinical symptoms as well as neuropathological hallmarks usually identified post-mortem, in recent years the field has shifted towards defining AD as a disease based on the presence of the hallmark proteinopathies of amyloid-β (Aβ) plaques and tau neurofibrillary tangles formalized in the National Institute on Aging and Alzheimer's Association (NIA-AA) research framework @jack_jr_nia-aa_2018. There are two main forms of AD known as Early Onset AD (EOAD) and Late Onset AD (LOAD), defined by an arbitrary diagnosis age cutoff of 60 years. It is estimated that 5-10% of AD cases are EOAD whereas LOAD is the much more typical form which accounts for 90-95% of all cases. While EOAD has an estimated heritability of 90-100%, and often caused by known autosomal dominant mutations in the APP, PSEN1, and PSEN2 genes, LOAD is a complex disorder with heritability estimated between 40-80% and is influenced by a combination of genetic, environmental, and lifestyle factors @chen_genetically_2021. Another important aspect of AD is the sex differences in the disease. There are signifcantly more women with AD however, however when accounting for overall higher mortality rates men and age-adjusted prevalence some studies have shown no difference, and some have shown a higher prevalence in women. Additionally there have been demonstrated differences in neuropathologies, Aβ and tau burden, and cognitive decline. The reasons for these differences are not well understood but are thought to be due to a combination of genetic, hormonal, and environmental factors @oveisgharan_sex_2018.


As mentioned the pathological hallmarks of AD are the extracellular deposition and accumulation of aggregated amyloid-β (Aβ) peptides into plaques and the intracellular accumulation of hyperphosphorylated tau protein into neurofibrillary tangles. The amyloid cascade hypothesis has been the dominant theory of AD pathogenesis since it was first formally proposed in 1992 @hardy_alzheimers_1992 which posits that Aβ deposition is the primary caustive event of AD. The hypothesis has been supported by a large body of evidence over decades and in light of the success of recent Aβ targeting antibody drugs @knopman_lecanemab_2023 the hypothesis has gained renewed strength. John Hardy, author of the 1992 paper, stated "...the trial, which achieved its primary and secondary end points, is proof that reducing brain amyloid has a clear clinical benefit, consistent with the amyloid hypothesis." @hardy_anti-amyloid_2023. Still the amyloid hypothesis has been very controversial and despite the success of anti-amyloid drugs remains unproven and faces critiques that it does not explain the full range of AD pathologies and observations. While the question of whether Aβ is the cause of AD is still debated, it is clear it plays a central role.

#note()[Maybe I need a quick sentence mentioning Ab42 vs Ab40?]


#note()[- One/two? paragraph on microglia
    - One to two sentences on microglia basic function 
    - One sentence of previous understanding of micrglial activation
    - Discuss newer understanding of microglial heterogeneity and activation states that has been revealed by scRNAseq, focus on Keren-Shaul et al. 2017 and Krasemen et al. 2017 first. Discuss role of 
    - Mention evidence for general aging microglial heterogeneity (Hajdarovic et al. 2022)
    - Mention of the importance of understanding microglial heterogeneity in AD
    - Discuss other studies that have also identified DAM/MgND. Olah et. al., Tuddenheim et. al., Wang et. al., Hammond et. al., Millet et. al., Mancuso et. al.
    - Sentences about non-DAM microglia subtypes (Wang & Millet & Mancuso)
      - Particular attentiion to Interferon, MHCII,Cytokine/Chemokine expressing microglia
    - Try to focus on genes that are implicated in the identified subtypes.
    - Describe difficulty in validating microglial expression profiles due to differences in human/mouse microglia, autopsy/surgery techniques(olah, Marsh. et al.), mouse amyloid model, age of mice/humans. ]


=== Microglia in Alzheimer's Disease <microglia-alzheimers-disease>
In addition to the two classic proteinopathies the disease is also charactrized by synaptic and neuronal loss, neuroinflammation, microgliosis, and astrocytosis. Microglia are the primary innate immune cells of the central nervous system (CNS) and serve many functions in both the developing and adult brain, including phagocytosis, synaptic pruning, neuronal surveillance, antigen presentation, and cytokine release. They are found in all regions of the brain, comprising approximately 5-12% of the total cell population of the adult mouse brain depending on the region(lawson 1990). #text(fill: red)[(should i add something about microglia/macrophage progenitor, self renewal, and general morphology?)] In his 1907 characterization of the disease Alois Alzheimer described alterations to the glia cells, however, they were not considered to have a central role in the disease until the 1980s when evidence showed that microglia were embedded within Aβ plaques in an "activated" state(Selkoe 1989) and later that they could internalize Aβ (Paresce 1996). These and other observations made microglia a focus of research in AD. 
Historically microglia were often labelled as either 'resting' or 'activated' based on their morphology and expression of surface markers. In the 2000s evidence challenging this simplistic classification emerged, showing that microglia are dynamic and constantly changing in response to their evironment(Nimmerjahn 2005). Classifying microglia as either 'M1' or 'M2' also became popular around this time, but in light of newer evidence has also been shown to be inadequate at describing the heterogeneity of microglia and use is now discouraged(Paolicelli 2022).

In the past 15 years microglia have become the subject of intense research thanks to the advent of new genetic technologies which have suggested a driving role for microglia in AD pathogenesis. In the early 2010s, a series of genome wide association studies (GWAS) identified several genetic risk factors for AD and found many of the high risk genetic variants were in genes expressed in microglia and associated with immune function including TREM2, CD33, ABCA7, MS4A, CR1, CLU, and BIN1 (Harold, Lambert, Seshadri, hollingworth, Naj, Guerrero, Jonsson). As new methods were developed for the investigation of the transcriptome of microglia (butovsky 2014, Bennett 2016, Satoh 2016) insights were gained into mechanisms of the previously identified genetic risk factors (Wang2015). It also allowed for the exploration of patterns of gene expression by age, brain region, and perturbations(Grabert2016). These analyses began to suggest that microglia were less homogenous than previously thought. The heterogeneity of microlgia, in particular in 'activated' states, became apparent with the advancement of single cell RNA seqeuencing techniques and the publication of the landmark papers by Keren-Shaul et al. 2017 and Krassmann et al. 2017. These papers leveraged the transcriptional profile of individual microglia to identify a distinct activation state in mouse models of AD that they termed Disease Associated Microglia (DAM) and microglial neurodegenerative phenotype (MgnD) respectively. These microglia exhibited downregulation of expression of genes associated with homeostatic microglia (Cx3cr1, Tmem119, P2ry12) and upregulation of genes associated with inflammation and phagocytosis (Trem2, Apoe, Cst7, Tyrobp, B2m, Fth1, Lyz2, Axl, Csf1, Lpl, Cd9, Itgax, Clec7a). Additionally Keren-Shaul et al 2017 used Trem2 knockout mice to show that the development of the DAM phenotype was dependent on Trem2. Numerous other studies have now also attempted to identify DAM/MgnD microglia(Hammond 2019, Mathys 2017, Li2019, Masuda2019, Sala frigero 2019, friedman, Ellwanger2020, Millet2024, Mancuso, Lee2023, Grubman2021). These studies were generally able to identify cells with expression profiles with striking similarity to DAM/MgnD, however it is important to keep in mind that this expression profile has not been formally and rigorously defined leading different groups to use slightly different gene marker lists, different names, and different criteria for defining subpopulations. There are also differences due to use of different mouse models of AD (some which model Aβ pathology, tau pathology, or both), alterations to the genetic background with knockout or knockin of genes, differences in age of the mice, and technical differences in isolation, dissociation, and sequencing. Despite these difficulties, the DAM/MgnD phenotype is now widely accepted as a microglial activation state associated with AD pathology in mouse models. These studies also identified other microglial subtypes with distinct expression profiles again however, it is difficult to compare the profiles. Some patterns of microglia expression have been identified across multiple studies, such as those associated with interferon response genes (Ifit2, irf7, Ifitm3, Cxcl10, Oasl1, etc), MHCII genes (H2-Aa, H2-Ab1, H2-Eb1, H2-DMa, etc), and G2/M phase genes (Top2a, Mki67, Cenpe, Mcm5, Birc5, etc) (chen2021). An issue with the DAM/MgnD expression profile has been the lack of clear equivalent in Human single nucleus RNA sequencing (snRNAseq) studies. Several studies using post mortem human brain tissue of AD patients have now been published, some have identified populations with expression profiles bearing similarities to the DAM/MgnD phenotype or interferon response microglia, however whether these constitute the same populations as those in mouse models is still debated(Olah2020, NSun2023, Prater2023, Mathys2019, Nguyen2020, tuddenheim2024, Safaiyan2021, Rosenzweig2024). The explosion of research in microglial hetergeneity using single cell approaches has led to the identification of many new microglial subtypes, but there is still a lack of consensus on nomenclature and the biological significance of these subtypes. Given the huge amount of data generated by these studies and the differences in methodology considerable effort needs to be made to create a coherent understanding of the microglial heterogeneity in AD.(Reid2024, Pettas2022). 

#text(fill: red)[(maybe this is all too review-y? perhaps should narrow focus on only the most important recent findings (Keren-Shaul, Krassmann, Olah, and maybe one or two others)?)]


#note()[- One paragraph on APOE
    - Mention of APOE as the strongest genetic risk factor for AD
    - discuss basic facts, lipid transport, isoforms, history of discovery
    - discuss studies that have shown APOE effect in AD, Astrocytes and microglial; Lee et al., Balu et al.]

=== Apolipoprotein E <apolipoprotein-e>
Prior to GWAS studies only one gene had been been identified as a risk factor for AD, the apolipoprotein E (APOE) gene. Since its indentification as a risk factor in 1993 by Corder et al it has remained the strongest genetic risk factor despite dozens of other genes being identified. The APOE gene is located on chromosome 19 and has three common alleles, APOE2, APOE3, and APOE4 which result from the two single nucleotide polymorphisms rs429358 and rs7412(Seripa2011). APOE3 is the most common allele accounting for an estimated 78% in the global population, while APOE4 is the second most common at 15% and APOE2 is the rarest at 6%(eisenberg2010). There is a fourth allele, APOE1 (sometimes called APOE3r), but has only been found in a small number of families(Seripa2011).
The three alleles result in 6 genotypes, 3 heterozygous and 3 homozygous, which produce three different proteoforms: apoE2, apoE3, and apoE4. The proteoforms differ at two amino acid positions, 112 and 158, with the E2 proteoform having cysteine at both positions, E3 having cysteine at position 112 and arginine at position 158, and E4 having arginine at both positions. The protein is 34kDa and composed of 299 amino acids, and is the primary lipid and cholesterol transporter of the CNS. It is primarily produced by astrocytes, but is also produced by microglia, oligodendrocytes, pericytes and neurons, and its expression by microglia increases signifcantly in the presence of AB.
APOE3 is considered to be the "neutral" allele, while APOE2 is protective and APOE4 is deleterious. APOE2 carriers have a 50% reduced risk of developing AD compared to homozygous APOE3 carriers(Li2020),while APOE4 carriers have a 3 fold increased risk for a single allele and 9-15 fold increased risk for homozygous carriers compared to APOE3 carriers(yamazaki2020). In addition to increased risk of AD, the apoE4 proteoform also induces increased Aβ plaque load and density(Schmechel1993), increased tau pathology(Baek2020), impairs autophagy(Parcon2018, Simonovitch2020), and reduces brain Aβ clearance(Castellano2012). These observations in addition to the upregulation of APOE in microglia associated with Aβ plaques has led to the hypothesis that apoE has a proteoform-specific effect on microglial function, which affets the pathogenesis of AD. A full review of APOE and its interaction with AD can't be provided here but is available elsewhere(fernandez-calle2022).
  
#note()[- Short paragraph intorudcing 5xFADE3/4 model (amyloid model with human APOE), mention mouse apoe.]

=== 5xFAD and EFAD <5xfad-efad>
Due to the difficulties posed by obtaining and studying human Alzheimer’s disease (AD) brain tissue, transgenic mouse models have become common for investigating disease pathogenesis. One widely used model is the 5xFAD mouse, which was engineered to express five familial AD (FAD)-linked mutations in the amyloid precursor protein (APP) and presenilin-1 (PSEN1) and presenilin-2 (PSEN2) genes (Vassar, 2006). These mutations drive an increased production of amyloid-β (Aβ), particularly Aβ42, leading to rapid cerebral Aβ accumulation. By two months of age, 5xFAD mice begin to show amyloid plaques and gliosis, and by six months, they display neuronal loss and impaired memory function, making them a valuable model for studying Aβ-related pathology.

However, 5xFAD mice have several important limitations. Most notably, they do not develop tau-related pathology, the second hallmark of AD, which must be considered when interpreting findings. Additionally, these mice express the murine Apoe gene, which differs significantly from human APOE: mouse Apoe shares only 199 out of 299 amino acids with its human counterpart and exists in a single proteoform (Tai, 2017). 

To address these issues, 5xFAD mice have been crossed with human APOE knock-in mice to generate the EFAD model, which expresses human APOE2, APOE3, or APOE4 alleles. EFAD mice better replicate the isoform-specific effects observed in humans: for example, E4FAD mice show increased Aβ deposition and greater neuroinflammation compared to E2FAD or E3FAD mice. The introduction of human APOE in place of murine APOE delays the accumulation of amyloid plaques from ~2 months in 5xFAD to ~6 months, which further demonstrates their functional differences(Tai, 2017). Although the EFAD model improves the study of APOE genotype effects on Aβ pathology, it remains an incomplete model of human AD, as it still lacks tau pathology and does not fully recapitulate the multifactorial aspects of sporadic late-onset AD. Nevertheless, EFAD mice are a critical tool for exploring the interaction between APOE genotype and Aβ pathology under controlled experimental conditions.


#note()[Maybe I need a little paragraph discussing the sex-related differences in AD, and interaction with APOE?]


=== Microglial Priming <microglial-priming>
Microglia have been shown to possess an innate immune memory (IIM), a phenomenon in which exposure to an inflammatory stimulus has long-lasting/sustained effects on the microglial response to subsequent stimuli(NeherCunningham2019). This has been demonstrated in human AD patients and in mouse models of AD. In an amyloid model mouse of AD, microglia injected with Lipopolysaccharide (LPS) at 3 months of age show distinct changes in microgial expression profiles and alter neuropathology(Wendeln2018). Whether the IIM effect is protective or detrimental is unclear, as studies have shown contradictory effects although it appears to be dependent on the timing and dose of the LPS injection. Previous published research by members of the lab that LPS injection at 6 weeks of age in the 5xFAD model lead to improved long term memory function, reduced microglia soma size, and increased microglial internalization of Aβ(Yang2023). (NEED a couple more sentences here... importance of being able to change microglial response given its recently appreciated role in AD pathogenesis, and sense how strong or how this effect is changed based on timing? something bout lecanamb AB clearing and how ppriming could be used to improve this approach? )


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
#note()[- Need a paragraph which essentially describes this (from wu et al. 2024): 
  - "The specific parameters used for UMAP and clustering were defined after assessment of a wide range of possible
   parameters, which were evaluated in light of cell-state annotations and differential expression results. We 
   start from underclustering and we progressively increase the resolution by identifying further heterogeneity in
    the data, but we prevent overclustering by assessing that high resolutions lead to the definition of extra 
    clusters that do not significantly differ in gene expression from the existing ones. Cell states were annotated
     by means of differential expression (FindAllMarkers() function for overall differential expression, 
     FindMarkers() for side by side comparison)"
- We start with a high resolition and merge clusters that are not significantly different in their top expressed genes. Comparing clusters to each other the entire dataset.
- We also compare the top DEGs to common cluster expression signatures found in  the literature; Keren-shaul, Olah, Prater, Mancuso, Sala-Frigero, Mathys, Millet, Hasselman, Sun N, Askew, give cluster names.]

Principal component analysis (PCA) was computed on normalized, log-transformed data. Neighborhood graphs were constructed using the top 25 principal components as input after analyzing the scree plot, and UMAP embeddings were calculated with a min_dist parameter of 0.5 to preserve both local and global structures.

Leiden clustering (Traag2019) was performed iteratively across a range of resolutions (0.1–1.9 in steps of 0.1). UMAPs colored by cluster identities were visually inspected at each resolution to identify regions of high hetergeneity and monitor for under- and overclustering. Annotation started with highest resolution clustering to capture fine-grained heterogeneity and then merged clusters that showed high similarity in their gene expression profiles. Overclustering was identified when newly created clusters did not significantly differ from neighboring clusters in their top differentially expressed genes (DEGs).

DEGs were identifiedusing the Wilcoxon rank-sum test (rank_genes_groups function). DEGs were calculated for each cluster against the total dataset and neighbouring clusters to combine clusters into with sufficiently different expression profiles. Expression profile annotations were assigned based on the combination of differential expression profiles, cluster relationships, and comparison to canonical gene signatures described in previous studies, including Keren-Shaul et al. (2017), Olah et al. (2020), Prater et al. (2020), Mancuso et al. (2019), Sala-Frigerio et al. (2019), Mathys et al. (2019), Millet et al. (2024), Hasselmann et al. (2019), Sun et al. (2021).

Clusters were annotated into 9 expression profiles; Homeostatic, DAM, MHC-II/Antigen Presentation DAM, Cytokine Response, Cycling(G2M/Cdk1+),  BAM-like, Neuronal surveillance/Ribosomal biogenesis, DAM with Ribosomal upregulation, and Immediate Early Gene (IEG) Upregulated. Naming of profiles was based on similarity to previously described profiles in the literature, though should not be considered a formal definition. 

== Cluster Proportions <cluster-proportions>
The proportion of each cluster was calculated as the number of cells in the cluster divided by the total number of cells in the sample. A three way ANOVA was performed to assess the effect of genotype (E3 vs E4), treatment (LPS vs Vehicle), and disease status (FAD vs WT) on the proportion of each cluster. An ordinary least squares (OLS) linear regression model was fit to the data with the proportion of each cluster as the dependent variable and the variable main effects and interactions as predictors.  Tukey’s Honest Significant Difference (HSD) test was used for post-hoc pairwise comparisons between conditions. All statistical analyses were performed using the statsmodels package in Python.

== Pseudobulk Analysis <pseudobulk-analysis>
#note()[- Describe Psuedobulk principal, considering expression together of all cells from an individual.
- Subset the data to specific genotypes (E3WT, E3FAD, E4WT, E4FAD) and use Decoupler to generate expression profiles for each individual using 'sum' method.
-   Deseq2 to identify differentially expressed genes between groups.
This method is recognized as more robust and accurate than the more common wilcoxon rank-sum method of identifying DEGs in scRNAseq data. Squair J.W., et al 2021]

Pseudobulk analysis was performed to more robustly identify DEGs between condition groups as it has been shown to significantly reduce the number of false positives and increase the accuracy of subsequent enrichment analyses as compared to the Wilcoxon rank-sum test (Squair et al., 2021). The Python package 'Decoupler' was used to generate pseudobulk profiles for each individual mouse using the 'sum' method. The pseudobulk profiles were then used to identify DEGs between groups using 'PyDeseq2', a python implementation of the DESeq2 algorithm(Muzellec2023, Love2014). The Wald test was used generate p-values with cook's filtering and correction for multiple testing using Benjamini-Hochberg method. The model was fit with disease staus (FAD vs WT), APOE genotype (E3 vs E4), and treatment (LPS vs Vehicle) as factors with their interactions, and specific comparisons were extracted using the 'results' function. Genes were considered differentially expressed if they had a adjusted p-value < 0.05 and a log2 fold change > |0.25|. The pseudobulk DEGs were then used for downstream analysis including Gene Ontology (GO) enrichment analysis and transcription factor inference.

== Transcription Factor Inference
Transcription factor activity was inferred using the Decoupler package (Muzellec et al., 2023).  For each LPS versus vehicle contrast, the DESeq2 Wald test statistic for every gene was formatted as a 1 × G vector.  The CollecTRI regulon collection (Müller-Dött et al., 2023) containing 1,100 transcription factors (TF) with signed TF - target genes was used in conjunction with the test statistics in a ULM (Univariate Linear Model) to infer TF activity in a pseudobulked condition compared to another. 

== Gene Set Enrichment Analysis (GSEA)
Gene Set Enrichment Analysis (GSEA) is a statistical method used to determine whether selected sets of genes show statistically significant differences between two biological states (Subramanian et al., 2005). Rather than focusing only on statistically significant DEGs, as other methods do, GSEA evaluates the entire set of genes ranked by a metric which reflects their different expression between conditions, thereby better capturing the pathway-level enrichment differences.

For each contrast of interest, genes were ranked by their Wald test statistic derived from the pseudobulk differential expression analysis. GSEA was then performed using the implementation provided in the Decoupler package (Muzellec et al., 2023), with the MSigDB Hallmark gene sets (Liberzon et al., 2015) as the reference collection. This approach enabled the identification of pathways significantly enriched in up- or down-regulated genes across our conditional contrasts of interest.



= Results <results>



=== Cluster Expression Profiles 
// #note()[- This section should simply describe the clustering and the expression profiles, the genes they were based on, and these patterns seen in the literature a bit. No discussion of the differences between conditions.] 

Following dimensionality reduction and Leiden clustering, we identified nine distinct microglial expression profiles based on differentially expressed genes (DEGs) between clusters and reference to canonical expression signatures described in previous studies. Cluster annotation was guided by comparisons to known microglial states in both mouse and human brain studies, including those by Keren-Shaul et al. (2017), Sala-Frigerio et al. (2019), Olah et al. (2020), Mathys et al. (2019), and others. The final annotated expression profiles were: Homeostatic, DAM, MHC-II/Antigen Presentation DAM, Cytokine Response, Cycling (G2M/Cdk1+), BAM-like, Neuronal Surveillance/Ribosomal Biogenesis, Ribosomal DAM, and Immediate Early Gene (IEG) Upregulated.

These expression profiles are visualized in Figure 1. The UMAP plots (Figure 1a) show cluster organization across the entire dataset and by condition. Figure 1b displays the top 10 DEGs for each cluster when compared to the rest, and Figure 1c shows density plots of cluster distributions per condition.

#figure(caption:[ Distinct Expression Clusters and Their Proportions by Condition \ #text(weight:"regular", size: 11pt)[*a.* UMAP visualizations of the dataset colored by manually annotated expression profiles, entire dataset, and by condition. *b.* Heatmap of the top 10 differentially expressed genes for each cluster expression profile showing scaled log-normalized expression values (z-scores) *c.* UMAP density plots for each condition.]],
image("figs/fig1.svg")
)

The Homeostatic cluster was characterized by high expression of canonical microglial maintenance genes, including P2ry12, Tmem119, and Cx3cr1, and served as a reference for identifying more activated or altered states. The IEG cluster shared elevated expression of homeostatic markers but also showed upregulation of immediate early genes such as Fos, Jun, and Egr1, a pattern that has been interpreted as both a stimulus-responsive microglial state and, alternatively, as a dissociation-induced artifact (Marsh et al., 2022; Millet 2024).

Two small clusters were also identified. The Cycling cluster was defined by elevated expression of cell cycle and proliferation-associated genes such as Stmn1, Top2a, and H2az1, consistent with previously reported populations of proliferative microglia in AD models and aged brain tissue (Sala-Frigerio et al., 2019; Sun et al., 2021; Ellwanger et al., 2020). The BAM-like (Border-Associated Macrophage-like) cluster was characterized by high expression of perivascular macrophage markers including Pf4, Mrc1, and Ms4a7. Although microglial enrichment was performed during FACS, this population likely reflects a small number of CNS-associated macrophages retained during sorting.

The DAM (Disease-Associated Microglia) cluster showed a transcriptional profile closely matching that originally described by Keren-Shaul et al. (2017), including high expression of Trem2, Apoe, Cst7, Tyrobp, and Clec7a, accompanied by suppression of homeostatic genes. Three additional clusters were identified with overlapping DAM features but distinct gene expression signatures. The Cytokine Response cluster exhibited elevated expression of inflammatory chemokines and NF-κB target genes such as Ccl3, Ccl4, and Il1b. The MHC-II/Antigen Presentation DAM cluster was enriched for antigen presentation genes including Cd74, H2-Ab1, H2-Aa, and H2-Eb1. The Ribosomal DAM cluster retained DAM-like features with additional upregulation of ribosomal genes.





=== Expression Profile Propotions Across Conditions
// #note()[- This Section should describe the differences in the proportions of the clusters between conditions. First mentioning the FAD/WT difference between homeostatic and DAM microglia. Mention Higher DAM in E4FAD than E3FAD. Lower homeostatic in E4WT than E3WT and E4FAD than E3FAD. MHCII cluster is higher in E4FAD than E3FAD. FAD have cytokine response clusters, and seems in E4FAD LPS seems to cause a reduction.]

After identifying different expression profiles in the dataset, the relative abundance of each cell type was assessed across experimental conditions. The proportions observed in the eight experimental groups are shown in Figure 2a, with absolute cell counts presented in Supplementary Figure 2 and corresponding values listed in Supplementary Table 2.
Disease versus wild-type is the main cause of variation in cell type composition, with the most pronounced differences observed in the homeostatic and DAM clusters. Homeostatic cells accounted for 39–56% of the population in WT mice, compared to only 3–17% in FAD mice, representing a statistically significant reduction (ANOVA p < 1e-8; linear model coefficient = +46.2%, p < 0.001). The APOE genotype was not associated with a statistically significant change in homeostatic cell proportions, nor was there a significant difference with LPS treatment. 

#figure(caption:[ Changes in Cluster Proportions by Condition \ #text(weight:"regular", size: 11pt)[*a.* Stacked Bar plot of the proportion of each cluster in each condition. *b.* Bar plots of the proportion of specific clusters in each condition, significant differences determine by three way anova and Tukey's post-hoc test and indicated by asterisks. ]],
image("figs/fig2.svg")
)

The DAM cluster was absent in all WT groups except for a small proportion (1.9%) in the E3WT LPS group and ranged from 23-52% in FAD mice. Additionally, E4FAD mice have a larger porportion of DAM microglia than E3FAD mice (interaction coefficient = +12.3%, p = 0.012), which is consistent with previous single cell studies. The three clusters representing more "activated" states similar to DAM (Cytokine Response, MHC-II/Antigen Presentation DAM, and Ribosomal DAM), also show zero or very low proportions in WT mice, and are significantly increased in FAD mice.

The Cytokine Response cluster shows an interesting pattern in which the E4FAD LPS group has a significantly lower proportion of cells compared to the E4FAD vehicle group (Tukey p-adj = 0.028), while the E3FAD LPS group is not significantly different from the E3FAD vehicle group (Tukey p-adj = 0.807). This suggests an APOE isoform-specific effect of LPS priming on this microglial population. This interpretation is supported by a statistically significant three-way interaction term in the ANOVA (p = 0.038) and a notable effect size in the linear model (coefficient = −15.1%, p = 0.038), indicating that LPS treatment reduces the proportion of Cytokine Response cells specifically in E4FAD mice.


=== E4FAD but Not E3FAD or WT Mice Exhibit Reduced NfkB-related Expression after LPS Priming

#note()[-  There is a reduction in the expression of pro-inflammatory andeet cytokine-related genes in E4FAD mice after LPS treatment. will talk about the fact that E4FAD has been identified in the past to have higher expression of these genes(hammond?). But the LPS priming reduces the expression of these genes. This reduction is not seen in the E3FAD mice. To support this result show the pseudbulk volcanos highlighting genes that are relevant. And also show significant difference between specific genes related to these pathways. ]

To more robustly identify differentially expressed genes (DEGs) between conditions, pseudobulk analysis was performed on the dataset. This involves summing the expression each gene in all cells from an individual mouse to create a single expression profile, aiding in overcoming the sampling variability in single cell RNA sequencing. The pseudobulk profiles can then be used to identify DEGs with methods typically used for bulk RNA sequencing data, such as DESeq2 which produces fewer false positives than the more common Wilcoxon rank-sum test used in single cell analysis. #note()[Do I need to include justification for pseudobulk analysis? Is this justification more appropriate in the methods section? or discussion? should I include a figure with PCAs/PC metadata association/diagnostic plots?] 
DEGs were identified between the LPS and vehicle treatment groups in the four genetic backgrounds, the results of which are visible in the volcano plots in Figure 3a. There are no more than 1 significant DEG which is shared between any of the comparisons. However, as a general pattern it can be seen that the FAD mice have a much larger number of downregulate genes in the LPS treatment group compared to the vehicle group, while the WT mice show a much smaller number of downregulated genes.
#figure(caption:[Differentially expressed Genes and Pathways in LPS versus Vehicle Treated Mice\ #text(weight:"regular", size: 11pt)[*a.* Volcano plots of pseudobulk DEGs in LPS versus vehicle comparisons in the four genotypes with top genes labelled. *b.* Correlation plot of LogFC values of genes in LPS versus vehicle treatment in E3FAD compared to E4FAD. Pearson's correlation computed on genes above the logFC threshold of 0.25 *c.* GSEA of significant pathways in the LPS versus vehicle comparisons computed on Wald test statistics.]],
image("figs/fig3.svg")
)
When the E3FAD and E4FAD LPS treatment group DEGs are compared, a statistically significant correlation is observed (Pearson's r = -0.44, p = 1.5e-5), indicating that genes that are upregulated in E3FAD mice are likely to be downregulated in E4FAD mice. A number of these genes are those associated with the NfkB pathway, including CCl4, Tlr2, and Nfkbia. 
Based on the pseudobulk DEGs, Gene Set Enrichment Analysis (GSEA) was performed to identify pathways that were significantly enriched in the LPS treatment groups. The results (Fig 3c) show that the top upregulated pathways in both E3WT and E3FAD groups are TNF-a signaling via NfkB, while it is the top downregulated pathway in the E4FAD group. It does not appear in the significant pathways in the E4WT group, however the top upregulated pathway in this group is Interferon gamma response. Additionally, the E3WT and E3FAD groups show signifcant upregulation of the Inflammatory response pathway and TGF beta signaling, while the E4FAD group shows downregulation of the Inflammatory response pathway and no upregulation of inflammatory pathways. 

We also used the pseudobulk DEGs to infer transcription factor activity using the CollecTRI regulon collection. 

#note()[- Add a figure and short section showing the top 6-12 trasncription factors that are upregulated on the basis of the pseudobulk analysis. Can include network plots. Can show that these top transcription factors are different between E3FAD and E4FAD, in being upstream of the NfkB pathway.]

// #figure(caption:[ Transcription Factor Analysis of Treatment effects in E3FAD and E4FAD mice \ #text(weight:"regular", size: 11pt)[*A.* Bar plot of the top 10 differentially expressed transcription factors in E3FAD mice. *B.* Bar plot of the top 10 differentially expressed transcription factors in E4FAD mice. *C.* Network plot of the top 5 transcription factors and their target genes in E3FAD mice. *D.* Network plot of the top 5 transcription factors and their target genes in E4FAD mice. ]],
// image("figs/fig4.svg")
// )

=== Cytokine Response Microglia
#note()[- Take a closer look at the expression of Cytokine Response Microglia, specifically look at the expression of specific genes.
- Talk about the chemokines and cytokines in the inflammatory response, whats unique about this cluster, the higher expression of Il1, TNF-a, NfkB related genes.  
- Maybe I put this section after the proportions section?]

Given that the dampening of NF-κB-related pathways in the E4FAD LPS group coincided with a reduced proportion of Cytokine Response microglia, the transcriptional profile of this cluster was examined to assess its potential contribution to the observed pseudobulk differences.

The Cytokine Response cluster was marked by elevated expression of genes associated with pro-inflammatory signaling and canonical NF-κB pathway activation. These included the chemokines Ccl3 and Ccl4, the interleukins Il1a and Il1b, and the TNF family member Tnf. Several key regulators and targets of NF-κB signaling, including Nfkbia, Nfkbiz, and Tnfaip3, were also upregulated. Additional genes enriched in this cluster were involved in stress response (Atf3, Gadd45b), immune regulation (Cd83, Zfp36, Bcl2a1b), pattern recognition (Tlr2), cytokine signaling (Csf1), and chemokine receptor activity (Ccrl2).

These expression patterns are shown in Figure 4. UMAP projections reveal that these transcripts are highly expressed within the Cytokine Response cluster, with limited expression in other microglial states. Violin plots further illustrate the cluster-specific enrichment of these genes.

The selective expression of inflammatory genes within this cluster, combined with its reduced abundance in E4FAD mice following LPS treatment, suggests that changes in this microglial state contribute to the largest observed pathway-level differences in the more conservative pseudobulk analysis. 

#figure(caption:[Genes in the Cytokine Response Cluster \ #text(weight:"regular", size: 11pt)[UMAP plots of the dataset colored by the expression of specific genes which are differentially upregulated in the Cytokine Response cluster and which are part of the NF-κB pathway paired with the violin plots in all clusters.]],
image("figs/fig4.svg")
)

In addition to the upregulation of pro-inflammatory signaling genes, the Cytokine Response cluster was found to show elevated expression of several immediate early genes (IEGs), including Fos, Jun, Egr1, and Dusp1. This overlap prompted further comparison with the IEG-enriched cluster, which also showed increased expression of these same IEGs relative to the rest of the dataset. Despite this similarity, clear differences were observed between the two clusters. The IEG-enriched cluster retained high expression of canonical homeostatic markers such as Tmem119, Cx3cr1, and P2ry12, which were largely absent in the Cytokine Response cluster. In contrast, the Cytokine Response cluster was enriched for DAM-associated genes including Apoe, Clec7a, Trem2, and Ctsd, both in comparison to the total population and when directly compared to the IEG-enriched cluster.

When compared to the DAM cluster, the Cytokine Response cluster remained distinct. The DAM cluster did not show enrichment for IEGs or NF-κB pathway-associated transcripts, suggesting that the Cytokine Response cluster represents a distinct activation state with overlapping features of both DAM and IEG-enriched microglia. These relationships are visualized in the dot plot in figure 5.

#figure(
  caption: [Comparison of Overlapping and Distinct Features Between Microglial States \ 
  #text(weight: "regular", size: 11pt)[Dot plot showing the expression of immediate early genes (IEGs), homeostatic markers, DAM-associated genes, and NF-κB pathway-related transcripts across the Cytokine Response, IEG-enriched, and DAM clusters. Dot size represents the fraction of cells expressing each gene, calculated using a threshold of log-normalized expression > 0 (equivalent to raw counts ≥ 1). Color indicates average expression level scaled across clusters (Z-score normalized).] ],
  image("figs/fig5.svg")
)

=== REeeallyllly could use a section on the proteomics data

x
= Discussion <discussion>

#note()[
  - Discuss Neuroinflammation and the impact on AD pathogenesis, microglial heterogeneity its importance and lack of clarity on mechanisms, transitions/tmeporal changes (dynamic) and biological effect.
    - 
  - short mention of APOE genotype biggest risk factor of ad and the seeming importance related to microglia heterogeinety and function (phagocytosis) + lps priming microglial/immune modulation could be a part of therapies if we can better figure out how it happens mechanistically and effects of priming.

  - Mention that this data is consistent with previous studies showing that FAD mice have low population of homeostatic microglia + that E4FAD mice have a higher proportion of DAM microglia.
    -
  - Quick talk about how DAM is thought to be an endstate of from homeo > DAM and there is in betwen stages (Dam 1 and DAM 2 from the original Keren-Shaul)

  - Mention the connection to our labs previous (unpublished) work showing that LPS treatment at P9 leads to (insert exact differences here) in (non-humanized APOE) 5xFAD mice.
    - "Unpublished data from other members from our lab which gave LPS treatment at postnatal day 9 to 5xFAD mice found changes in transcriptional levels of particular cytokines and a population cluster of microglia not seen in vehicle treated mice. these findings prompted this analysis to wonder how human APOE isoforms may modulate this type of effect. 

  - discuss that lee 2024 used an acute LPS treatment in E3/E4FAD and found cytokine response up, and that hammond also find CR microglia in aged mice, mancuso finds cytokine repsonse cluster in human microglia transplated into mice.
    - It is not clear whether these microglia are damaging or protective and how they might influence uptak, synapse loss, Ab uptake and plaque deposition.

  - Discuss that the data (appears to be consisent with/could support) as of yet unpublished data from the same mice showing that LPS treatment at P9 showed an improvment in long term memory function in E4FAD mice but not E3FAD mice.
    -
  - Discuss that if disregulated neuroinflammation/increase of proinflammatory microglia (seen in vehicle EFADs) is detrmintal, reduced expression of pro-inflammatory genes in E4FAD mice after LPS treatment could be protect (in line with the above memory function data)
    -
  - Discuss IEG population in our dataset (including that IEGs are transcription factores in inflammatory responses), and the elavated level in Cytokine response mg. Then mention the Marsh et al. (+olah+Li2020) implications on millet (+mancuso) interpretation. Add in Gern2024, not just enzymatic tissue dissociation but also FACS sorting, not just IEGs but also cytokine response genes.
    - I guess at somepoint, either here or in the "Cytokine Response Microglia" results section I should write about the fact that it is also expressing IEG genes, and maybe compare and contrast to IEG-enriched and DAM its two closest populations. 

    - "Briefly mentioned in the results, there is a subcluster within the dataset called IEG enriched which shows upregulation of the immediate early genes, primarily Jun and Fos family which also has similar expression of the 'homeostatic' gene signature as well. this expression pattern in microglia has been identified previously in Marsh et al 2022 as being caused by ex vivo stress caused by the use of enzymatic dissociaiton during the processing of tissue for single cell. They find that this expression pattern can be avoided by using mechanical dissociation or by adding transcription inhibtors during the enzymatic dissociation procedure, and advocate for the addition of inhibitors (actinomycin D) for single cell transcriptomics experiments. Mancuso generated a single cell microglia transcriptomic dataset using the App^Nl-G-F amyloid mouse model with human microglia with APOE 2,3,4 transplanted and included the inhibitors specifically to address the results reported in marsh, however they report a cluster called Transitioning Cytokine Response Microglia (tCRM) which they describe as "...show high levels of homeostatic genes but also express CRM markers". The top 15 DEGs for this cluster are FOS, JUN, DUSP1, KLF2, HSP1A1A, IER2, IER3, RHOB, JUNB, Ch25H, HSP1AB, JUND, FOSB, CEBPD, RGS1, primarily immediate early genes and closely matching the ex vivo stress signature identified by Marsh, despite their use of inhibitors. They do not report any differences in this cluster between, the APOE isoforms, but that it is significantly reduced in APOE knockout mice. This suggests that despite the opportunity for enzymatic dissociation to alter the expression this type of signature represents a true biological profile. Further evidence for a similar expression pattern is seen in Millet et al 2024 who using E3/E4 x 5xFAD mice aged to 96 weeks find a cluster enriched in the 96 week mice they name "Terminally Inflammatory Microglia" (TIM) which they describe as "Marked by concomitant expression of inflammatory genes such as S100a8and S100a9and immediate early response genes such as Fos, Jun, and Egr1". In their murine samples they do not use inhibitors and acknowledge resemblance to the ex vivo stress profile described by Marsh, but justify TIMs as a true microglial state by the fact that they find a similar population in human snRNAseq microglial data produced both with mechanical and enzymatic dissociation and further that aged brains are significantly more likely to showing this microglial state and it is modulated by APOE genotype, suggesting it is not solely ex vivo stress. 
    A further study published by Mulenge 2024 builds upon the work described by Marsh showing that the use of transcriptional inhibitors reduced these activation signatures Furthermore they report an important finding that FACS sorting to enrich for microglia also induces an ex vivo activation signature with elevated expression of Immediate Early genes and they identify Zfp36, Dusp1, Jun, Ccl4, Ccl3, Tnf, Socs3, Hspa1a as being consistently upregulated across sorted datasets. However they note that it is a consistent upregulation seen across the entire dataset, not in particular clusters. As most microglial single cell transcriptomics studies use cell sorting to increase microglial yield from samples including Millet and Mancuso's experimental procedure and this analysis use sorted microglia caution should be used when interpreting the levels of these markers in particular. However, in this analysis the upregulation was not consistent across the entire dataset and the relevant genes can be seen to be uniquely upregulated in the Cytokine Response cluster, suggesting that this signal is not simply due to cell sorting stress. 

  - Maybe Discuss Malat1 presence and interpretations? Li et al. 2024/Cell ranger FAQ (actually probably not important)
]

Neuroinflammation has become recognized as a central phenomenon in the pathogenesis of AD, both regulating and responding to the Ab production, plaque formation, p-tau deposition, neuronal and synapse loss, in a complex and dynamic system which we have yet to understand fully (heneka2024). While it is clear all cells of the CNS are involved, evidence has accumulated highlighting microglia, the primary innate immune cells of the brain, as key to all these facets. As the amount of research into microglia's role in AD has exploded so to has our appreciation for their many functions and regulation of these function though our understanding of mechanisms of transition remain unclear. 

Advances in laboratory techniques and in particular single-cell transcriptomics have given us an entirely new view of microglial heterogeneity. Our understanding has evolved from an outdated model of binary resting and activation states to many distinct activation states influenced by the disease state, brain region, age, local environment, and genome. This has been evidenced by many studies describing specific transcriptional programs, perhaps most influential among these the DAM described by Keren-Shaul in 2017. Described as a response to amyloid pathology following a two step activation process dependent on signalling by TREM2. Subsequent studies have built upon this finding identifying a host of different microglial profiles which affect their cycling microglia, antigen presenting microglia, phagocytic microglia, and various inflammatory subtypes. The current challenge in this area of research is integrating a rapidly expanding body of evidence produced with varied model systems, methodologies, nomenclatures and sometimes conflicting interpretations.

It has long been known that APOE is the strongest genetic risk factor for late-onset AD, in the context of microglial diversity it is becoming even more significant. Evidence is pointing to it being a key mediator of microglial phenotype in response to AB. It is highly expressed by microglial in response to amyloid pathology and has been shown to modulate phagocytosis, lipid metabolism, and inflammatory signalling (source here). 


LPS priming of microglia has been shown to have an affect on microglial function both in short time periods and long after the insult(Perry2014). Unpublished data from members of our lab have demonstrated an effect of LPS priming at post natal day 9 in microglial response to amyloid pathology. They observed changes in transcriptional levels of particular cytokines and a population cluster of microglia not seen in vehicle treated mice(preprint; yang2024). these findings combined with the known effect of APOE prompted this analysis to investigate how human APOE isoforms may modulate this type of effect.

The data presented here provide additional evidence of the impact of APOE genotype on microglial phenotypes and the influence of LPS priming. As expected there is a plurality of microglia in the non-disease mice in the homeostatic cluster and near-zero in the disease associated transcriptomic profile. The FAD mice however have large roportions of their populations in the disease associated clusters, and in increase in the related activated states. Consistent with previous studies with mice with humanized APOE, an increase in the proportion of the DAM cluster was observed in the APOE4 mice verus their E3 counterparts. 

The analysis shows a genotype specific priming effect on the expression of particular cytokines and inflammatory genes, primarily those in the canonical track of the NF-κB signaling pathway. In mice who were treated with LPS there is a downregulation of these genes (fig3) specifically in E4FAD mice while an upregulation was seen in LPS treated E3FAD mice. The genes related to these pathways are not expressed at consistent levels across the population is not due to a global change, rather a specific population of microglia express these genes. This population is not seen in the non-disease mice, indicating it is a response to amyloid pathology. The analysis of the clusters suggests that the pathway level effect is driven by a reduction of the proportion of microglia in the Cytokine Response cluster in the LPS treated E4FAD mice. 

A Cytokine Response cluster, or a population with increased cytokine, chemokine, and NfkB-related genes has been observed in other studies. Lee2024 used E3FAD and E4FAD mice challenged with LPS and single cell sequencing performed 24 hours after and observed a cluster with increased expression NfκB-related genes, which was enriched in E3 versus E4. Hammond2024 finds a cluster defined by these same genes in their model, which is signifcantly enriched by amyloid beta, but no significant differences between APOE genotypes. Hammond2019 analysed healthy mice at various ages and also found clusters of microglia with similar chemokine and cytokine expression, and found an increase of this population with age. This microglial phenotype appears across a variety of models and conditions, however it is not yet known what affect they might have on surrounding cells, or their role in AD as protective or damaging.

It will be important to continue toresearch these types of microglia and their function, how they develop, progress, and transition within different regions and in response to disease states. THe above analysis shows that this population can be modulated by pre-symptomatic treatment which could provide a possible avenue for therapy. This transcriptomic analysis is part of a larger study of the lps priming and APOE interaction. Unpublished data from the study indicate improved short term memory outcomes in behaviour tests of LPS treated E4FAD mice not seen in the E3 counterparts. Combined with the analysis above indicates that modulating microglial phenotype at an early timepoint could have beneficial effects on the symptoms of AD. 


One of the subclusters identified in our dataset shows upregulation of immediate early genes (IEGs), most prominently from the Fos and Jun families, while also retaining expression of canonical homeostatic microglial markers. This pattern, referred to here as IEG-enriched, has been previously reported by Marsh et al. (2022) who showed that it can be induced by enzymatic dissociation during single-cell processing. They found that the use of mechanical dissociation or the addition of transcriptional inhibitors, such as actinomycin D, during enzymatic dissociation could prevent the emergence of this expression pattern. Based on these findings, they concluded that this IEG-rich profile represents an ex vivo stress response, and advocated for the use of transcriptional inhibitors as a standard practice for microglial single-cell experiments.

However, more recent studies suggest that this signature may not be purely artifactual. Mancuso et al. generated a single-cell microglial transcriptomic dataset using the App^NL-G-F amyloid mouse model with transplanted human microglia expressing APOE2, APOE3, or APOE4, and included transcriptional inhibitors specifically to address the ex vivo activation signature reported by Marsh. Despite this, they still identified a cluster they termed “Transitioning Cytokine Response Microglia” (tCRM), which they describe as “…show high levels of homeostatic genes but also express CRM markers.” The top 15 DEGs for this cluster (FOS, JUN, DUSP1, KLF2, HSPA1A1A, IER2, IER3, RHOB, JUNB, Ch25H, HSPA1AB, JUND, FOSB, CEBPD, RGS1) are primarily immediate early genes and closely resemble the ex vivo stress signature identified by Marsh. While Mancuso et al. did not report any APOE isoform-specific differences in this population, they found that it was significantly reduced in APOE knockout mice, suggesting that despite the use of inhibitors, this expression profile may reflect a real APOE-dependent microglial state.

Further evidence for the biological relevance of this type of expression pattern comes from Millet et al. (2024), who performed scRNA-seq on microglia from E3 and E4 × 5xFAD mice aged to 96 weeks. They identified a population they termed “Terminally Inflammatory Microglia” (TIM), which was enriched in aged animals and which they describe as “marked by concomitant expression of inflammatory genes such as S100a8 and S100a9 and immediate early response genes such as Fos, Jun, and Egr1.” Although their protocol did not include transcriptional inhibitors and they acknowledge resemblance to the ex vivo stress profile reported by Marsh, they argue that TIMs represent a bona fide microglial state. They support this by showing that similar populations are seen in human snRNAseq datasets prepared using both mechanical and enzymatic dissociation, and that the TIM cluster is more prevalent in aged brains and modulated by APOE genotype.

A more recent study by Mulenge et al. (2024) extended the work of Marsh by evaluating the effects of both dissociation and sorting. They confirmed that transcriptional inhibitors can suppress IEG expression but also showed that FACS sorting itself induces an ex vivo activation signature. Specifically, they found that Zfp36, Dusp1, Jun, Ccl4, Ccl3, Tnf, Socs3, and Hspa1a were consistently upregulated across sorted datasets. However, they noted that this expression was generally elevated across the entire dataset and not confined to specific clusters.

In our data, by contrast, these IEG and cytokine-related genes are not uniformly upregulated across the full microglial population but are instead restricted to the Cytokine Response and IEG-enriched clusters. This spatial restriction, along with parallels to expression profiles seen in Mancuso and Millet even in the presence of transcriptional inhibitors, supports the interpretation that these clusters represent true biological microglial states, not merely artifacts of tissue processing.


= Limitations and Future Direction

#note()[Discuss limitations of the study.
- 5xFAD model is not a perfect model of AD.
- Only females
- Only 1 time point(?)
- Possible effects due to tissue processing as seen in Marsh
Future directions.
- We know microglia are heterogenous across regions, spatial transcriptomic methods could be used to see differences in areas of the brain.
- 
- ]

= Supporting information <supporting-information>
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