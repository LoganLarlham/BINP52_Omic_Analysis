#import "@local/wordometer:0.1.5": word-count, total-words
#show: word-count.with(exclude: heading)

#set text(font: "TeX Gyre Termes", size: 12pt)
#show heading.where(level: 1): set text(size: 18pt)
#show heading.where(level: 2): set text(size: 16pt)
#show figure.caption: set text(weight: "bold")
#show figure.caption: set par(leading: 1em)
#show figure.caption: set align(left)
#show <n>: set text(fill: red)



#show figure.where(
  kind: table
): set block(breakable: true)
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


#v(20pt)
Alzheimer's disease (AD) is a neurodegenerative disorder that is characterized by progressive cognitive decline and memory loss, and the foremost cause of dementia worldwide affecting an estimated 55.2 million people in 2019 @noauthor_global_2021. While historically AD has been defined based on clinical symptoms as well as neuropathological hallmarks usually identified post-mortem, in recent years the field has shifted towards defining AD as a disease based on the presence of the hallmark proteinopathies of amyloid-β (Aβ) plaques and tau neurofibrillary tangles formalized in the National Institute on Aging and Alzheimer's Association (NIA-AA) research framework @jack_jr_nia-aa_2018. There are two main forms of AD known as Early Onset AD (EOAD) and Late Onset AD (LOAD), defined by an arbitrary diagnosis age cutoff of 60 years. It is estimated that 5-10% of AD cases are EOAD whereas LOAD is the much more typical form which accounts for 90-95% of all cases. While EOAD has an estimated heritability of 90-100%, and often caused by known autosomal dominant mutations in the APP, PSEN1, and PSEN2 genes, LOAD is a complex disorder with heritability estimated between 40-80% and is influenced by a combination of genetic, environmental, and lifestyle factors @chen_genetically_2021. Another important aspect of AD is the sex differences in the disease. There are signifcantly more women with AD however, however when accounting for overall higher mortality rates men and age-adjusted prevalence some studies have shown no difference, and some have shown a higher prevalence in women. Additionally there have been demonstrated differences in neuropathologies, Aβ and tau burden, and cognitive decline. The reasons for these differences are not well understood but are thought to be due to a combination of genetic, hormonal, and environmental factors @oveisgharan_sex_2018.


As mentioned the pathological hallmarks of AD are the extracellular deposition and accumulation of aggregated amyloid-β (Aβ) peptides into plaques and the intracellular accumulation of hyperphosphorylated tau protein into neurofibrillary tangles. The amyloid cascade hypothesis has been the dominant theory of AD pathogenesis since it was first formally proposed in 1992 @hardy_alzheimers_1992 which posits that Aβ deposition is the primary caustive event of AD. The hypothesis has been supported by a large body of evidence over decades and in light of the success of recent Aβ targeting antibody drugs @knopman_lecanemab_2023 the hypothesis has gained renewed strength. John Hardy, author of the 1992 paper, stated "...the trial, which achieved its primary and secondary end points, is proof that reducing brain amyloid has a clear clinical benefit, consistent with the amyloid hypothesis." @hardy_anti-amyloid_2023. Still the amyloid hypothesis has been very controversial and despite the success of anti-amyloid drugs remains unproven and faces critiques that it does not explain the full range of AD pathologies and observations. While the question of whether Aβ is the cause of AD is still debated, it is clear it plays a central role.


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


In addition to the two classic proteinopathies the disease is also charactrized by synaptic and neuronal loss, neuroinflammation, microgliosis, and astrocytosis. Microglia are the primary innate immune cells of the central nervous system (CNS) and serve many functions in both the developing and adult brain, including phagocytosis, synaptic pruning, neuronal surveillance, antigen presentation, and cytokine release. They are found in all regions of the brain, comprising approximately 5-12% of the total cell population of the adult mouse brain depending on the region(lawson 1990). #text(fill: red)[(should i add something about microglia/macrophage progenitor, self renewal, and general morphology?)] In his 1907 characterization of the disease Alois Alzheimer described alterations to the glia cells, however, they were not considered to have a central role in the disease until the 1980s when evidence showed that microglia were embedded within Aβ plaques in an "activated" state(Selkoe 1989) and later that they could internalize Aβ (Paresce 1996). These and other observations made microglia a focus of research in AD. 
Historically microglia were often labelled as either 'resting' or 'activated' based on their morphology and expression of surface markers. In the 2000s evidence challenging this simplistic classification emerged, showing that microglia are dynamic and constantly changing in response to their evironment(Nimmerjahn 2005). Classifying microglia as either 'M1' or 'M2' also became popular around this time, but in light of newer evidence has also been shown to be inadequate at describing the heterogeneity of microglia and use is now discouraged(Paolicelli 2022).

In the past 15 years microglia have become the subject of intense research thanks to the advent of new genetic technologies which have suggested a driving role for microglia in AD pathogenesis. In the early 2010s, a series of genome wide association studies (GWAS) identified several genetic risk factors for AD and found many of the high risk genetic variants were in genes expressed in microglia and associated with immune function including TREM2, CD33, ABCA7, MS4A, CR1, CLU, and BIN1 (Harold, Lambert, Seshadri, hollingworth, Naj, Guerrero, Jonsson). As new methods were developed for the investigation of the transcriptome of microglia (butovsky 2014, Bennett 2016, Satoh 2016) insights were gained into mechanisms of the previously identified genetic risk factors (Wang2015). It also allowed for the exploration of patterns of gene expression by age, brain region, and perturbations(Grabert2016). These analyses began to suggest that microglia were less homogenous than previously thought. The heterogeneity of microlgia, in particular in 'activated' states, became apparent with the advancement of single cell RNA seqeuencing techniques and the publication of the landmark papers by Keren-Shaul et al. 2017 and Krassmann et al. 2017. These papers leveraged the transcriptional profile of individual microglia to identify a distinct activation state in mouse models of AD that they termed Disease Associated Microglia (DAM) and microglial neurodegenerative phenotype (MgnD) respectively. These microglia exhibited downregulation of expression of genes associated with homeostatic microglia (Cx3cr1, Tmem119, P2ry12) and upregulation of genes associated with inflammation and phagocytosis (Trem2, Apoe, Cst7, Tyrobp, B2m, Fth1, Lyz2, Axl, Csf1, Lpl, Cd9, Itgax, Clec7a). Additionally Keren-Shaul et al 2017 used Trem2 knockout mice to show that the development of the DAM phenotype was dependent on Trem2. Numerous other studies have now also attempted to identify DAM/MgnD microglia(Hammond 2019, Mathys 2017, Li2019, Masuda2019, Sala frigero 2019, friedman, Ellwanger2020, Millet2024, Mancuso, Lee2023, Grubman2021). These studies were generally able to identify cells with expression profiles with striking similarity to DAM/MgnD, however it is important to keep in mind that this expression profile has not been formally and rigorously defined leading different groups to use slightly different gene marker lists, different names, and different criteria for defining subpopulations. There are also differences due to use of different mouse models of AD (some which model Aβ pathology, tau pathology, or both), alterations to the genetic background with knockout or knockin of genes, differences in age of the mice, and technical differences in isolation, dissociation, and sequencing. Despite these difficulties, the DAM/MgnD phenotype is now widely accepted as a microglial activation state associated with AD pathology in mouse models. These studies also identified other microglial subtypes with distinct expression profiles again however, it is difficult to compare the profiles. Some patterns of microglia expression have been identified across multiple studies, such as those associated with interferon response genes (Ifit2, irf7, Ifitm3, Cxcl10, Oasl1, etc), MHCII genes (H2-Aa, H2-Ab1, H2-Eb1, H2-DMa, etc), and G2/M phase genes (Top2a, Mki67, Cenpe, Mcm5, Birc5, etc) (chen2021). An issue with the DAM/MgnD expression profile has been the lack of clear equivalent in Human single nucleus RNA sequencing (snRNAseq) studies. Several studies using post mortem human brain tissue of AD patients have now been published, some have identified populations with expression profiles bearing similarities to the DAM/MgnD phenotype or interferon response microglia, however whether these constitute the same populations as those in mouse models is still debated(Olah2020, NSun2023, Prater2023, Mathys2019, Nguyen2020, tuddenheim2024, Safaiyan2021, Rosenzweig2024). The explosion of research in microglial hetergeneity using single cell approaches has led to the identification of many new microglial subtypes, but there is still a lack of consensus on nomenclature and the biological significance of these subtypes. Given the huge amount of data generated by these studies and the differences in methodology considerable effort needs to be made to create a coherent understanding of the microglial heterogeneity in AD.(Reid2024, Pettas2022). 

#text(fill: red)[(maybe this is all too review-y? perhaps should narrow focus on only the most important recent findings (Keren-Shaul, Krassmann, Olah, and maybe one or two others)?)]


#note()[- One paragraph on APOE
    - Mention of APOE as the strongest genetic risk factor for AD
    - discuss basic facts, lipid transport, isoforms, history of discovery
    - discuss studies that have shown APOE effect in AD, Astrocytes and microglial; Lee et al., Balu et al.]


Prior to GWAS studies only one gene had been been identified as a risk factor for AD, the apolipoprotein E (APOE) gene. Since its indentification as a risk factor in 1993 by Corder et al it has remained the strongest genetic risk factor despite dozens of other genes being identified. The APOE gene is located on chromosome 19 and has three common alleles, APOE2, APOE3, and APOE4 which result from the two single nucleotide polymorphisms rs429358 and rs7412(Seripa2011). APOE3 is the most common allele accounting for an estimated 78% in the global population, while APOE4 is the second most common at 15% and APOE2 is the rarest at 6%(eisenberg2010). There is a fourth allele, APOE1 (sometimes called APOE3r), but has only been found in a small number of families(Seripa2011).
The three alleles result in 6 genotypes, 3 heterozygous and 3 homozygous, which produce three different proteoforms: apoE2, apoE3, and apoE4. The proteoforms differ at two amino acid positions, 112 and 158, with the E2 proteoform having cysteine at both positions, E3 having cysteine at position 112 and arginine at position 158, and E4 having arginine at both positions. The protein is 34kDa and composed of 299 amino acids, and is the primary lipid and cholesterol transporter of the CNS. It is primarily produced by astrocytes, but is also produced by microglia, oligodendrocytes, pericytes and neurons. 
APOE3 is considered to be the "neutral" allele, while APOE2 is protective and APOE4 is deleterious. APOE2 carriers have a 50% reduced risk of developing AD compared to homozygous APOE3 carriers(Li2020),while APOE4 carriers have a 3 fold increased risk for a single allele and 9-15 fold increased risk for homozygous carriers compared to APOE3 carriers(yamazaki2020). In addition to increased risk of AD, the apoE4 proteoform also induces increased Aβ plaque load and density(Schmechel1993), increased tau pathology(Baek2020), impairs autophagy(Parcon2018, Simonovitch2020), and reduces brain Aβ clearance(Castellano2012). These observations in addition to the upregulation of APOE in microglia associated with Aβ plaques has led to the belief that apoE proteforms effect on microglial function is a key mechanism in AD pathogenesis. A full review of APOE and its interaction with AD can't be provided here but is available elsewhere(fernandez-calle2022).
  
#note()[- Short paragraph intorudcing 5xFADE3/4 model (amyloid model with human APOE), mention mouse apoe.]



#note()[- SHort Paragraph about LPS activation 
    - LPS activation of microglia as a model for neuroinflammation
    - Mention Unpublished data from lab that early life inflammation can lead to long term changes in microglia
    - This study aims to this early life inflammation model to study the effects of APOE genotype on the microglial response to AD (amyloid) pathology.]
- discuss/justify use of only female mice; more severe phenotype in model, higher prevalence of AD in humans, and Wu et al. showing female specific clusters.



= Materials and methods <materials-and-methods>
== Animals  <animals>
- Describe the mouse model, number of samples, number per Group, animal care details

#figure(
  caption:[Summary of Genotype Treatment Combinations \ #text(weight:"regular")[]],
table(
  stroke: none,
  align: center,
  columns: 3,
  fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },

  table.header(
    [#lorem(1)],
    [#lorem(1)],
    [#lorem(1)],
  ),
  table.hline(),
  table.cell(rowspan: 2, align: horizon)[#lorem(1)], [#lorem(1)],
  [4], [#lorem(1)], [5],
  table.hline(stroke: 0.5pt),
  table.cell(rowspan: 2, align: horizon, fill:none)[#lorem(1)], [#lorem(1)],
  [4], [#lorem(1)], [5],
  table.hline(stroke: 0.5pt),
  table.cell(rowspan: 2, align: horizon)[#lorem(1)], [#lorem(1)],
  [4], [LPS], [3],
  table.hline(stroke: 0.5pt),
  table.cell(rowspan: 2, align: horizon, fill: none)[#lorem(1)], [#lorem(1)],
  [5], [LPS], [5],
  table.hline(stroke: 0.5pt),
))


== Experimental procedure and sample processing <experimental-procedure>
- Describe LPS injection
- Describe Tissue collection (brain dissection)
- Describe tissue processing, i.e. method of dissociation, microglial isolation
- Mention FACS sorting
== Cell Sequencing <cell-sequencing>
- single cell library preparation
- Mention of 10x genomics platform and chemistry used
- Mention of sequencing depth/ quality metrics of raw data
- Mention of Cell Ranger pipeline/ Alignment to mm10
- Hashing/HTO of samples and use of Hashsolo to demultiplex


== Quality Control, Normalization, and Data Correction <qc-normalization-dc>
- Individual libraries are filtered to remove cells that are found to be empty or doublets by Hashsolo
- Individual libraries are filtered to remove cells that have less than 2000 reads, more than 5% mitochondrial reads(find the paper which mentions this), less than 5% ribosomal reads, cells with fewer than 300 genes expressed, and genes with fewer than 3 cells expressing them.
- Create supplemetary tables of QC metrics; before and after
- Libraries are merged and then normalized using the Scanpy normalize_total and log1p functions.
- Merged dataset is analyzed and based on expression of canonical marker genes and clustering, any infilitrating non-microglia cells are removed from the dataset.
- the dataset is re-normalized after removal of non-microglia cells.
- Highly Variable Genes (HVGs) are identified using the Scanpy highly_variable_genes function and are used for principal component analysis.
- Merged dataset is analyzed for batch effects by principal component analysis and various batch correction methods are tested; Harmony, BBKNN, and ComBat which represent 3 different approaches to batch correction.
- High variability of mitochondrial genes was identified and the data was regressed for mitochondrial gene expression.


#figure(
  caption:[Summary of Datasets Before and After Quality Control\ #text(weight:"regular")[]],
table(
  stroke: none,
  align: center,
  columns: 5,
  fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },

  
  table.header(
    table.cell()[Dataset],
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
  [D1], [5500], [33989], [5357], [15849],
  [D2], [3683], [33989], [3572], [15311], 
  [D3], [4717], [33989], [4614], [16152],
  [D4], [3968], [33989], [3819], [15238],
  table.hline(stroke: 0.5pt),
  [Merged],[-],[-],[17362],[17089]
))

#lorem(30)


== Clustering and Annotation <clustering-annotation>
- Need a paragraph which essentially describes this (from wu et al. 2024): 
  - "The specific parameters used for UMAP and clustering were defined after assessment of a wide range of possible
   parameters, which were evaluated in light of cell-state annotations and differential expression results. We 
   start from underclustering and we progressively increase the resolution by identifying further heterogeneity in
    the data, but we prevent overclustering by assessing that high resolutions lead to the definition of extra 
    clusters that do not significantly differ in gene expression from the existing ones. Cell states were annotated
     by means of differential expression (FindAllMarkers() function for overall differential expression, 
     FindMarkers() for side by side comparison)"
- We start with a high resolition and merge clusters that are not significantly different in their top expressed genes. 
- We also compare the top DEGs to common cluster expression signatures found in  the literature; Keren-shaul, Olah, Prater, Mancuso, Sala-Frigero, Mathys, Millet, Hasselman, Sun N, Askew give cluster names.


== Pseudobulk Analysis <pseudobulk-analysis>
- Describe Psuedobulk principal, considering expression together of all cells from an individual.
- Subset the data to specific genotypes (E3WT, E3FAD, E4WT, E4FAD) and use Decoupler to generate expression profiles for each individual using 'sum' method.
-   Deseq2 to identify differentially expressed genes between groups.
This method is recognized as more robust and accurate than the more common wilcoxon rank-sum method of identifying DEGs in scRNAseq data. Squair J.W., et al 2021

== Transcription Factor Inference
- Small description of what transcription factors are and how they regulate gene expression
- How they can be inferred from scRNAseq data (Univariate linear model)
- Use of Decoupler on pseudobulk DEGs for VEH to LPS comparisons to infer TFs

== Gene Ontology Analysis
- Description of Gene ontology analysis and how it can be used to identify enriched biological processes
- Gprofiler tool used with GO databases. 
- Statistical method of Gprofiler
- Using both pseudobulk DEGs a and using Wilcoxon rank-sum DEGs to identify enriched biological processes
- Also find enriched biological process based on (wilcoxon) DEGs by cluster(v rest) 

== Statistical difference testig of individual gene expression 
- Ddescribe the use of two way anova to identify gene clusters which are proportionally different between groups
- Use of two way anova two find differences of individual gene mean expression between groups using n=3-5 mice per group



// #figure(
//   caption: [#lorem(1) \ #text(weight:"regular", size: 11pt)[#lorem(15)]],
//   image("Leiden_clusters.png")
// )




// #figure(caption: [#lorem(1) \ #text(weight:"regular", size: 11pt)[#lorem(10)]],
// image("Sample_Genotype.png")
// )

#lorem(30)



// #figure(caption: [#lorem(1) \ #text(weight:"regular")[#lorem(10)]],
// image("Decoupler_pred.png")
// )

// #figure(caption: [#lorem(1)s \ #text(weight:"regular", size: 11pt)[#lorem(10).]],
// align(image("Marker_genes.png")
// ))

#lorem(20)
#lorem(20)

// #figure(caption: [Heatmaps of Genes Expression Grouped by Leiden Clusters \ #text(weight:"regular", size: 11pt)[Expression levels are log-transformed normalized values. Clusters are generated from the Leiden algorithm with resolution 1.4.]],
// scale(x:105%, y:105%, reflow: true)[#image("Ribo_heatmap.png")],
// )


When the analysis of canonical marker genes, top differentially expressed genes, and results of the automated annotation by Decoupler and ACT were taken together the following annotated map of the dataset was produced (Figure 6.). The map includes 6 cell types which have been constructed by combining clusters identified by the Leiden algorithm, starting with 20 individual clusters identified with a resolution of 1.4. Clusters were combined on the basis of shared expression profiles. 

// #figure(caption: [Expression levels of Marker genes \ #text(weight:"regular", size: 11pt)[UMAP visualization of data clustered into manually annotated clusters representative of gene expression profiles. Legend includes top DEGs for each cluster compared to rest of dataset.]],
// image("Annotated_Umap.png")
// )


#lorem(20)

// #figure(caption:[Differentially Expressed Genes between DAM and High ribosomal DAM \ #text(weight:"regular", size: 11pt)[Top differentially expressed genes by Wilcoxon rank-sum statistical test between DAM and High ribosomal DAM with all other cell types excluded to identify the most different genes between them]],
// image("DAM_rank_genes.png")
// )

#lorem(20)

// #figure(caption:[ Bar plot of Cell Cycle Phases of cell populations\ #text(weight:"regular", size: 11pt)[Cell populations are manually annotated as in Fig 6. ]],
// image("Cellcycle_plot.png")
// )




= Results <results>
- As has been previously shown the 5xFAD model in mice does show a "DAM" microglial cluster which is characterized by high expression of Apoe, Trem2, Cst7, and other genes.
- E4FAD mice show a higher proportion of DAM microglia than E3FAD mice. regardless of LPS treatment.
- found 11 different clusters with different genetic expression profiles.
- The largest differences between cluster proportions are found between the disease and non disease groups. 
- Within the disease groups, the largest differences are found between the E3FAD and E4FAD groups.
- Profiles with an enriched cytokine response (ccl2/3/4 expression) and interferon response itfml3 etc. are found in the E4FAD group.
- Identificaiton of an MHCII expressing cluster, discuss similarity to BAM expresion profile but notabole lack of specific genes  like MRc1, Pf4, Stab1, Lyz2, Lyve1 (Sun Jiang 2024)
- Reduction of BAM in all FAD mice. 


= Discussion <discussion>

#note()[
  - Discuss the Marsh et al. (+olah) implications on millet interpretation.
  - Discuss Malat1 interpretations? Li et al. 2024
]
= Limitations and Future Direction

#lorem(80)

#lorem(40)

#lorem(20)

#lorem(50)

= Supporting information <supporting-information>
#show figure.where(
  kind: table
): set figure(supplement:[Supporting Table])

#show figure.where(
  kind: image
): set figure(supplement:[Supporting Figure])

#counter(figure.where(kind: table)).update(0)
#figure(
  caption:[#lorem(6)\ #text(weight:"regular")[#lorem(6)]],
table(
  stroke: none,
  columns: 3,
  fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },

   table.header(
    table.cell()[Cluster],
    table.vline(),
    table.cell(colspan: 1, align: horizon)[Cell Type],
    table.vline(),
    table.cell(colspan: 1, align: horizon)[Number of Variable Genes],
  ),
  table.hline(),
  [Group_0], [Monocyte-derived dendritic cell], [5],
  table.hline(),
  [Group_0], [Microglial cell], [15],
  table.hline(),
  [Group_1], [Microglial cell], [5],
  table.hline(),
  [Group_1], [Homeostatic microglial cell], [15],
  table.hline(),
  [Group_10], [NA], [5],
  table.hline(),
  [Group_10], [Dendritic cell], [15],
  table.hline(),
  [Group_11], [Microglial cell], [5],
  table.hline(),
  [Group_11], [Microglial cell], [15],
  table.hline(),
  [Group_12], [Microglial cell], [5],
  table.hline(),
  [Group_12], [Microglial cell], [15],
  table.hline(),
  [Group_13], [Tissue-resident macrophage], [5],
  table.hline(),
  [Group_13], [Tissue-resident macrophage], [15],
  table.hline(),
  [Group_14], [Tissue-resident macrophage], [5],
  table.hline(),
  [Group_14], [Tissue-resident macrophage], [15],
  table.hline(),
  [Group_15], [Microglial cell], [5],
  table.hline(),
  [Group_15], [Astrocyte], [15],
  table.hline(),
  [Group_16], [Monocyte-derived dendritic cell], [5],
  table.hline(),
  [Group_16], [Monocyte], [15],
  table.hline(),
  [Group_17], [Microglial cell], [5],
  table.hline(),
  [Group_17], [Microglial cell], [15],
  table.hline(),
  [Group_18], [Epithelial cell], [5],
  table.hline(),
  [Group_18], [Microglial cell], [15],
  table.hline(),
  [Group_2], [Homeostatic microglial cell], [5],
  table.hline(),
  [Group_2], [Homeostatic microglial cell], [15],
  table.hline(),
  [Group_3], [Homeostatic microglial cell], [5],
  table.hline(),
  [Group_3], [Homeostatic microglial cell], [15],
  table.hline(),
  [Group_4], [Microglial cell], [5],
  table.hline(),
  [Group_4], [Microglial cell], [15],
  table.hline(),
  [Group_5], [NA], [5],
  table.hline(),
  [Group_5], [T cell], [15],
  table.hline(),
  [Group_6], [Microglial cell], [5],
  table.hline(),
  [Group_6], [Microglial cell], [15],
  table.hline(),
  [Group_7], [Microglial cell], [5],
  table.hline(),
  [Group_7], [Microglial cell], [15],
  table.hline(),
  [Group_8], [NA], [5],
  table.hline(),
  [Group_8], [NA], [15],
  table.hline(),
  [Group_9], [Cardiac muscle cell], [5],
  table.hline(),
  [Group_9], [Megakaryocyte], [15],
))


#pagebreak()
#show bibliography: set par(leading: 1em, first-line-indent: -0.5em)
#show bibliography: set block(spacing:2em)
#bibliography("./ref.bib", style: "vancouver")




