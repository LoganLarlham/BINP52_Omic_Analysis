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
  Supervisor email: #underline[#lorem(2)] \ Experimental Neuroinflammation Laboratory, Department of Experimental Medical Sciences, Faculty of Medicine, Lund University, Lund, Sweden])
#pagebreak()

#show par: set block(spacing: 2em)
#set par(leading: 2em,  justify: true, first-line-indent: 2em)
#set page(numbering: "1/1", number-align: right,)
#counter(page).update(1) 

= Abstract <abstract>

- Two sentence mentioning AD and APOE
- One/two sentence of importance of microglia in AD
- One sentence emphasizing the characterization of microglia (with scRNAseq) in AD
- One sentence describing aim of study
- One sentence mentioining the exprimental design
- One sentence summarizing the results


= Introduction <introduction>

- One Paragraph on AD
    - low level information of AD
    - Mention of previous understanding of AD 
    - Mention of newer understanding of AD in which microglia play a central role (De strooper?)
- One paragraph on APOE
    - Mention of APOE as the strongest genetic risk factor for AD
    - Discuss previous roles of APOE 

#lorem(50)

#lorem(40)

#lorem(100)

#lorem(50)



= Materials and methods <materials-and-methods>
== Animals  <animals>

#lorem(80)

== Experimental procedure and sample processing <experimental-procedure>

#lorem(50)

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


== Cell Sequencing <cell-sequencing>
#lorem(40)


== Quality Control, Normalization, and Data Correction <qc-normalization-dc>
#lorem(30)

#lorem(30)



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
#lorem(30)

// #figure(
//   caption: [#lorem(1) \ #text(weight:"regular", size: 11pt)[#lorem(15)]],
//   image("Leiden_clusters.png")
// )


#lorem(50)

#lorem(50)

#lorem(50)

#lorem(30)

= Results <results>

#lorem(30)

#lorem(31)



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




= Discussion <discussion>
#lorem(20)

#lorem(20)

#lorem(40)

#lorem(20)

#lorem(50)



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
// #bibliography("project_bib.bib", style: "apa")



