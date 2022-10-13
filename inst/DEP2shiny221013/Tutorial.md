# **Welcome to the tutorial**

<hr>

### **Contents**

#### 1. [Overview of the analysis](#Overview)

#### 2. [Differential protein expression for label-free quantification (DEP-LFQ) or other quantification](#DEP-LFQ)

####   2.1 [Overview](#overview_DEP) {#overview}

####   2.2 [Analytical steps](#analytical_steps_DEP)

####   2.3 [Output](#output_DEP)

####    2.3.1 [Result tables](#output_DEP_result_tables)

####    2.3.2 [Visualization of the results](#output_DEP_result_plot)

####     2.3.2.1 [Barplot of a protein of interest](#output_DEP_result_plot_single_protein)

####     2.3.2.2 [Heatmap of all significant proteins](#output_DEP_result_plot_heatmap)

####     2.3.2.3 [Volcano plots of specific contrasts](#output_DEP_result_plot_volcano)

####     2.3.2.4 [Customed volcano plots of specific contrasts](#output_DEP_result_plot_custom_volcano)

####     2.3.2.5 [Statistical plot](#output_DEP_result_plot_statistical)

####     2.3.2.6 [Normal distribution - Fit](#output_DEP_result_plot_normal_distribution)

####     2.3.2.7 [PCA plot](#output_DEP_result_plot_pca)

####     2.3.2.8 [Uniform manifold approximation and projection (UMAP) plot](#output_DEP_result_plot_umap)

####     2.3.2.9 [t-SNE plot](#output_DEP_result_plot_tsne)

####     2.3.2.10 [Scaling and variance stabilization](#output_DEP_result_plot_meansd)

####     2.3.2.11 [Correlation matrix](#output_DEP_result_plot_pearson)

####     2.3.2.12 [Distance matrix](#output_DEP_result_plot_gower)

####     2.3.2.13 [Sample CVs Plots](#output_DEP_result_plot_cv)

####     2.3.2.14 [Barplot of protein numbers](#output_DEP_result_plot_protein_number)

####     2.3.2.15 [Sample coverage](#output_DEP_result_plot_sample_coverage)

####     2.3.2.16 [Normalization plot](#output_DEP_result_plot_normalization)

####     2.3.2.17 [Missing values- Quant](#output_DEP_result_plot_missing_value_quant)

####     2.3.2.18 [Missing values- Heatmap](#output_DEP_result_plot_missing_value_heatmap)

####     2.3.2.19 [Imputation](#output_DEP_result_plot_imputation)

#### 3. [Differential gene expression for count data of RNAseq (DEG-RNAseq)](#DEG_RNAseq)

####   3.1 [Overview](#overview_RNAseq)

####   3.2 [Analytical steps](#analytical_steps_RNAseq)

####   3.3 [Output](#output_RNAseq)

####    3.3.1 [Result tables](#output_RNAseq_result_tables)

####    3.3.2 [Visualization of the results](#output_RNAseq_result_plot)

####     3.3.2.1 [Dispersion plot](#output_RNAseq_result_plot_dispersion)

####     3.3.2.2 [PCA plot](#output_RNAseq_result_plot_pca)

####     3.3.2.3 [Uniform manifold approximation and projection (UMAP) plot](#output_RNAseq_result_plot_umap)

####     3.3.2.4 [t-SNE plot](#output_RNAseq_result_plot_tsne)

####     3.3.2.5 [Correlation matrix](#output_RNAseq_result_plot_pearson)

####     3.3.2.6 [Distance matrix](#output_RNAseq_result_plot_gower)

####     3.3.2.7 [Heatmap of all significant genes](#output_RNAseq_result_plot_heatmap)

####     3.3.2.8 [Customed volcano plots of specific contrasts](#output_RNAseq_result_plot_custom_volcano)

####     3.3.2.9 [MA plot](#output_RNAseq_result_plot_MA)

####     3.3.2.10 [Barplot of a gene of interest](#output_RNAseq_result_single_gene)

#### 4. [Gene list tool](#Genelist_tool)

####   4.1 [Overview](#overview_genelist_tool)

####   4.2 [Analytical steps](#analytical_steps_genelist_tool)

####   4.3 [Output](#output_genelist_tool)

####    4.3.1 [Result tables](#output_genelist_tool_result_tables)

####    4.3.2 [Visualization of the results](#output_genelist_tool_result_plot)

####     4.3.2.1 [Venn diagram](#output_genelist_tool_result_plot_venn)

#### 5. [Gene annotation (Annotation)](#annotation)

####   5.1 [Overview](#overview_annotation)

####   5.2 [Analytical steps](#analytical_steps_annotation)

####   5.3 [Output](#output_annotation)

#### 6. [Over Representation Analysis (ORA)](#ORA_header)

####   6.1 [Overview](#overview_ORA)

####   6.2 [Analytical steps](#analytical_steps_ORA)

####   6.3 [Output](#output_ORA)

####    6.3.1 [Result tables](#output_ORA_result_tables)

####    6.3.2 [Visualization of the results](#output_ORA_result_plot)

####     6.3.2.1 [Bar plot](#output_ORA_result_plot_bar)

####     6.3.2.2 [Dot plot](#output_ORA_result_plot_dot)

####     6.3.2.3 [Optimized dot plot](#output_ORA_result_plot_dot_opt)

####     6.3.2.4 [Heatmap-like functional classification](#output_ORA_result_plot_heat)

####     6.3.2.5 [Gene-Concept Network](#output_ORA_result_plot_cnet)

####     6.3.2.6 [Enrichment Map](#output_ORA_result_plot_emap)

####     6.3.2.7 [goplot](#output_ORA_result_plot_goplot)

####     6.3.2.8 [GO graph](#output_ORA_result_plot_gograph)

#### 7. [Gene Set Enrichment Analysis (GSEA)](#GSEA_header)

####   7.1 [Overview](#overview_GSEA)

####   7.2 [Analytical steps](#analytical_steps_GSEA)

####   7.3 [Output](#output_GSEA)

####    7.3.1 [Result tables](#output_GSEA_result_tables)

####    7.3.2 [Visualization of the results](#output_GSEA_result_plot)

####     7.3.2.1 [Bar plot](#output_GSEA_result_plot_bar)

####     7.3.2.2 [Dot plot](#output_GSEA_result_plot_dot)

####     7.3.2.3 [Heatmap-like functional classification](#output_GSEA_result_plot_heat)

####     7.3.2.4 [Gene-Concept Network](#output_GSEA_result_plot_cnet)

####     7.3.2.5 [Enrichment Map](#output_GSEA_result_plot_emap)

####     7.3.2.6 [Running score and preranked list of GSEA result](#output_GSEA_result_plot_gseaplot)

#### 8. [String protein-protein interaction network (PPITools)](#PPI)

####   8.1 [Overview](#overview_PPI)

####   8.2 [Analytical steps](#analytical_steps_PPI)

####   8.3 [Output](#output_PPI)

####    8.3.1 [Result table](#output_PPI_result_tables)

####    8.3.2 [Visualization of the results](#output_PPI_result_plot)

####     8.3.2.1 [Network](#output_PPI_result_plot_network)

#### 9. [Heatmap of proteomics and RNAseq data together with the same order (RR-Heatmap)](#PR_heatmap)

####   9.1 [Overview](#overview_PR_heatmap)

####   9.2 [Analytical steps](#analytical_steps_PR_heatmap)

####   9.3 [Output](#output_PR_heatmap)

####    9.3.1 [Result table](#output_PR_heatmap_result_tables)

####    9.3.2 [Visualization of the results](#output_PR_heatmap_result_plot)

####     9.3.2.1 [Heatmaps of proteomics and RNAseq data together with the same order](#output_PR_heatmap_result_plot_heatmap)

#### 10. [Some useful links for convenience (Extralinks)](#extralinks)

#### 11. [Updates news](#updates_news)

<hr>

<a name="Overview"></a> \#\#\# **1. Overview of the analysis**

<p style="font-size:1.0em;text-align:justify">

This shiny application provides a workflow for mass spectrometry proteomics data such as label-free or other quantification data processed by MaxQuant or other software such as IsobarQuant. The full tutorial is available online and embedded within the application for users to access. Default parameters are provided to guide unfamiliar users through the analytical steps. Generally, the order of the analytical steps are from the top to the bottom, and default parameters are proper, also if desired, could be changed based on user preference. This application aims to build a whole workflow to help users to explore data.It supports: </br> </br>   1.   Differential analysis of proteomics data based on R package <a href="https://www.bioconductor.org/packages/release/bioc/html/limma.html">Limma</a> and <a href="https://www.bioconductor.org/packages/release/bioc/html/DEP.html">DEP</a> </br>

  2.   Differential analysis of RNAseq countmatrix data based on R package <a href="https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html">DEseq2</a> </br>

  3.   Gene list tool supporting to generate, import and compare gene lists which also could be imported to some other parts of this application </br>

  4.   Gene annotation including gene description, entrez id, the name of <a href="http://geneontology.org/">gene ontology</a> BP(biological process), CC(cellular component), MF(molecular function), the name of <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome pathway</a>, and the name of <a href="http://pfam.xfam.org/">pfam</a> </br>

  5.   ORA(over-representation analysis) including <a href="http://geneontology.org/">gene ontology</a>, <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome</a> pathway enrichment based on the most widely used and popular R package <a href="https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html">clusterProfiler</a> and <a href="https://bioconductor.org/packages/release/bioc/html/ReactomePA.html">ReactomePA</a> </br>

  6.   GSEA(Gene Set Enrichment Analysis) including <a href="http://geneontology.org/">gene ontology</a>, <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome</a> pathway enrichment based on the most widely used and popular R package <a href="https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html">clusterProfiler</a> and <a href="https://bioconductor.org/packages/release/bioc/html/ReactomePA.html">ReactomePA</a> </br>

  7.   PPI(protein--protein interactions) based on <a href="https://string-db.org/">STRING</a> </br>

  8.   Compare of differential analysis of proteomics and countmatrix of RNAseq data with a heatmap together </br>

  9.   some extra links for convenience </br> </br> Finally, considering the flexibility, this appalication is designed: </br> </br>   1.   Every part could either integrate with the results of differential analysis of proteomics data and RNAseq countmatrix data, or run separately. It allows user to have more flexible options of analytical steps </br>

  2.   Every part has a series of parameters that could be modified to help experienced user. There will also appear the detailed explanation when you hover on a parameter to help understand </br>

  3.   The results of every part including tables and figures that help explore data could be downloaded. Each of the main figures has a series of parameters that could be modified to allow user to custom such as heatmap, volcano, pca, pearson correlation, Gowers's distance plot and so on. And every figure has width and height parameters to allow user to get proper and nice pdf file </br> </br> Summary, this application provides a whole workflow for differential analysis of mass spectrometry proteomics data such as label-free or other quantification data processed by MaxQuant or other software such as IsobarQuant, optional compared with differential analysis of RNAseq countmatrix data considering that the level of protein is not always consistent with the level of RNA with a heatmap together. And it provides gene list tool that supporting to generate, import and compare gene lists which also could be imported to some other parts of this application. It could also give user insight into gene annotation, <a href="http://geneontology.org/">gene ontology</a>, <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome</a> pathway enrichment including both ORA and GSEA, as well as protein--protein interactions to explore data such as biological mechanism. Flexibly, each part could either be an independent functional module or integrate with the results of differential analysis of proteomics data and RNAseq countmatrix data. Finally, each part has a series of parameters that could be modified and appear detailed explanations when you hover on them, and also has a series of results including tabels and figures to help explore data which could always be downloaded properly by modifying some parameters. So, no matter you want to do the differential analysis of proteomics data, or the differential analysis of RNAseq countmatrix data, or gene annotation, or ORA, or GSEA, or PPI separately, or integrate differential analysis of proteomics data and RNAseq countmatrix data with each other or the other parts, this application will be an ideal choice.

</p>

<hr>

<a name="DEP-LFQ"></a> \#\#\# **2. Differential protein expression for label-free quantification (DEP-LFQ) or other quantification** <a name="overview_DEP"></a> \#\#\#\# **2.1 Overview** R package <a href="https://www.bioconductor.org/packages/release/bioc/html/limma.html">Limma</a> is the best generalized approach for intensity data such as those used in proteomics. And R package <a href="https://www.bioconductor.org/packages/release/bioc/html/DEP.html">DEP</a> provides an analysis workflow of mass spectrometry proteomics data for differential protein expression based on protein-wise linear models and empirical Bayes statistics using limma. We use these two packages to build the flexible and interactive shiny framework. The main steps include data upload, filtering, variance normalization and imputation of missing values, as well as statistical testing of differentially enriched / expressed proteins. Although R package DEP contains a shiny application, our framework is more flexible, elegant, ideal and friendly to users who don't know how to use R language. It provides more options for statistical test such as the expression columns, the threshold for the allowed number of missing values in at least one condition, the imputation method, the pvalue adjustment method, the threshold method, the type of contrasts and so on. Meanwhile, It supports not only label-free quantification data, but also other quantification such as TMT(Tandem Mass Tag) data. It also provides more figures such as pca plot and each figure has a series of parameters such as color helping custom and explore the data better. Finally, each figure could be downloaded as a proper pdf file through parameters **width** and **height**.

<a name="analytical_steps_DEP"></a> \#\#\#\# **2.2 Analytical steps** `Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br>  1) Please go to the `DEP-LFQ options` on the sidebar </br>  2) Go to `Files`, and the file requirements and examples could be found when click on the question mark on the right </br>   2.1) Upload file **ProteinGroups.txt** processed by MaxQuant or other software such as IsobarQuant </br>   2.2) Optional, Upload file **ExperimentalDesign.txt** when the **Sample annotation** is set to **Use Experimental Design** </br>  3) Go to `Columns` </br>   3.1) If your ProteinGroups.txt file contains Gene.names column, **Name column** is set automatically. If it contains Protein.IDs column, **ID column** is set automatically. **Delimiter** separating the name and id column is set to ";" by default. If it contains LFQ intensity columns, the **Expression columns** are set to LFQ intensity columns automatically. If it contains Peptides column, **peptide column** is set automatically. If it contains Reverse and Potential.contaminant columns, **Filter on columns** is set automatically. Of course, You can also set them by yourself. </br> `Note:` ensure that at least one of your name column and id column is non-empty; the peptide column and filter on columns can be empty. </br>   3.2) If desired, **allowed max.miss.num at least one condition** is allowed to set </br>   3.3) `Contrasts` on the top-right of the body is allowed to set to control or all or manual. When it is set to control or manual, there will appear parameter **Control** or **Manual test** that allows to set on the sidebar </br>  4) Go to `imputation options`, and select an imputation type. The explanations of the imputation methods could be found when click on the **Detailed information** </br>  5) Go to `Type of FDR correction`, and select the method of pvalue adjustment. We strongly suggest that keep the settings as default. If you have the mathematical and statistical knowledge, you could see paper [A unified approach to false discovery rate estimation](http://www.biomedcentral.com/1471-2105/9/303/), and [fdrtool: a versatile R package for estimating local and tail area- based false discovery rates](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461) to attempt to change the method. Note that: when you set parameter **correct p values by** to p values, and **type of FDR value** to fdr by BH, it is the BH correction method </br>  6) Go to `Threshold method`, and select the method that the cutoff of significant proteins based on. **intersect** means by adjusted pvalue and log2 fold change cutoff and **curve** means by curvature and x0 fold cutoff. And if you want to know more about the method 'curve', please move to paper </br>  7) Click on `Analyze` button </br>  8) A series of parameters are allowed to set such as **adj.P value**, **Log2 fold change**, **curvature**, **x0 fold** and so on </br>  9) Optional but Strongly recommended, click on **save result RData** on the sidebar to save the main result object in order to load the RData and ProteinGroups.txt together next time for reproducibility

<a name="output_DEP"></a> \#\#\#\# **2.3 Output**

<a name="output_DEP_result_tables"></a> \#\#\#\# **2.3.1 Result tables** Users could download the result tabels on the sidebar

-   **results:** Result of all proteins. Includes Gene names, Protein Ids, p-values (each pairwise comparisons), Adjusted p-values (each pairwise comparisons), boolean values for significance (each pairwise comparisons, note that the significant column is TRUE if any pairwise comparisons is TRUE), Log2 fold changes/ ratios (each pairwise comparisons), centered protein intensity (log2 transformed) in each sample
-   **significant_proteins:** Result of significant proteins according to the adjusted pvalue and log2 fold change cutoff. The columns are the same as **results**
-   **displayed_subset:** Result of the displayed subset according to parameters **Select direct comparisons** and **Exclude direct comparisons** on the top of the result table. The columns are the same as **results**
-   **full_dataset:** Result of all proteins. The columns not only include columns of **results**, but also the columns of ProteinGroups.txt to allow user to see more information of proteins
-   **Result RData:** Users could download the result RData on the sidebar

Strongly recommended, click on **save result RData** on the sidebar to save the main result object in order to load the RData and ProteinGroups.txt together next time for reproducibility. Because, some imputation methods not always get the same result such as MinProb

<a name="output_DEP_result_plot"></a> \#\#\#\# **2.3.2 Visualization of the results** Each plot has a series of parameters such as color helping custom and explore the data better, and could be downloaded as a proper pdf file through parameters **width** and **height** </br> `Note:` </br> - The box containing parameters may be collapsed in order to save space. You can click on the plus sign (+) on the right of the box to expand the corresponding category - If not otherwise specified, in the following examples, the parameter **Threshold method** is set to intersect

<a name="output_DEP_result_plot_single_protein"></a> \#\#\#\# **2.3.2.1 Barplot of a protein of interest** In the shiny application, it is named as **Selected Protein**. It is a barplot of a protein that you are interested in. x-axis represents the conditions / groups, y-axis represents log2 centered intensity. Width and height could be modified in order to get a proper pdf file </br> <img src="DEP_selected_protein.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_heatmap"></a> \#\#\#\# **2.3.2.2 Heatmap of all significant proteins** In the shiny application, it is named as **Heatmap**. The heatmap representation gives an overview of all significant/differentially expressed proteins (rows) in all samples (columns). It means that the rows of the heatmap are those proteins significant/differentially expressed in at least one contrast. This visualization allows the identfication of general trends such as if one sample or replicate is highly different compared to the others and might be considerd as an outlier. Additionally, the hierarchical clustering of samples (columns) indicates how related the different samples are and hierarchicalclustering of proteins (rows) identifies similarly behaving proteins. The proteins can be clustered by k-means clustering (kmeans argument) and the number of clusters can be defined by argument **Kmeans**. There are a series of parameters that could be modified to allow to get a customed figure. Here, we will give some examples. In order to get a heatmap that row names are readable and do not take up space, we set the log2 fold change cut off very large to 5 </br> - **example 1:** The default parameters </br> <img src="DEP_heatmap_exampe1.png" style="margin-left:0px;" width="40%"/> </br> - **example 2:** **Color** parameter could be modified to change the color limit of the colorbar. Here the color is lighter than example 1 </br> <img src="DEP_heatmap_exampe2.png" style="margin-left:0px;" width="40%"/> </br> - **example 3:** If you want to order the clusters of the rows from top to bottom, you should check parameter **If mysplit** under the collapsed box, and set parameter **my split**. Also , you can set the order of the columns by unchecking **Cluster columns** and setting **Custom columns order**. Here, we set the order of clusters of the rows 1 to 6 from top to bottom, the order of the columns A to B from left to right </br> <img src="DEP_heatmap_exampe3.png" style="margin-left:0px;" width="40%"/> </br> - **example 4:** We provide parameter **row color** by checking it to show colored row names to check the number of peptides associated with each protein group(purple: Peptides == 1, blue: Peptides == 2). If no color of row names is purple or blue, it means that the number of peptides associated with all protein groups is lager than 2. So, user could check the accuracy of the quantification. The accuracy of the quantification of a protein having one or two peptides may be less than those having more than two peptides </br> <img src="DEP_heatmap_exampe4.png" style="margin-left:0px;" width="40%"/> </br> <a name="output_DEP_result_plot_heatmap_example_5"></a> - **example 5:** If your contrasts is more than one, and want to show specified contrasts, the parameters **Contrast** and **Manual heatmap** are provided. And it supports one or more specified contrasts. One specified contrast represents to show significant/differentially expressed proteins (rows) in this contrast. Two or more specified contrasts represents to show significant/differentially expressed proteins (rows) in at least one of these specified contrasts. What you should do is that, select your interested contrasts by setting parameter **Contrast**, and then check **Manual heatmap**. Here we choose contrast B_vs_A from B_vs_A and C_vs_A </br> <img src="DEP_heatmap_exampe5.png" style="margin-left:0px;" width="40%"/> </br> - **example 6:** Color palettes from R package RColorBrewer are allowed to specify by setting parameter **colorbar**. Here is the color palette RdGy </br> <img src="DEP_heatmap_exampe6.png" style="margin-left:0px;" width="40%"/> </br> - **example 7:** Alternatively, a heatmap can be plotted using the contrasts by setting **Data presentation** to contrast </br> <img src="DEP_heatmap_exampe7.png" style="margin-left:0px;" width="40%"/> </br> - **example 8:** You could also only display proteins which you are interested in by selecting proteins through parameter **selected proteins** first, and then check parameter **row selected** </br> <img src="DEP_heatmap_exampe8.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_volcano"></a> \#\#\#\# **2.3.2.3 Volcano plots of specific contrasts** In the shiny application, it is named as **Volcano plot**. Volcano plots can be used to visualize a specific contrast (comparison between two samples). It is a graphical visualization by plotting the "Fold Change (Log2)" on the x-axis versus the --log10 of the "p-value" on the y-axis. The contrast could be set under **Contrast**. Protein labels could be added by checking parameter **Display names**. **Font size** could be set for labels. And the labels should be omitted when there are too many labels. Y-axis could be adjusted p values by checking parameter **Adjusted p values**. The width from 0 on x axis could be the same by checking parameter **Same width**. And by selecting **Mybreaks** and checking **My breaks** under the collapsed box, the breaks of x axis could be customed. Here is an example </br> <img src="DEP_volcano_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_custom_volcano"></a> \#\#\#\# **2.3.2.4 Customed volcano plots of specific contrasts** In the shiny application, it is named as **Custom Volcano**. It is more flexible and suggested than **Volcano plot**. It has a series of parameters that can be set. For example, the contrast could be set under parameter **Contrast**. **Font size** and **dot size** could be set for labels. Y-axis could be adjusted p values by checking parameter **Adjusted p values**. The width from 0 on x axis could be the same by checking parameter **Same width**. Here, we will give some examples. `Note that:` the width from 0 on x axis is already the same in the following examples, so we don't check the parameter **Same width** </br> - **example 1:** The default parameters. By default, it shows the top 20 significant proteins </br> <img src="DEP_custom_volcano_example1.png" style="margin-left:0px;" width="40%"/> </br> - **example 2:** It could show the top 20 significant proteins of only up-regulated or down-regulated by setting **label way** to up or down. The example sets it to up </br> <img src="DEP_custom_volcano_example2.png" style="margin-left:0px;" width="40%"/> </br> - **example 3:** The black line around the labeled points can be removed by setting **point outside width** to 0 </br> <img src="DEP_custom_volcano_example3.png" style="margin-left:0px;" width="40%"/> </br> - **example 4:** It also could show specified genes that you are interested in. First, choose your interested genes under parameter **selected proteins** which is under the collapsed box, and then set **label way** to selected proteins. Note, to save time, the order of the steps should not be reversed </br> <img src="DEP_custom_volcano_example4.png" style="margin-left:0px;" width="40%"/> </br> - **example 5:** To be more readable, the labels could be surrounded by a rectangular by checking **label with rectangle** under the collapsed box </br> <img src="DEP_custom_volcano_example5.png" style="margin-left:0px;" width="40%"/> </br> - **example 6:** It could be no any labels or all labels of significant proteins by setting **show number** under the collapsed box to 0 or a lager number than the number of significant proteins. It is not recommended to show all labels when there are too many labels. Here is an example of no labels </br> <img src="DEP_custom_volcano_example6.png" style="margin-left:0px;" width="40%"/> </br> - **example 7:** Just like the heatmap plot, we provide parameter **peptide color** by checking it to show colored points to check the number of peptides associated with each protein group(purple: Peptides == 1, blue: Peptides == 2, you can also modify the color by yourself). If no color of points is purple or blue, it means that the number of peptides associated with all protein groups is lager than 2. So, user could check the accuracy of the quantification. The accuracy of the quantification of a protein having one or two peptides may be less than those having more than two peptides </br> <img src="DEP_custom_volcano_example7.png" style="margin-left:0px;" width="40%"/> </br> - **example 8:** The color of up-regulated, down-regulated and not significant proteins can be customed by parameters **up**, **down** and **not significant** </br> <img src="DEP_custom_volcano_example8.png" style="margin-left:0px;" width="40%"/> </br> - **example 9:** Here is an example that when parameter **Threshold method** is set to curve. And parameter **curvature** is set to 0.6 as default, parameter **x0 fold** is set to 1 as default. The larger the value of the two parameters, the stricter the filtering conditions </br> <img src="DEP_custom_volcano_example9.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_statistical"></a> \#\#\#\# **2.3.2.5 Statistical plot** In the shiny application, it is named as **Statistical plot**. Here we provides the plot based on -log10.padj, -log10.pval, padj, pval, and t statistic to help explore the data for users who have mathematical and statistical knowledge. It supports two **plot type**: x-y represents that X-axis represents your selected x, and y-axis represents your selected y; histogram represents the histogram plot of your selected x. This is an example of 'x-y' type, and the contrast is set to B_vs_A </br> <img src="DEP_statistical_plot_example1.png" style="margin-left:0px;" width="40%"/> </br> You can set the **plot type** to histogram </br> <img src="DEP_statistical_plot_example2.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_normal_distribution"></a> \#\#\#\# **2.3.2.6 Normal distribution - Fit** In the shiny application, it is named as **Normal distribution - Fit**. It is only for **threshold method** curve . Because this threshold method needs log2 fold change of contrast to be a normal distribution. So, we give the histogram plot with a fitted normal distribution curve on it to help users to check. Here the contrast is set to B_vs_A </br> <img src="DEP_norm_distribution_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_pca"></a> \#\#\#\# **2.3.2.7 PCA plot** In the shiny application, it is named as **Pca plot**. Principal component analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset. This can be very useful to observe batch effects, such as clear differences between replicates. The more similar 2 samples are, the closer they cluster together </br> <img src="DEP_pca_example1.png" style="margin-left:0px;" width="40%"/> </br> A square picture is allowed by checking **if square** </br> <img src="DEP_pca_exampe2.png" style="margin-left:0px;" width="40%"/> </br> Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition </br> <img src="DEP_pca_example3.png" style="margin-left:0px;" width="40%"/> </br> You could set the number of top variable proteins to consider by parameter **number of top variable proteins to consider**. The default is 500, and here we use all the proteins </br> <img src="DEP_pca_example4.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_umap"></a> \#\#\#\# **2.3.2.8 Uniform manifold approximation and projection (UMAP) plot** In the shiny application, it is named as **UMAP plot**. Uniform manifold approximation and projection ([UMAP](%22https://github.com/tkonopka/umap%22)) is a technique for dimensional reduction. It is similar to PCA. The more similar 2 samples are, the closer they cluster together </br> <img src="DEP_UMAP_example1.png" style="margin-left:0px;" width="40%"/> </br> A square picture is allowed by checking **if square** </br> <img src="DEP_UMAP_example2.png" style="margin-left:0px;" width="40%"/> </br> Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition </br> <img src="DEP_UMAP_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_tsne"></a> \#\#\#\# **2.3.2.9 t-SNE plot** In the shiny application, it is named as **t-SNE plot**. [t-SNE](%22chrome-extension://ikhdkkncnoglghljlkmcimlnlhkeamad/pdf-viewer/web/viewer.html?file=https%3A%2F%2Fcran.r-project.org%2Fweb%2Fpackages%2FRtsne%2FRtsne.pdf%22) is a method for constructing a low dimensional embedding of high-dimensional data, distances or similarities. It is also similar to PCA. The more similar 2 samples are, the closer they cluster together </br> <img src="DEP_TSNE_example1.png" style="margin-left:0px;" width="40%"/> </br> A square picture is allowed by checking **if square** </br> <img src="DEP_TSNE_example2.png" style="margin-left:0px;" width="40%"/> </br> Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition </br> <img src="DEP_TSNE_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_meansd"></a> \#\#\#\# **2.3.2.10 Scaling and variance stabilization** In the shiny application, it is named as **meanSdPlot**. The data is scaled and variance stabilized using [vsn](%22https://bioconductor.org/packages/3.14/bioc/html/vsn.html%22) </br> <img src="DEP_meansdplot_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_pearson"></a> \#\#\#\# **2.3.2.11 Correlation matrix** In the shiny application, it is named as **Pearson correlation**. A correlation matrix can be plotted as a heatmap to visualize the Pearson correlations between the different samples. We provide parameters **lower** and **upper**. Setting the lower correlation limit can compress the plottable data and give greater contrast to the plot. It is rare that two replicates will vary greatly in their correlation. Here we set **lower** to 0.7 to see a good contrast </br> <img src="DEP_pearson_example1.png" style="margin-left:0px;" width="40%"/> </br> You can display the correlation values on the plot by checking **Add values** under the collapsed box. The decimal point digits and font size of the values also can be set by **Value digits** and **Value size**. Here we set the decimal point digits of the values to 4 </br> <img src="DEP_pearson_example2.png" style="margin-left:0px;" width="40%"/> </br> The palette of the color could be changed by **color panel** which provides the palettes from R package RColorBrewer, and if you want to get a reversed colorbar, please check **pal rev**. Here, we set the palette to Blues </br> <img src="DEP_pearson_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_gower"></a> \#\#\#\# **2.3.2.12 Distance matrix** In the shiny application, it is named as **Gower's distance**. A distance matrix can be plotted as a heatmap to visualize the Gower's distances between the different samples </br> <img src="DEP_gower_example1.png" style="margin-left:0px;" width="40%"/> </br> You can display the distance values on the plot by checking **Add values** under the collapsed box. The decimal point digits and font size of the values also can be set by **Value digits** and **Value size**. Here we set the decimal point digits of the values to 2 </br> <img src="DEP_gower_example2.png" style="margin-left:0px;" width="40%"/> </br> The palette of the color could be changed by **color panel** which provides the palettes from R package RColorBrewer, and if you want to get a reversed colorbar, please check **pal rev**. Here, we set the palette to Reds </br> <img src="DEP_gower_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_cv"></a> \#\#\#\# **2.3.2.13 Sample CVs Plots** In the shiny application, it is named as **Sample CVs**. It is a plot representing distribution of protein level coefficient of variation for each condition. Each plot also contains a vertical line representing median CVs percentage within that condition </br> <img src="DEP_cvplot_example.png" style="margin-left:0px;" width="60%"/> </br>

<a name="output_DEP_result_plot_protein_number"></a> \#\#\#\# **2.3.2.14 Barplot of protein numbers** In the shiny application, it is named as **Protein Numbers**. It is a barplot representing number of proteins identified and quantified in each sample </br> <img src="DEP_protein_number_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_sample_coverage"></a> \#\#\#\# **2.3.2.15 Sample coverage** In the shiny application, it is named as **Sample coverage**. It is a plot highlighting overlap between identified proteins across all samples in the experiment </br> <img src="DEP_protein_coverage_example.png" style="margin-left:0px;" width="20%"/> </br>

<a name="output_DEP_result_plot_normalization"></a> \#\#\#\# **2.3.2.16 Normalization plot** In the shiny application, it is named as **Normalization**. The data is normalized by variance stabilizing transformation (vsn). The normalization can be inspected by checking the distributions of the samples before and after normalization </br> <img src="DEP_normalization_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_missing_value_quant"></a> \#\#\#\# **2.3.2.17 Missing values- Quant** In the shiny application, it is named as **Missing values- Quant**. the densities and cumulative fractions are plotted for proteins with and without missing values to check whether missing values are biased to lower intense proteins </br> <img src="DEP_missing_values_quant_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_missing_value_heatmap"></a> \#\#\#\# **2.3.2.18 Missing values- Heatmap** In the shiny application, it is named as **Missing values- Heatmap**. Only proteins with at least one missing value are visualized. It could assist users to identify patterns of missing values (i.e., MAR: Missing values at random or MNAR: Missing values not at random) and guide the choice of imputation method, if desired. If you want to order the columns by yourself, please select your ordered column names through parameter **Custom columns order**, and then uncheck **Cluster columns** </br> <img src="DEP_missing_value_heatmap_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_DEP_result_plot_imputation"></a> \#\#\#\# **2.3.2.19 Imputation** In the shiny application, it is named as **Imputation**. The effect of the imputation on the distributions can be visualized through a desity plot of protein intensity (log2) distrubution for each condition after and before missing value imputation being performed </br> <img src="DEP_imputation_example.png" style="margin-left:0px;" width="40%"/> </br>

<hr>

<a name="DEG_RNAseq"></a> \#\#\# **3. Differential gene expression for count data of RNAseq (DEG-RNAseq)** <a name="overview_RNAseq"></a> \#\#\#\# **3.1 Overview** A basic task in the analysis of count data from RNA-seq is the detection of differentially expressed genes. We use the most widely used and popular R package [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to perform differential gene expression analysis for count data of RNAseq. For more information about DEseq2, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). It needs a count matrix and an experimental design table as input. It provides tables to allow to download, a series of figures and parameters that allow user to change. Finally, each figure could be downloaded as a proper pdf file through parameters **width** and **height**.

<a name="analytical_steps_RNAseq"></a> \#\#\#\# **3.2 Analytical steps** `Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br>  1) Please go to the `DEG-RNAseq options` on the sidebar </br>  2) Go to `Files`, and the file requirements and examples could be found when click on the question mark on the right </br>   2.1) Upload file **Countmatrix.txt** of RNAseq data </br>   2.2) Upload file **ExperimentalDesign.txt** </br>  3) Go to `RNAseq_settings` </br>   3.1) Select your design, and it can be multi-factor by parameter **Select design**. This is normally a subset of the variables in the experimental design table, which constitute the main source of information. For more information, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start) </br>   3.2) Choose the contrast to build upon by parameter **Build the contrast upon** </br>    3.2.1) When the contrast of interest has three or more levels available, the likelihood ratio test can be used instead of the Wald test, to allow for an ANOVA-like analysis across groups. Parameters **perform a LRT test** and **Select the reduced model** will appear automatically. For more information, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test) </br>   3.3) `Contrasts` on the top-right of the body is allowed to set to control or all or manual. When it is set to control or manual, there will appear parameter **Control** or **Manual test** that allows to set </br>   3.4) Set the threshold on the row sums of the counts by parameter **Threshold row sums**. For more information, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering) </br>  4) Go to `ID  Transformation`, and choose whether to transform your id to gene symbol by checking or unchecking parameter **ID transformation**. It is not needed when your input id type is already gene symbol. And it is strongly recommended when your input id type is `ENSEMBL` id in order to make the result more readable </br>   4.1) When you check parameter **ID transformation**, it will appear parameter **No-mapping to rowname** automatically which allows you to set according to your own need. Then select the species by parameter **Select species**. Next **select your id type** </br> `Note:` if you don't have the R package of genome wide annotation for the selected species, it will give a message below parameter **Select species** on the sidebar, please download the corresponding R package following it first </br>  5) Go to `result options`, and check or uncheck **Independent filtering**, **Shrink lfc**, **Use IHW**

-   **Independent filtering:** whether independent filtering should be applied automatically. For more information, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt)
-   **Shrink lfc:** whether shrink the log fold change for the contrast of interest. For more information, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#log-fold-change-shrinkage-for-visualization-and-ranking)
-   **Use IHW:** Whether use Independent Hypothesis Weighting (IHW) as a filtering function. For more information, here is a [link](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#independent-hypothesis-weighting)

 6) Choose how many cores to use for computing by slider **Cores to use** depending on the available resources </br>  7) Click on `Analyze` button </br>  8) A series of parameters are allowed to set such as **adj.P value**, **Log2 fold change** and so on

<a name="output_RNAseq"></a> \#\#\#\# **3.3 Output** <a name="output_RNAseq_result_tables"></a> \#\#\#\# **3.3.1 Result tables** Users could download the result tabels on the sidebar - **full_results:** Result of all genes after filtering, it includes the result of differential gene expression analysis and the normalized count data - **significant_results:** Significant result of differential gene expression analysis according to the adjusted pvalue and Log2 fold change cutoff

<a name="output_RNAseq_result_plot"></a> \#\#\#\# **3.3.2 Visualization of the results** For every figure, you can change **width** and **height** in order to get a proper pdf file. <a name="output_RNAseq_result_plot_dispersion"></a> \#\#\#\# **3.3.2.1 Dispersion plot** In the shiny application, it is named as **Diagno dispests**. Plotting the dispersion estimates is a useful diagnostic. The dispersion plot below is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates </br> <img src="RNAseq_diagno_dispests_example.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_pca"></a> \#\#\#\# **3.3.2.2 PCA plot** In the shiny application, it is named as **Pca plot**. Principal component analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset. This can be very useful to observe batch effects, such as clear differences between replicates. The more similar 2 samples are, the closer they cluster together </br> <img src="RNAseq_pca_exampe1.png" style="margin-left:0px;" width="40%"/> </br>

A square picture is allowed by checking **if square** </br> <img src="RNAseq_pca_exampe2.png" style="margin-left:0px;" width="40%"/> </br>

Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition </br> <img src="RNAseq_pca_exampe3.png" style="margin-left:0px;" width="40%"/> </br>

You could set the number of top variable genes to consider by parameter **number of top variable genes to consider**. The default is 500, and here we use all the genes </br> <img src="RNAseq_pca_exampe4.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_umap"></a> \#\#\#\# **3.3.2.3 Uniform manifold approximation and projection (UMAP) plot** In the shiny application, it is named as **UMAP plot**. Uniform manifold approximation and projection ([UMAP](%22https://github.com/tkonopka/umap%22)) is a technique for dimensional reduction. It is similar to PCA. The more similar 2 samples are, the closer they cluster together </br> <img src="RNAseq_UMAP_example1.png" style="margin-left:0px;" width="40%"/> </br> A square picture is allowed by checking **if square** </br> <img src="RNAseq_UMAP_example2.png" style="margin-left:0px;" width="40%"/> </br> Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition </br> <img src="RNAseq_UMAP_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_tsne"></a> \#\#\#\# **3.3.2.4 t-SNE plot** In the shiny application, it is named as **t-SNE plot**. [t-SNE](%22chrome-extension://ikhdkkncnoglghljlkmcimlnlhkeamad/pdf-viewer/web/viewer.html?file=https%3A%2F%2Fcran.r-project.org%2Fweb%2Fpackages%2FRtsne%2FRtsne.pdf%22) is a method for constructing a low dimensional embedding of high-dimensional data, distances or similarities. It is also similar to PCA. The more similar 2 samples are, the closer they cluster together </br> <img src="RNAseq_TSNE_example1.png" style="margin-left:0px;" width="40%"/> </br> A square picture is allowed by checking **if square** </br> <img src="RNAseq_TSNE_example2.png" style="margin-left:0px;" width="40%"/> </br> Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition </br> <img src="RNAseq_TSNE_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_pearson"></a> \#\#\#\# **3.3.2.5 Correlation matrix** In the shiny application, it is named as **Pearson correlation**. A correlation matrix can be plotted as a heatmap to visualize the Pearson correlations between the different samples. We provide parameters **lower** and **upper**. Setting the lower correlation limit can compress the plottable data and give greater contrast to the plot. It is rare that two replicates will vary greatly in their correlation. Here we set **lower** to 0.97 to see a good contrast </br> <img src="RNAseq_pearson_example1.png" style="margin-left:0px;" width="40%"/> </br>

You can display the correlation values on the plot by checking **Add values** under the collapsed box. The decimal point digits and font size of the values also can be set by **Value digits** and **Value size**. Here we set the decimal point digits of the values to 4 </br> <img src="RNAseq_pearson_example2.png" style="margin-left:0px;" width="40%"/> </br>

The palette of the color could be changed by **color panel** which provides the palettes from R package RColorBrewer, and if you want to get a reversed colorbar, please check **pal rev**. Here, we set the palette to Blues </br> <img src="RNAseq_pearson_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_gower"></a> \#\#\#\# **3.3.2.6 Distance matrix** In the shiny application, it is named as **Gower's distance**. A distance matrix can be plotted as a heatmap to visualize the Gower's distances between the different samples\
</br> <img src="RNAseq_gower_example1.png" style="margin-left:0px;" width="40%"/> </br>

You can display the distance values on the plot by checking **Add values** under the collapsed box. The decimal point digits and font size of the values also can be set by **Value digits** and **Value size**. Here we set the decimal point digits of the values to 2 </br> <img src="RNAseq_gower_example2.png" style="margin-left:0px;" width="40%"/> </br>

The palette of the color could be changed by **color panel** which provides the palettes from R package RColorBrewer, and if you want to get a reversed colorbar, please check **pal rev**. Here, we set the palette to Reds </br> <img src="RNAseq_gower_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_heatmap"></a> \#\#\#\# **3.3.2.7 Heatmap of all significant genes** In the shiny application, it is named as **Heatmap**. The heatmap representation gives an overview of all significant/differentially expressed genes (rows) in all samples (columns). It means that the rows of the heatmap are those genes significant/differentially expressed in at least one contrast. This visualization allows the identfication of general trends such as if one sample or replicate is highly different compared to the others and might be considerd as an outlier. Additionally, the hierarchical clustering of samples (columns) indicates how related the different samples are and hierarchicalclustering of genes (rows) identifies similarly behaving genes. The genes can be clustered by k-means clustering (kmeans argument) and the number of clusters can be defined by argument **Kmeans**. There are a series of parameters that could be modified to allow to get a customed figure. Here, we will give some examples. In order to get a heatmap that row names are readable and do not take up space, we set the log2 fold change cut off very large to 4

-   **example 1:** The default parameters </br> <img src="RNAseq_heatmap_exampe1.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 2:** **Color** parameter could be modified to change the color limit of the colorbar. Here the color is lighter than example 1 </br> <img src="RNAseq_heatmap_exampe2.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 3:** If you want to order the clusters of the rows from top to bottom, you should check parameter **If mysplit** under the collapsed box, and set parameter **my split**. Also , you can set the order of the columns by unchecking **Cluster columns** and setting **Custom columns order**. Here, we set the order of clusters of the rows 2 to 1 from top to bottom, the order of the columns untrt to trt from left to right </br> <img src="RNAseq_heatmap_exampe3.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 4:** If your contrasts is more than one, and want to show specified contrasts, the parameters **Contrast** and **Manual heatmap** are provided. And it supports one or more specified contrasts. One specified contrast represents to show significant/differentially expressed genes (rows) in this contrast. Two or more specified contrasts represents to show significant/differentially expressed genes (rows) in at least one of these specified contrasts. What you should do is that, select your interested contrasts by setting parameter **Contrast**, and then check **Manual heatmap**. Because our example has only one contrast, you can see the [example](#output_DEP_result_plot_heatmap_example_5) in the `Differential protein expression for label-free quantification (DEP-LFQ) or other quantification` part

-   **example 5:** Color palettes from R package RColorBrewer are allowed to specify by setting parameter **colorbar**. Here is the color palette RdGy </br> <img src="RNAseq_heatmap_exampe5.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 6:** Alternatively, a heatmap can be plotted using the raw normalized count by setting **Data presentation** to raw </br> <img src="RNAseq_heatmap_exampe6.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 7:** You could also only display genes which you are interested in by selecting genes through parameter **selected genes** first, and then check parameter **row selected** </br> <img src="RNAseq_heatmap_exampe7.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_custom_volcano"></a> \#\#\#\# **3.3.2.8 Customed volcano plots of specific contrasts** In the shiny application, it is named as **Custom Volcano**. Volcano plots can be used to visualize a specific contrast (comparison between two samples). It is a graphical visualization by plotting the "Fold Change (Log2)" on the x-axis versus the --log10 of the "p-value" on the y-axis. It has a series of parameters that can be set. For example, the contrast could be set under parameter **Contrast**. **Font size** and **dot size** could be set for labels. Y-axis could be adjusted p values by checking parameter **Adjusted p values**. It is flexible. Here, we will give some examples

-   **example 1:** The default parameters. By default, it shows the top 20 significant genes </br> <img src="RNAseq_custom_volcano_example1.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 2:** The width from 0 on x axis could be the same by checking **Same width** </br> <img src="RNAseq_custom_volcano_example2.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 3:** It could show the top 20 significant genes of only up-regulated or down-regulated by setting **label way** to up or down. The example sets it to up </br> <img src="RNAseq_custom_volcano_example3.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 4:** The black line around the labeled points can be removed by setting **point outside width** to 0 </br> <img src="RNAseq_custom_volcano_example4.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 5:** It also could show specified genes that you are interested in. First, check **labeled** which is under the collapsed box. Then, parameter **selected genes** will appear, and choose your interested genes under it. Finally, set **label way** to selected genes. Note, to save time, the order of the steps should not be reversed </br> <img src="RNAseq_custom_volcano_example5.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 6:** To be more readable, the labels could be surrounded by a rectangular by checking **label with rectangle** under the collapsed box </br> <img src="RNAseq_custom_volcano_example6.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 7:** It could be no any labels or all labels of significant genes by setting **show number** under the collapsed box to 0 or a lager number than the number of significant genes. It is not recommended to show all labels when there are too many labels. Here is an example of no labels </br> <img src="RNAseq_custom_volcano_example7.png" style="margin-left:0px;" width="40%"/> </br>
-   **example 8:** The color of up-regulated, down-regulated and not significant genes can be customed by parameters **up**, **down** and **not significant** </br> <img src="RNAseq_custom_volcano_example8.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_plot_MA"></a> \#\#\#\# **3.3.2.9 MA plot** In the shiny application, it is named as **MA plot**. MA plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if it met the adjusted p value and log2 fold change cutoff. It has a series of parameters that can be set. For example, the contrast could be set under parameter **Contrast**. Here, we will give some examples

-   **example 1:** The default parameters </br> <img src="RNAseq_MA_plot_example1.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 2:** A rug plot is a compact visualisation designed to supplement a 2d display with the two 1d marginal distributions. And you can get a rug plot by checking parameter **Add rug** </br> <img src="RNAseq_MA_plot_example2.png" style="margin-left:0px;" width="40%"/> </br>

-   **example 3:** It also could show specified genes that you are interested in. First, check parameter **Selected genes**. Then, parameter **selected genes** will appear, and choose your interested genes under it. Finally, check parameter **Labeled**. Note, to save time, the order of the steps should not be reversed </br> <img src="RNAseq_MA_plot_example3.png" style="margin-left:0px;" width="40%"/> </br>

<a name="output_RNAseq_result_single_gene"></a> \#\#\#\# **3.3.2.10 Barplot of a gene of interest** In the shiny application, it is named as **Selected gene**. It is a barplot of a gene that you are interested in. x-axis represents the conditions / groups, y-axis represents normalized counts. You can choose a gene that you are interested in under parameter **selected genes**. Add labels for the points by checking **Add labels**. The parameter **Labels repel** represents that whether the text labels repel away from each other and away from the data points. Y-axis can start at 0 by checking parameter **ylimZero**. And you can select an interesting group by parameter **Color by**. Finally, width and height could be modified in order to get a proper pdf file </br> <img src="RNAseq_selected_gene.png" style="margin-left:0px;" width="40%"/> </br>

<hr>

<a name="Genelist_tool"></a> \#\#\# **4. Gene list tool** <a name="overview_genelist_tool"></a> \#\#\#\# **4.1 Overview** To help users explore the significant results of DEP-LFQ, DEG-RNAseq parts or the imported gene lists better, we provide the gene list tool. It supports to generate, import and compare gene lists. And the gene lists could be imported to some other parts of this application including annotation, ORA, and PPI parts. There are three panels:

-   **Generate gene list by cutoff:** the significant results of DEP-LFQ or DEG-RNAseq parts will appear automatically based on all contrasts under existing gene lists. Of course, users could generate the significant result of DEP-LFQ or DEG-RNAseq based on some selectable options such as contrast, adjusted pvalue cutoff, log2 fold change cutoff or significant type (up, down or up and down) and so on. And if you adjust the threshold parameters such as adjusted pvalue cutoff, log2 fold change cutoff of DEP-LFQ or DEG-RNAseq, the gene lists will update automatically
-   **Import gene list:** you can import gene list which you are interested in by pasting them in the paste box. And you can name it. The imported gene lists will be under existing imported gene lists
-   **Compare gene lists:** all the significant results of DEP-LFQ or DEG-RNAseq parts or imported gene lists will appear in the **List drag input** box. You can drag them to boxes A, B, C or D. And you can get the venn diagram and corresponding tables by clicking on **Compare** button

<a name="analytical_steps_genelist_tool"></a> \#\#\#\# **4.2 Analytical steps**  1) Please go to the `Gene list tool options` on the sidebar </br>  2) There are three panels: </br>   2.1) **Generate gene list by cutoff** </br> `Note that:` to use this part, ensure that you have completed at least one of DEP-LFQ and DEG-RNAseq parts </br>    2.1.1) Please set **adjusted P cutoff**, **log2 fold change cutoff** and select **Generate gene list of** DEP (for proteomics) or DEG (for RNAseq), comparison (contrast), and **Significant type** (up, down or up and down) </br> `Note that:` when the threshold method is set to curve of DEP-LFQ part, the parameters **adjusted P cutoff** and **log2 fold change cutoff** will automatically become **curvature cutoff** and **x0 fold cutoff** </br>    2.1.2) Click on the **Generate** button, and it will give a message of your generated gene list including the saved name and the significant number of it. The name of generated gene list will be under **existing gene lists** </br>   2.2) **Import gene list** </br>    2.2.1) Paste the genes that you are interested in into the paste box </br>    2.2.2) Set the name of your imported gene list </br>    2.2.3) Click on the **Import** button, and it will give a message of your imported gene list including the total number of it. At the same time, the list of imported genes will be displayed below the button and the name of it will be under **existing imported gene lists** </br>   2.3) **Compare gene lists** </br>    2.3.1) The **List drag input** box contains all the significant results of DEP-LFQ or DEG-RNAseq parts or imported gene lists. Drag the gene lists that you are interested in to boxes A, B, C or D </br> `Note that:` you should select at least two gene lists to compare, and if you don't want a gene list, you can drag it back to the list drag input box </br>    2.3.2) Click on the **Compare** button, the venn diagram and corresponding tables will appear. And you can download them

<a name="output_genelist_tool"></a> \#\#\#\# **4.3 Output** <a name="output_genelist_tool_result_tables"></a> \#\#\#\# **4.3.1 Result tables**

-   **Input panel:** the genes of your selected gene lists
-   **Output panel:** the tabel corresponding to the venn diagram including name, count, ratio, count_ratio, and text

<a name="output_genelist_tool_result_plot"></a> \#\#\#\# **4.3.2 Visualization of the results** <a name="output_genelist_tool_result_plot_venn"></a> \#\#\#\# **4.3.2.1 Venn diagram** - **example 1:** Here we select four gene lists as example </br> <img src="genelist_tool_example1.png" style="margin-left:0px;" width="40%"/> </br> - **example 2:** If you hover on an area of the venn diagram, the corresponding contents will appear </br> <img src="genelist_tool_example2.png" style="margin-left:0px;" width="40%"/> </br>

<hr>

<a name="annotation"></a> \#\#\# **5. Gene annotation (Annotation)** <a name="overview_annotation"></a> \#\#\#\# **5.1 Overview** It is usually necessary to do a gene annotation analysis when you get a list of differentially enriched / expressed genes or genes that you are interested in. So you could know more information about these genes so as to understand the relevant biological significance. Here we provide the annotations of gene description, entrez id, the name of [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), the name of [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway, and the name of [pfam](http://pfam.xfam.org/). SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported as input. We support not only pasting a gene list directly, but also importing from gene list tool to make this application flexible. Finally, we provide the result table of gene annotation to allow user to download. Note that, genes that have no entrez ids will not appear in the result table.

<a name="analytical_steps_annotation"></a> \#\#\#\# **5.2 Analytical steps**  1) Please go to the `Annotation options` on the sidebar </br>  2) You can choose to paste a gene list directly or import from gene list tool </br>   2.1) If you have a gene list that you are interested in, please paste them directly in the box under **Please paste your gene list** </br>   2.2) If you want to import from gene list tool, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part or you have imported some gene lists from gene list tool, and then check **Import from gene list tool**. It will appear parameter **Choose gene list** which provides gene lists obtained from gene list tool for you to choose. And select one </br>  3) Select the organism </br>  4) Click on the **Analyze** button </br> `Note:` if you don't have the R package of genome wide annotation for the selected organism, it will give a message on the body of this shiny application, please download the corresponding R package following it first </br>  5) A table of the result of gene annotation will appear and allow user to download through the **Save table** button on the sidebar

<a name="output_annotation"></a> \#\#\#\# **5.3 Output** Users could download the result tabel on the sidebar - **Result table:** Result of gene annotation. Includes gene name, gene description, entrez id, the name of [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), the name of [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway, and the name of [pfam](http://pfam.xfam.org/)

<hr>

<a name="ORA_header"></a> \#\#\# **6. Over Representation Analysis (ORA)** <a name="overview_ORA"></a> \#\#\#\# **6.1 Overview** Over Representation Analysis (ORA) is a widely used approach to determine whether known biological functions or processes are over-represented (= enriched) in an experimentally-derived gene list, e.g. a list of differentially expressed genes (DEGs). For further information, here is a [link](http://yulab-smu.top/clusterProfiler-book/chapter2.html#over-representation-analysis). We use the most widely used and popular R package [clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) to perform over representation analysis. It needs a gene list (eg: DEGs) as input, and SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported. It supports [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway enrichment. It provides tables to allow to download, a series of figures and parameters that allow user to change. And we support not only pasting a gene list directly, but also importing from gene list tool to make this application flexible. Except that gene ontology enrichment has a few more parameters such as **Ontology** and **removed redundancy of enriched GO terms** and figures output **Goplot** and **Gograph**, their design and analytical steps are exactly the same. So, in this tutorial, we will just take gene ontology enrichment as an example. Finally, each figure could be downloaded as a proper pdf file through parameters **width** and **height**.

<a name="analytical_steps_ORA"></a> \#\#\#\# **6.2 Analytical steps** `Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br>  1) Please go to the `ORA options` on the sidebar </br>  2) You can choose to paste a gene list directly or import from gene list tool </br>   2.1) If you have a gene list that you are interested in, please paste them directly in the box under **Please paste your gene list** </br>   2.2) If you want to import from gene list tool, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part or you have imported some gene lists from gene list tool, and then check **Import from gene list tool**. It will appear parameter **Choose gene list** which provides gene lists obtained from gene list tool for you to choose. And select one </br>  3) Select database GO, KEGG, or Reactome. You can select one or more </br>  4) Select the organism </br>  5) If you check **Import from gene list tool**. It will appear parameter **If with log2 fold change**. Check or uncheck it, it means that whether your gene list contains log2 fold change (total two columns), and it is used only in **Cnetplot** and **Emaplot** to fill the color. If you paste a gene list directly, it will automatically recognize whether the gene list that you paste has log2 fold change </br>  6) Click on the **Analyze** button and go to the panel corresponding to the database that you selected on the body of this shiny application </br> `Note:` if you don't have the R package of genome wide annotation for the selected organism, it will give a message on the body of this shiny application, please download the corresponding R package following it first </br>  7) A series of parameters are allowed to set such as **Ontology**, **adj. P value** and so on

<a name="output_ORA"></a> \#\#\#\# **6.3 Output** <a name="output_ORA_result_tables"></a> \#\#\#\# **6.3.1 Result tables** Users could download the result tabels at the top left of the table on the body of this shiny application

-   **full_results:** Result of all terms
-   **significant_results:** Result of significant terms according to the pvalue, adjusted pvalue and qvalue cutoff. Generally speaking, the default parameters are appropriate

<a name="output_ORA_result_plot"></a> \#\#\#\# **6.3.2 Visualization of the results** For every figure, you can change **width** and **height** in order to get a proper pdf file. <a name="output_ORA_result_plot_bar"></a> \#\#\#\# **6.3.2.1 Bar plot** In the shiny application, it is named as **Bar plot**. Bar plot is the most widely used method to visualize enriched terms. Here, it depicts the enrichment scores (p values or adjusted p values) and gene count as bar color and height. You can change the bar color by parameter **colorBy**. And you can change the number of top significant terms to show by setting **ShowCategory** </br> <img src="ORA_barplot_example1.png" style="margin-left:0px;" width="60%"/> </br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together </br> <img src="ORA_barplot_example2.png" style="margin-left:0px;" width="60%"/> </br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each BP, CC and MF </br> <img src="ORA_barplot_example3.png" style="margin-left:0px;" width="60%"/> </br>

<a name="output_ORA_result_plot_dot"></a> \#\#\#\# **6.3.2.2 Dot plot** In the shiny application, it is named as **Dot plot**. Dot plot is similar to bar plot with the capability to encode another score as dot size. Here, the x-axis represents gene ratio, the color represents p values or adjusted p values by setting **colorBy**, and the dot size represents gene count. And you can change the number of top significant terms to show by setting **ShowCategory** </br> <img src="ORA_dotplot_example1.png" style="margin-left:0px;" width="60%"/> </br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together </br> <img src="ORA_dotplot_example2.png" style="margin-left:0px;" width="60%"/> </br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each BP, CC and MF </br> <img src="ORA_dotplot_example3.png" style="margin-left:0px;" width="60%"/> </br>

<a name="output_ORA_result_plot_dot_opt"></a> \#\#\#\# **6.3.2.3 Optimized dot plot** In the shiny application, it is named as **Dot opt**. It is simplified relative to the dot plot. The color represents p values or adjusted p values by setting **colorBy**, and the dot size represents gene count. And you can change the number of top significant terms to show by setting **ShowCategory** </br> <img src="ORA_dot_opt_example.png" style="margin-left:0px;" width="80%"/> </br>

<a name="output_ORA_result_plot_heat"></a> \#\#\#\# **6.3.2.4 Heatmap-like functional classification** In the shiny application, it is named as **Heatplot**. It is similar to cnetplot, while displaying the relationships as a heatmap. The gene-concept network may become too complicated if user want to show a large number significant terms. The heatplot can simplify the result and more easy to identify expression patterns. The rows represent enriched items, the columns represent genes, and the color represents the log2 fold change when parameter **If with log2 fold change** is checked. From this plot, you can see genes of your gene list that the significant items include. And you can change the number of top significant terms to show by setting **ShowCategory** </br> This is an example that parameter **If with log2 fold change** is unchecked </br> </br> </br> <img src="ORA_heatplot_example1.png" style="margin-left:0px;" width="80%"/> </br> </br> </br>

This is an example that parameter **If with log2 fold change** is checked </br> </br> </br> <img src="ORA_heatplot_example2.png" style="margin-left:0px;" width="80%"/> </br> </br> </br>

<a name="output_ORA_result_plot_cnet"></a> \#\#\#\# **6.3.2.5 Gene-Concept Network** In the shiny application, it is named as **Cnetplot**. Both the barplot and dotplot only displayed most significant enriched terms, while users may want to know which genes are involved in these significant terms. The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network. The yellow point represents items, the colors of the gene points are filled with log2 fold change when parameter **If with log2 fold change** is checked. If a gene belongs to a item, then there will be a line between them. You can change the number of top significant terms to show by setting **ShowCategory** </br> This is an example that parameter **If with log2 fold change** is unchecked </br> </br> </br> <img src="ORA_cnetplot_example1.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

This is an example that parameter **If with log2 fold change** is checked </br> </br> </br> <img src="ORA_cnetplot_example2.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

This is an example that parameter **Circular** is unchecked to change the layout </br> </br> </br> <img src="ORA_cnetplot_example3.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

<a name="output_ORA_result_plot_emap"></a> \#\#\#\# **6.3.2.6 Enrichment Map** In the shiny application, it is named as **Emaplot**. Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module. You can change the number of top significant terms to show by setting **ShowCategory** </br> </br> </br> <img src="ORA_emaplot_example.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

<a name="output_ORA_result_plot_goplot"></a> \#\#\#\# **6.3.2.7 goplot** In the shiny application, it is named as **Goplot**. It can accept output of enrichGO and visualized the enriched GO induced graph. It shows the hierarchical relationships of terms </br> </br> </br> <img src="ORA_goplot_example.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

<a name="output_ORA_result_plot_gograph"></a> \#\#\#\# **6.3.2.8 GO graph** In the shiny application, it is named as **GOgraph**. It is another form of showing hierarchical relationships of terms. `Notice that`, it has its own action button, please click on the **plot** button when you want to plot it </br> </br> </br> <img src="ORA_gograph_example.png" style="margin-left:0px;" width="80%"/>

<hr>

<a name="GSEA_header"></a> \#\#\# **7. Gene Set Enrichment Analysis (GSEA)** <a name="overview_GSEA"></a> \#\#\#\# **7.1 Overview** Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states(e.g. phenotypes). For further information, here are some links, [link1](http://www.gsea-msigdb.org/gsea/index.jsp), [link2](http://yulab-smu.top/clusterProfiler-book/chapter2.html#gene-set-enrichment-analysis). And the difference between GSEA and ORA is that:

-   ORA mainly focuses on a few genes such as genes that are significantly up-regulated or down-regulated. So it is easy to miss some genes that are not significantly differentially expressed but have important biological significance.
-   GSEA is not limited to significantly differentially expressed genes, all genes can be used in GSEA. From the perspective of gene set enrichment, it is theoretically easier to include the impact of slight but coordinated changes on biological pathways, especially when gene sets with not too large differential multiples.

We use the most widely used and popular R package [clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) to perform gene set enrichment analysis. It needs a gene list containing log2 fold change (total two columns) as input, and SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported. It supports [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway enrichment. It provides tables to allow to download, a series of figures and parameters that allow user to change. And we support not only pasting a gene list containing log2 fold change (total two columns) directly, but also importing the result of DEP-LFQ or DEG-RNAseq part to make this application flexible. Except that gene ontology enrichment has a few more parameters such as **Ontology** and **removed redundancy of enriched GO terms**, their design and analytical steps are exactly the same. So, in this tutorial, we will just take gene ontology enrichment as an example. Finally, each figure could be downloaded as a proper pdf file through parameters **width** and **height**.

<a name="analytical_steps_GSEA"></a> \#\#\#\# **7.2 Analytical steps** `Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br>  1) Please go to the `GSEA options` on the sidebar </br>  2) You can choose to paste a gene list containing log2 fold change (total two columns) directly or import the result of DEP-LFQ or DEG-RNAseq part </br>   2.1) If you have a gene list containing log2 fold change (total two columns) that you are interested in, please paste them directly in the box under **Please paste your gene list** </br>   2.2) If you want to import the result of DEP-LFQ or DEG-RNAseq part, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part, and then check **Import from DEP-LFQ or DEG-RNAseq**. It will appear parameter **Choose gene list** which provides gene lists for you to choose. And select one </br>  3) Select database gseGO, gseKEGG, or gseReactome. You can select one or more </br>  4) Select the organism </br>  5) Click on the **Analyze** button and go to the panel corresponding to the database that you selected on the body of this shiny application </br> `Note:` if you don't have the R package of genome wide annotation for the selected organism, it will give a message on the body of this shiny application, please download the corresponding R package following it first </br>  6) A series of parameters are allowed to set such as **Ontology**, **adj. P value**, **NES**, **Phenotype** and so on

<a name="output_GSEA"></a> \#\#\#\# **7.3 Output** <a name="output_GSEA_result_tables"></a> \#\#\#\# **7.3.1 Result tables** Users could download the result tabels at the top left of the table on the body of this shiny application

-   **full_results:** Result of all terms
-   **significant_results:** Result of significant terms according to the pvalue, adjusted pvalue and NES (the absolute value of normalized enrichment score) cutoff. Generally speaking, the default parameters are appropriate

<a name="output_GSEA_result_plot"></a> \#\#\#\# **7.3.2 Visualization of the results** The parameter **Phenotype** means that the phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed. Default is both activated and suppressed </br> For every figure, you can change **width** and **height** in order to get a proper pdf file. <a name="output_GSEA_result_plot_bar"></a> \#\#\#\# **7.3.2.1 Bar plot** In the shiny application, it is named as **Bar plot**. Bar plot is the most widely used method to visualize enriched terms. Here, it depicts the enrichment scores (p values or adjusted p values) and NES as bar color and height. You can change the bar color by parameter **colorBy**. And you can change the number of top significant terms to show by setting **ShowCategory** </br> Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed </br> <img src="GSEA_barplot_example1.png" style="margin-left:0px;" width="50%"/> </br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to activated </br> <img src="GSEA_barplot_example2.png" style="margin-left:0px;" width="50%"/> </br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to suppressed </br> <img src="GSEA_barplot_example3.png" style="margin-left:0px;" width="50%"/> </br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed </br> <img src="GSEA_barplot_example4.png" style="margin-left:0px;" width="50%"/> </br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each phenotype that you selected based on each BP, CC and MF. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed </br> <img src="GSEA_barplot_example5.png" style="margin-left:0px;" width="60%"/> </br>

<a name="output_GSEA_result_plot_dot"></a> \#\#\#\# **7.3.2.2 Dot plot** In the shiny application, it is named as **Dot plot**. Dot plot is similar to bar plot with the capability to encode another score as dot size. Here, the x-axis represents NES, the color represents p values or adjusted p values by setting **colorBy**, and the dot size represents gene count. And you can change the number of top significant terms to show by setting **ShowCategory** </br> Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed </br> <img src="GSEA_dotplot_example1.png" style="margin-left:0px;" width="50%"/> </br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to activated </br> <img src="GSEA_dotplot_example2.png" style="margin-left:0px;" width="50%"/> </br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to suppressed </br> <img src="GSEA_dotplot_example3.png" style="margin-left:0px;" width="50%"/> </br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed </br> <img src="GSEA_dotplot_example4.png" style="margin-left:0px;" width="50%"/> </br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each phenotype that you selected based on each BP, CC and MF. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed </br> <img src="GSEA_dotplot_example5.png" style="margin-left:0px;" width="60%"/> </br>

<a name="output_GSEA_result_plot_heat"></a> \#\#\#\# **7.3.2.3 Heatmap-like functional classification** In the shiny application, it is named as **Heatplot**. It is similar to cnetplot, while displaying the relationships as a heatmap. The gene-concept network may become too complicated if user want to show a large number significant terms. The heatplot can simplify the result and more easy to identify expression patterns. The rows represent enriched items, the columns represent genes, and the color represents the log2 fold change. From this plot, you can see genes of your gene list that the significant items include. And you can change the number of top significant terms to show by setting **ShowCategory** </br> </br> </br> <img src="GSEA_heatplot_example.png" style="margin-left:0px;" width="100%"/> </br> </br> </br>

<a name="output_GSEA_result_plot_cnet"></a> \#\#\#\# **7.3.2.4 Gene-Concept Network** In the shiny application, it is named as **Cnetplot**. Both the barplot and dotplot only displayed most significant enriched terms, while users may want to know which genes are involved in these significant terms. The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network. The yellow point represents items, the colors of the gene points are filled with log2 fold change. If a gene belongs to a item, then there will be a line between them. You can change the number of top significant terms to show by setting **ShowCategory** </br> </br> </br> <img src="GSEA_cnetplot_example1.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

This is an example that parameter **Circular** is unchecked to change the layout </br> </br> </br> <img src="GSEA_cnetplot_example2.png" style="margin-left:0px;" width="60%"/> </br> </br> </br>

<a name="output_GSEA_result_plot_emap"></a> \#\#\#\# **7.3.2.5 Enrichment Map** In the shiny application, it is named as **Emaplot**. Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module. You can change the number of top significant terms to show by setting **ShowCategory** </br> </br> </br> <img src="GSEA_emaplot_example.png" style="margin-left:0px;" width="50%"/> </br> </br> </br>

<a name="output_GSEA_result_plot_gseaplot"></a> \#\#\#\# **7.3.2.6 Running score and preranked list of GSEA result** In the shiny application, it is named as **Gseaplot**. Running score and preranked list are traditional methods for visualizing GSEA result. It is a graphical visualization of the distribution of the gene set and the enrichment score. You can select the term that you want to show by parameter **Term**. The pvalue, adjusted pvalue and NES are displayed as a table on the plot </br> </br> </br> <img src="GSEA_gseaplot_example.png" style="margin-left:0px;" width="50%"/> </br> </br> </br>

<hr>

<a name="PPI"></a> \#\#\# **8. String protein-protein interaction network (PPITools)** <a name="overview_PPI"></a> \#\#\#\# **8.1 Overview** Proteins and their functional interactions form the backbone of the cellular machinery. Their connectivity network needs to be considered for the full understanding of biological phenomena. [String](https://string-db.org/) is a curated database with known and predicted protein-protein interactions and detailed description of their software can be found in their latest [manuscript](https://pubmed.ncbi.nlm.nih.gov/30476243/). We perform protein--protein interactions based on STRING. It needs a gene list as input, and SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported. It provides table and network to allow to download, and a series of parameters that allow user to change. We support not only pasting a gene list directly, but also importing from gene list tool to make this application flexible.

<a name="analytical_steps_PPI"></a> \#\#\#\# **8.2 Analytical steps** `Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br>  1) Please go to the `PPITools options` on the sidebar </br>  2) You can choose to paste a gene list directly or import from gene list tool </br>   2.1) If you have a gene list that you are interested in, please paste them directly in the box under **Please paste your gene list** </br>   2.2) If you want to import the result of DEP-LFQ or DEG-RNAseq part, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part or you have imported some gene lists from gene list tool, and then check **Import from gene list tool**. It will appear parameter **Choose gene list** which provides gene lists obtained from gene list tool for you to choose. And select one </br>  3) Select the organism </br>  4) Select concerned scores, it means that select which type of evidence will contribute to the prediction of the score under the active interaction sources </br>  5) Set the minimum required interaction score by parameter **scores cutoff**. The minimum required interaction score puts a threshold on the confidence score, such that only interaction above this score are included in the predicted network. Confidence limits are as follows:

-   low confidence - 0.15 (or better)
-   medium confidence - 0.4
-   high confidence - 0.7
-   highest confidence - 0.9

 6) Click on the **String it!** button </br> `Note:` if you don't have the R package of genome wide annotation for the selected organism, it will give a message on the body of this shiny application, please download the corresponding R package following it first. And if you don't have the string data of the selected organism, it will give some messages on the body of this shiny application, please download files from string following them first </br>  7) A series of parameters are allowed to set such as **input node color**, **input line color**, **choose node shape**, **choose layout style of network** and so on

<a name="output_PPI"></a> \#\#\#\# **8.3 Output** <a name="output_PPI_result_tables"></a> \#\#\#\# **8.3.1 Result table** Users could download the result tabel on the sidebar by clicking on **Save table** button which is named as **stringResult.txt** automatically - **stringResult.txt:** Data for the interaction network

<a name="output_PPI_result_plot"></a> \#\#\#\# **8.3.2 Visualization of the results** <a name="output_PPI_result_plot_network"></a> \#\#\#\# **8.3.2.1 Network** Users could download the network on the sidebar by clicking on **Save network** button which is named as **stringNetwork.html** automatically. And there are a series of parameters that could be modified to allow to get a customed figure such as **input node color(like:\#2EA9DF)**, **input line color(like: \#ADD8E6)**, **fonts size**, **choose layout style of network** and so on. For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br> </br> </br> This is an example </br> <img src="stringNetwork_example1.png" style="margin-left:0px;" width="30%"/> </br> This is an another example </br> <img src="stringNetwork_example2.png" style="margin-left:0px;" width="30%"/> </br>

<hr>

<a name="PR_heatmap"></a> \#\#\# **9. Heatmap of proteomics and RNAseq data together with the same order (RR-Heatmap)** <a name="overview_PR_heatmap"></a> \#\#\#\# **9.1 Overview** Mass spectrometry (MS)-based proteomics is powerful method for the quantitative profiling of proteins and enables increasingly comprehensive insights into changes of the proteome. RNASeq (named as an abbreviation of "RNA sequencing") is a sequencing technique which uses next-generation sequencing (NGS) to reveal the presence and quantity of RNA in a biological sample. Both of them are methods commonly used by biologists to explore biological mechanisms. Acknowledgedly, the expression level of a gene at the protein and mRNA level is not always consistent. Here we provide a heatmap of proteomics and RNAseq data together with the same order in order to help explore the expression level of a gene at the protein and mRNA level. It needs matrices of proteome and RNAseq normalized count, and a gene list that you are interested in as input. We support not only uploading files directly, but also importing the results of DEP-LFQ and DEG-RNAseq parts to make this application flexible. We provide table and heatmap to allow to download and a series of parameters that allow user to change. Finally, the heatmap could be downloaded as a proper pdf file through parameters **width** and **height**.

<a name="analytical_steps_PR_heatmap"></a> \#\#\#\# **9.2 Analytical steps** `Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter </br>  1) Please go to the `PR-Heatmap options` on the sidebar </br>  2) You can choose to upload files directly or import the results of DEP-LFQ and DEG-RNAseq parts </br>   2.1) If you have matrices of proteome and RNAseq normalized count, and a gene list that you are interested in, please upload them directly. The file requirements and examples could be found when click on the question mark on the right of **Proteinmatrix.txt** </br>   2.2) If you want to import the results of DEP-LFQ and DEG-RNAseq parts, please ensure that you have completed the analysises of DEP-LFQ and DEG-RNAseq parts, and ensure that the contrasts of these two analysises are the same. Then click on the menuitem **Import from** and check **DEP-LFQ and DEG-RNAseq**. It will appear the parameter **Contrast**, and then choose the contrast that you are interested in </br>  3) Click on the **Start** button </br>  4) A series of parameters are allowed to set such as **color for protein**, **color for rna** and so on

<a name="output_PR_heatmap"></a> \#\#\#\# **9.3 Output** <a name="output_PR_heatmap_result_tables"></a> \#\#\#\# **9.3.1 Result table** Users could download the result tabel on the sidebar - **Result table:** The table that the heatmap plot is based on

<a name="output_PR_heatmap_result_plot"></a> \#\#\#\# **9.3.2 Visualization of the results** <a name="output_PR_heatmap_result_plot_heatmap"></a> \#\#\#\# **9.3.2.1 Heatmaps of proteomics and RNAseq data together with the same order** In the shiny application, it is named as **Heatmap**. It is a heatmap of proteomics and RNAseq data together with the same order in order to help explore the expression level of a gene at the protein and mRNA level. It has a series of parameters that can be set. For example, **row font size**, **col font size**, **title size**, **Cluster columns**, **color for protein**, and **color for rna** could be set. And you can change **width** and **height** in order to get a proper pdf file. Here, we will give some examples

**`For uploading files directly:`**

**Note that:** if a row name shows grey color, it means that there is no corresponding protein or rnaseq data of this gene </br> </br> The parameter **Choose the genes to show:** </br> - **all:** all genes that you uploaded </br> - **both:** genes that both have RNAseq data and protein data of your uploaded </br> - **at least one:** genes that at least have one of RNAseq and protein data of your uploaded </br> </br> **example 1:** set parameter **Choose the genes to show** to all </br> <img src="PR_heatmap_exampe1.png" style="margin-left:0px;" width="40%"/> </br> **example 2:** set parameter **Choose the genes to show** to both </br> <img src="PR_heatmap_exampe2.png" style="margin-left:0px;" width="40%"/> </br> **example 3:** set parameter **Choose the genes to show** to at least one </br> <img src="PR_heatmap_exampe3.png" style="margin-left:0px;" width="40%"/> </br> **example 4:** Color palettes from R package RColorBrewer are allowed to specify by setting parameter **colorbar**. Here is the color palette RdGy </br> <img src="PR_heatmap_exampe4.png" style="margin-left:0px;" width="40%"/> </br> **example 5:** You could also only display genes which you are interested in by selecting genes through parameter **selected genes** first, and then check parameter **row selected** which are under the plus sign (+) on the top right of this shiny application body </br> <img src="PR_heatmap_exampe5.png" style="margin-left:0px;" width="40%"/> </br>

**`For importing the results of DEP-LFQ and DEG-RNAseq parts:`**

**Note that:** </br> (1) it is based on the union of your selected RNAseq and protein significantly differential genes </br> (2) if a row name shows grey color, it means that there is no corresponding protein or rnaseq data of this gene </br> (3) if a row name shows red color, it means that this gene is significant according to the analysis of DEP-LFQ or DEG-RNAseq part </br> </br> The parameter **Choose the genes to show:** </br> - **at least one:** genes that at least have one of RNAseq and protein data </br> - **both:** genes that both have RNAseq data and protein data </br> - **both significant:** both RNAseq data and protein data are significant </br> </br> **example 1:** set parameter **Choose the genes to show** to at least one </br> <img src="PR_heatmap_exampe_1.png" style="margin-left:0px;" width="40%"/> </br> **example 2:** set parameter **Choose the genes to show** to both </br> <img src="PR_heatmap_exampe_2.png" style="margin-left:0px;" width="40%"/> </br> **example 3:** set parameter **Choose the genes to show** to both significant </br> <img src="PR_heatmap_exampe_3.png" style="margin-left:0px;" width="40%"/> </br> **example 4:** Color palettes from R package RColorBrewer are allowed to specify by setting parameter **colorbar**. Here is the color palette RdGy </br> <img src="PR_heatmap_exampe_4.png" style="margin-left:0px;" width="40%"/> </br> **example 5:** You could also only display genes which you are interested in by selecting genes through parameter **selected genes** first, and then check parameter **row selected** which are under the plus sign (+) on the top right of this shiny application body </br> <img src="PR_heatmap_exampe_5.png" style="margin-left:0px;" width="40%"/> </br>

<hr>

<a name="extralinks"></a> \#\#\# **10. Some useful links for convenience (Extralinks)** Please go to the `Extralinks` on the sidebar. For convenience, here we have listed some links commonly used by biologists to explore biological issues. It includes [GeneCards](https://www.genecards.org/), [UniProt](https://www.uniprot.org/), [PubMed](https://pubmed.ncbi.nlm.nih.gov/), [HPA](https://www.proteinatlas.org/), [UALCAN](http://ualcan.path.uab.edu/index.html), [GEPIA](http://gepia.cancer-pku.cn/), and [cBioPortal](http://www.cbioportal.org/). For a brief description of each link, please see the `Extralinks` part of this application

<hr>

<a name="updates_news"></a> \#\#\# **11. Updates news** Please go to the `Updates news` on the sidebar. We will list the updates news here
