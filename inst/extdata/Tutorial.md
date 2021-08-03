# **Welcome to the tutorial**
<hr>
### **Contents**
#### 1. [Overview of the analysis](#Overview)
#### 2. [Differential protein expression for label-free quantification (DEP-LFQ)](#DEP-LFQ)
#### &emsp; 2.1 [Overview](#overview_DEP)
#### &emsp; 2.2 [Analytical steps](#analytical_steps_DEP)
#### &emsp; 2.3 [Output](#output_DEP)
#### &emsp;&emsp; 2.3.1 [Result tables](#output_DEP_result_tables)
#### &emsp;&emsp; 2.3.2 [Visualization of the results](#output_DEP_result_plot)
#### &emsp;&emsp;&emsp; 2.3.2.1 [Barplot of a protein of interest](#output_DEP_result_plot_single_protein)
#### &emsp;&emsp;&emsp; 2.3.2.2 [Heatmap of all significant proteins](#output_DEP_result_plot_heatmap)
#### &emsp;&emsp;&emsp; 2.3.2.3 [Volcano plots of specific contrasts](#output_DEP_result_plot_volcano)
#### &emsp;&emsp;&emsp; 2.3.2.4 [Customed volcano plots of specific contrasts](#output_DEP_result_plot_custom_volcano)
#### &emsp;&emsp;&emsp; 2.3.2.5 [PCA plot](#output_DEP_result_plot_pca)
#### &emsp;&emsp;&emsp; 2.3.2.6 [Correlation matrix](#output_DEP_result_plot_pearson)
#### &emsp;&emsp;&emsp; 2.3.2.7 [Distance matrix](#output_DEP_result_plot_gower)
#### &emsp;&emsp;&emsp; 2.3.2.8 [Sample CVs Plots](#output_DEP_result_plot_cv)
#### &emsp;&emsp;&emsp; 2.3.2.9 [Barplot of protein numbers](#output_DEP_result_plot_protein_number)
#### &emsp;&emsp;&emsp; 2.3.2.10 [Sample coverage](#output_DEP_result_plot_sample_coverage)
#### &emsp;&emsp;&emsp; 2.3.2.11 [Normalization plot](#output_DEP_result_plot_normalization)
#### &emsp;&emsp;&emsp; 2.3.2.12 [Missing values- Quant](#output_DEP_result_plot_missing_value_quant)
#### &emsp;&emsp;&emsp; 2.3.2.13 [Missing values- Heatmap](#output_DEP_result_plot_missing_value_heatmap)
#### &emsp;&emsp;&emsp; 2.3.2.14 [Imputation](#output_DEP_result_plot_imputation)
#### 3. [Gene annotation (Annotation)](#annotation)
#### &emsp; 3.1 [Overview](#overview_annotation)
#### &emsp; 3.2 [Analytical steps](#analytical_steps_annotation)
#### &emsp; 3.3 [Output](#output_annotation)
#### 4. [Over Representation Analysis (ORA)](#ORA)
#### &emsp; 4.1 [Overview](#overview_ORA)
#### &emsp; 4.2 [Analytical steps](#analytical_steps_ORA)
#### &emsp; 4.3 [Output](#output_ORA)
#### &emsp;&emsp; 4.3.1 [Result tables](#output_ORA_result_tables)
#### &emsp;&emsp; 4.3.2 [Visualization of the results](#output_ORA_result_plot)
#### &emsp;&emsp;&emsp; 4.3.2.1 [Bar plot](#output_ORA_result_plot_bar)
#### &emsp;&emsp;&emsp; 4.3.2.2 [Dot plot](#output_ORA_result_plot_dot)
#### &emsp;&emsp;&emsp; 4.3.2.3 [Optimized dot plot](#output_ORA_result_plot_dot_opt)
#### &emsp;&emsp;&emsp; 4.3.2.4 [Heatmap-like functional classification](#output_ORA_result_plot_heat)
#### &emsp;&emsp;&emsp; 4.3.2.5 [Gene-Concept Network](#output_ORA_result_plot_cnet)
#### &emsp;&emsp;&emsp; 4.3.2.6 [Enrichment Map](#output_ORA_result_plot_emap)
#### &emsp;&emsp;&emsp; 4.3.2.7 [goplot](#output_ORA_result_plot_goplot)
#### &emsp;&emsp;&emsp; 4.3.2.8 [GO graph](#output_ORA_result_plot_gograph)
#### 5. [Gene Set Enrichment Analysis (GSEA)](#GSEA)
#### &emsp; 5.1 [Overview](#overview_GSEA)
#### &emsp; 5.2 [Analytical steps](#analytical_steps_GSEA)
#### &emsp; 5.3 [Output](#output_GSEA)
#### &emsp;&emsp; 5.3.1 [Result tables](#output_GSEA_result_tables)
#### &emsp;&emsp; 5.3.2 [Visualization of the results](#output_GSEA_result_plot)
#### &emsp;&emsp;&emsp; 5.3.2.1 [Bar plot](#output_GSEA_result_plot_bar)
#### &emsp;&emsp;&emsp; 5.3.2.2 [Dot plot](#output_GSEA_result_plot_dot)
#### &emsp;&emsp;&emsp; 5.3.2.3 [Heatmap-like functional classification](#output_GSEA_result_plot_heat)
#### &emsp;&emsp;&emsp; 5.3.2.4 [Gene-Concept Network](#output_GSEA_result_plot_cnet)
#### &emsp;&emsp;&emsp; 5.3.2.5 [Enrichment Map](#output_GSEA_result_plot_emap)
#### &emsp;&emsp;&emsp; 5.3.2.6 [Running score and preranked list of GSEA result](#output_GSEA_result_plot_gseaplot)
#### 6. [String protein-protein interaction network (PPITools)](#PPI)
#### &emsp; 6.1 [Overview](#overview_PPI)
#### &emsp; 6.2 [Analytical steps](#analytical_steps_PPI)
#### &emsp; 6.3 [Output](#output_PPI)
#### &emsp;&emsp; 6.3.1 [Result table](#output_PPI_result_tables)
#### &emsp;&emsp; 6.3.2 [Visualization of the results](#output_PPI_result_plot)
#### &emsp;&emsp;&emsp; 6.3.2.1 [Network](#output_PPI_result_plot_network)






<a name="Overview"></a>
### **1. Overview of the analysis**
<p style="font-size:1.0em;text-align:justify">
This shiny application provides a workflow for mass spectrometry proteomics data(label-free quantification data) processed by MaxQuant. The full tutorial is available online and embedded within the application for users to access. Default parameters are provided to guide unfamiliar users through the analytical steps. Generally, the order of the analytical steps are from the top to the bottom, and default parameters are proper, also if desired, could be changed based on user preference. This application aims to build a whole workflow to help users to explore data.It supports:
</br>
</br>
&emsp; 1. &nbsp;
Differential analysis of proteomics data based on R package <a href="https://www.bioconductor.org/packages/release/bioc/html/limma.html">Limma</a> and <a href="https://www.bioconductor.org/packages/release/bioc/html/DEP.html">DEP</a> 
</br>

&emsp; 2. &nbsp;
Gene annotation including gene description, entrezid, the name of <a href="http://geneontology.org/">gene ontology</a> BP(biological process), CC(cellular component), MF(molecular function), the name of <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome pathway</a>, and the name of <a href="http://pfam.xfam.org/">pfam</a>
</br>

&emsp; 3. &nbsp;
ORA(over-representation analysis) including <a href="http://geneontology.org/">gene ontology</a>, <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome</a> pathway enrichment based on the most widely used and popular R package <a href="https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html">clusterProfiler</a>
</br>

&emsp; 4. &nbsp;
GSEA(Gene Set Enrichment Analysis) including <a href="http://geneontology.org/">gene ontology</a>, <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome</a> pathway enrichment based on the most widely used and popular R package <a href="https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html">clusterProfiler</a>
</br>

&emsp; 5. &nbsp;
PPI(protein–protein interactions) based on <a href="https://string-db.org/">STRING</a>
</br>

&emsp; 6. &nbsp;
Differential analysis of RNAseq countmatrix data based on R package <a href="https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html">DEseq2</a>
</br>

&emsp; 7. &nbsp;
Compare of differential analysis of proteomics and countmatrix of RNAseq data with a heatmap together
</br>

&emsp; 8. &nbsp;
some extra links for convenience
</br>
</br>
Finally, considering the flexibility, this appalication is designed:
</br>
</br>
&emsp; 1. &nbsp;
Every part could either integrate with the results of differential analysis of proteomics data and RNAseq countmatrix data, or run separately. It allows user to have more flexible options of analytical steps
</br>

&emsp; 2. &nbsp;
Every part has a series of parameters that could be modified to help experienced user. There will also appear the detailed explanation when you hover on a parameter to help understand
</br>

&emsp; 3. &nbsp;
The results of every part including tables and figures that help explore data could be downloaded. Each of the main figures has a series of parameters that could be modified to allow user to custom such as heatmap, volcano, pca, pearson correlation, Gowers's distance plot and so on. And every figure has width and height parameters to allow user to get proper and nice pdf file
</br>
</br>
Summary, this application provides a whole workflow for differential analysis of mass spectrometry proteomics data(label-free data) processed by MaxQuant, optional compared with differential analysis of RNAseq countmatrix data considering that the level of protein is not always consistent with the level of RNA with a heatmap together. It could also give user insight into gene annotation, <a href="http://geneontology.org/">gene ontology</a>, <a href="https://www.genome.jp/kegg/">KEGG</a> and <a href="https://reactome.org/">reactome</a> pathway enrichment including both ORA and GSEA, as well as protein–protein interactions to explore data such as biological mechanism. Flexibly, each part could either be an independent functional module or integrate with the results of differential analysis of proteomics data and RNAseq countmatrix data. Finally, each part has a series of parameters that could be modified and appear detailed explanations when you hover on them, and also has a series of results including tabels and figures to help explore data which could always be downloaded properly by modifying some parameters. So, no matter you want to do the differential analysis of proteomics data, or the differential analysis of RNAseq countmatrix data, or gene annotation, or ORA, or GSEA, or PPI separately, or integrate differential analysis of proteomics data and RNAseq countmatrix data with each other or the other parts, this application will be an ideal choice.
</p>
<hr>



<a name="DEP-LFQ"></a>
### **2. Differential protein expression for label-free quantification  (DEP-LFQ)**
<a name="overview_DEP"></a>
#### **2.1 Overview**
R package <a href="https://www.bioconductor.org/packages/release/bioc/html/limma.html">Limma</a> is the best generalized approach for intensity data such as those used in proteomics. And R package <a href="https://www.bioconductor.org/packages/release/bioc/html/DEP.html">DEP</a> provides an analysis workflow of mass spectrometry proteomics data for differential protein expression based on protein-wise linear models and empirical Bayes statistics using limma. We use these two packages to build the flexible and interactive shiny framework. The main steps include data upload, filtering, variance normalization and imputation of missing values, as well as statistical testing of differentially enriched / expressed proteins. Although R package DEP contains a shiny application, our framework is more flexible, elegant, ideal and friendly to users who don't know how to use R language. It provides more options for statistical test such as the threshold for the allowed number of missing values in at least one condition, the imputation method, the pvalue adjustment method, the type of contrasts. It also provides more figures such as pca plot and each figure has a series of parameters such as color helping custom and explore the data better. Finally, each figure could be downloaded as a proper pdf file.

<a name="analytical_steps_DEP"></a>
#### **2.2 Analytical steps**
`Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter
</br>
&emsp;1) Please go to the `DEP-LFQ options` on the sidebar
</br>
&emsp;2) Go to `Files`, and the file requirements and examples could be found when click on the question mark on the right
</br>
&emsp;&emsp;2.1) Upload file **ProteinGroups.txt** processed by MaxQuant
</br>
&emsp;&emsp;2.2) Optional, Upload file **ExperimentalDesign.txt** when the **Sample annotation** is set to **Use Experimental Design**
</br>
&emsp;3) Go to `Columns`
</br>
&emsp;&emsp;3.1) **Name column**, **ID column**, **Filter on columns** are set automatically
</br>
&emsp;&emsp;3.2) If desired, **allowed max.miss.num at least one condition** is allowed to set
</br>
&emsp;&emsp;3.3) `Contrasts` on the top-right of the body is allowed to set to control or all or manual. When it is set to control or manual, there will appear parameter **Control** or **Manual test** that allows to set
</br>
&emsp;4) Go to `imputation options`, and select an imputation type. The explanations of the imputation methods could be found when click on the **Detailed information**
</br>
&emsp;5) Go to `Type of FDR correction`, and select the method of pvalue adjustment
</br>
&emsp;6) Click on `Analyze` button
</br>
&emsp;7) A series of parameters are allowed to set such as **adj.P value**, **Log2 fold change** and so on
</br>
&emsp;8) Optional but Strongly recommended, click on **save result RData** on the sidebar to save the main result object in order to load the RData and ProteinGroups.txt together next time for reproducibility 

<a name="output_DEP"></a>
#### **2.3 Output**

<a name="output_DEP_result_tables"></a>
#### **2.3.1 Result tables**
Users could download the result tabels on the sidebar

- **results:** Result of all proteins. Includes Gene names, Protein Ids, p-values (each pairwise comparisons), Adjusted p-values (each pairwise comparisons), boolean values for significance (each pairwise comparisons, note that the significant column is TRUE if any pairwise comparisons is TRUE), Log2 fold changes/ ratios (each pairwise comparisons), centered protein intensity (log2 transformed) in each sample 
- **significant_proteins:** Result of significant proteins according to the adjusted pvalue and log2 fold change cutoff. The columns are the same as **results**
- **displayed_subset:** Result of the displayed subset according to parameters **Select direct comparisons** and **Exclude direct comparisons** on the top of the result table. The columns are the same as **results**
- **full_dataset:** Result of all proteins. The columns not only include columns of **results**, but also the columns of ProteinGroups.txt to allow user to see more information of proteins
- **Result RData:** Users could download the result RData on the sidebar

Strongly recommended, click on **save result RData** on the sidebar to save the main result object in order to load the RData and ProteinGroups.txt together next time for reproducibility. Because, some imputation methods not always get the same result such as MinProb

<a name="output_DEP_result_plot"></a>
#### **2.3.2 Visualization of the results**
Each plot has a series of parameters such as color, width and height helping custom and explore the data better, and could be downloaded as a proper pdf file. `Note:` The box containing parameters may be collapsed in order to save space. You can click on the plus sign (+) on the right of the box to expand the corresponding category

<a name="output_DEP_result_plot_single_protein"></a>
#### **2.3.2.1 Barplot of a protein of interest**
In the shiny application, it is named as **Selected Protein**. Barplot of a protein of interest. x-axis represents the conditions / groups, y-axis represents log2 centered intensity. Width and height could be modified in order to get a proper pdf file.
</br>
<img src="selected_protein.png" style="margin-left:0px;" width=40% />
</br>

<a name="output_DEP_result_plot_heatmap"></a>
#### **2.3.2.2 Heatmap of all significant proteins**
In the shiny application, it is named as **Heatmap**. The heatmap representation gives an overview of all significant/differentially expressed proteins (rows) in all samples (columns). It means that the rows of the heatmap are those proteins significant/differentially expressed in at least one contrast. This visualization allows the identfication of general trends such as if one sample or replicate is highly different compared to the others and might be considerd as an outlier. Additionally, the hierarchical clustering of samples (columns) indicates how related the different samples are and hierarchicalclustering of proteins (rows) identifies similarly behaving proteins. The proteins can be clustered by k-means clustering (kmeans argument) and the number of clusters can be defined by argument **Kmeans**. There are a series of parameters that could be modified to allow to get a customed figure. Here, we will give some examples. In order to get a heatmap that row names are readable and do not take up space, we set the log2 fold change cut off very large to 5.
</br>
- **example 1:** The default parameters
</br>
<img src="DEP_heatmap_exampe1.png" style="margin-left:0px;" width=40% />
</br>
- **example 2:** **Color** parameter could be modified to change the color limit of the colorbar. Here the color is lighter than example 1
</br>
<img src="DEP_heatmap_exampe2.png" style="margin-left:0px;" width=40% />
</br>
- **example 3:** If you want to order the clusters of the rows from top to bottom, you should check parameter **If mysplit** under the collapsed box, and set parameter **my split**. Also , you can set the order of the columns by unchecking **Cluster columns** and setting **Custom columns order**. Here, we set the order of clusters of the rows 1 to 6 from top to bottom, the order of the columns A to B from left to right
</br>
<img src="DEP_heatmap_exampe3.png" style="margin-left:0px;" width=40% />
</br>
- **example 4:** We provide parameter **rowname color** by checking it to show colored row names to check the number of peptides associated with each protein group(purple: Peptides == 1, blue: Peptides == 2). If no color of row names is purple or blue, it means that the number of peptides associated with all protein groups is lager than 2. So, user could check the accuracy of the quantification. The accuracy of the quantification of a protein having one or two peptides may be less than those having more than two peptides
</br>
<img src="DEP_heatmap_exampe4.png" style="margin-left:0px;" width=40% />
</br>
- **example 5:** If your contrasts is more than one, and want to show specified contrasts, the parameters **Contrast** and **Manual heatmap** are provided. And it supports one or more specified contrasts. One specified contrast represents to show significant/differentially expressed proteins (rows) in this contrast. Two or more specified contrasts represents to show significant/differentially expressed proteins (rows) in at least one of these specified contrasts. What you should do is that, select your interested contrasts by setting parameter **Contrast**, and then check **Manual heatmap**. Here we choose contrast B_vs_A from B_vs_A and C_vs_A
</br>
<img src="DEP_heatmap_exampe5.png" style="margin-left:0px;" width=40% />
</br>
- **example 6:** Color palettes from R package RColorBrewer are allowed to specify by setting parameter **colorbar**. Here is the color palette RdGy
</br>
<img src="DEP_heatmap_exampe6.png" style="margin-left:0px;" width=40% />
</br>
- **example 7:** Alternatively, a heatmap can be plotted using the contrasts by setting **Data presentation** to contrast
</br>
<img src="DEP_heatmap_exampe7.png" style="margin-left:0px;" width=40% />
</br>

<a name="output_DEP_result_plot_volcano"></a>
#### **2.3.2.3 Volcano plots of specific contrasts**
In the shiny application, it is named as **Volcano plot**. Volcano plots can be used to visualize a specific contrast (comparison between two samples). It is a graphical visualization by plotting the "Fold Change (Log2)" on the x-axis versus the –log10 of the "p-value" on the y-axis. The contrast could be set under **Contrast**. Protein labels could be added by checking parameter **Display names**. **Font size** could be set for labels. And it should be omitted when there are too many labels. Y-axis could be adjusted p values by checking parameter **Adjusted p values**. The width from 0 on x axis could be the same by checking **Same width**. And by selecting **Mybreaks** and checking **My breaks** under the collapsed box, the breaks of x axis could be customed. Here is an example
</br>
<img src="DEP_volcano_example.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_custom_volcano"></a>
#### **2.3.2.4 Customed volcano plots of specific contrasts**
In the shiny application, it is named as **Custom Volcano**. It is more flexible and suggested than **Volcano plot**. Here, we will give some examples. 
</br>
- **example 1:** The default parameters. By default, it shows the top 20 significant proteins
</br>
<img src="DEP_custom_volcano_example1.png" style="margin-left:0px;" width=40% />
</br>
- **example 2:** It could show the top 20 significant proteins of only up-regulated or down-regulated by setting **label way** to up or down. The example sets it to up
</br>
<img src="DEP_custom_volcano_example2.png" style="margin-left:0px;" width=40% />
</br>
- **example 3:** The black line around the labeled points can be removed by setting **point outside width** to 0
</br>
<img src="DEP_custom_volcano_example3.png" style="margin-left:0px;" width=40% />
</br>
- **example 4:** It also could show specified genes that you are interested in. First, choose your interested genes under parameter **selected proteins** which is under the collapsed box, and then set **label way** to selected proteins. Note, to save time, the order of the steps should not be reversed.
</br>
<img src="DEP_custom_volcano_example4.png" style="margin-left:0px;" width=40% />
</br>
- **example 5:** To be more readable, the labels could be surrounded by a rectangular by checking **label with rectangle** under the collapsed box
</br>
<img src="DEP_custom_volcano_example5.png" style="margin-left:0px;" width=40% />
</br>
- **example 6:** It could be no any labels or all labels of significant proteins by setting **show number** under the collapsed box to 0 or a lager number than the number of significant proteins. It is not recommended to show all labels when there are too many labels. Here is an example of no labels
</br>
<img src="DEP_custom_volcano_example6.png" style="margin-left:0px;" width=40% />
</br>
- **example 7:** Just like the heatmap plot, we provide parameter **peptide color** by checking it to show colored points to check the number of peptides associated with each protein group(purple: Peptides == 1, blue: Peptides == 2, you can also modify the color by yourself). If no color of points is purple or blue, it means that the number of peptides associated with all protein groups is lager than 2. So, user could check the accuracy of the quantification. The accuracy of the quantification of a protein having one or two peptides may be less than those having more than two peptides.
</br>
<img src="DEP_custom_volcano_example7.png" style="margin-left:0px;" width=40% />
</br>
- **example 8:** The color of up-regulated, down-regulated and not significant proteins can be customed by parameters **up**, **down** and **not significant**
</br>
<img src="DEP_custom_volcano_example8.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_pca"></a>
#### **2.3.2.5 PCA plot**
In the shiny application, it is named as **Pca plot**. Principal component analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset. This can be very useful to observe batch effects, such as clear differences between replicates. The more similar 2 samples are, the closer they cluster together.
</br>
<img src="DEP_pca_example1.png" style="margin-left:0px;" width=40% />
</br>
A square picture is allowed by checking **if square**
</br>
<img src="DEP_pca_exampe2.png" style="margin-left:0px;" width=40% />
</br>
Facet is allowed by selecting the third option of **Color and shape**, Here we select the third option Condition
</br>
<img src="DEP_pca_example3.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_pearson"></a>
#### **2.3.2.6 Correlation matrix**
In the shiny application, it is named as **Pearson correlation**. A correlation matrix can be plotted as a heatmap to visualize the Pearson correlations between the different samples. We provide parameters **lower** and **upper**. Setting the lower correlation limit can compress the plottable data and give greater contrast to the plot. It is rare that two replicates will vary greatly in their correlation. Here we set **lower** to 0.7 to see a good contrast.
</br>
<img src="DEP_pearson_example1.png" style="margin-left:0px;" width=40% />
</br>
You can display the correlation values on the plot by checking **Add values** under the collapsed box. The decimal point digits and font size of the values also can be set by **Value digits** and **Value size**. Here we set the decimal point digits of the values to 4.
</br>
<img src="DEP_pearson_example2.png" style="margin-left:0px;" width=40% />
</br>
The palette of the color could be changed by **color panel** which provides the palettes from R package RColorBrewer, and if you want to get a reversed colorbar, please check **pal rev**. Here, we set the palette to Blues.
</br>
<img src="DEP_pearson_example3.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_gower"></a>
#### **2.3.2.7 Distance matrix**
In the shiny application, it is named as **Gower's distance**. A distance matrix can be plotted as a heatmap to visualize the Gower's distances between the different samples. 
</br>
<img src="DEP_gower_example1.png" style="margin-left:0px;" width=40% />
</br>
You can display the distance values on the plot by checking **Add values** under the collapsed box. The decimal point digits and font size of the values also can be set by **Value digits** and **Value size**. Here we set the decimal point digits of the values to 2.
</br>
<img src="DEP_gower_example2.png" style="margin-left:0px;" width=40% />
</br>
The palette of the color could be changed by **color panel** which provides the palettes from R package RColorBrewer, and if you want to get a reversed colorbar, please check **pal rev**. Here, we set the palette to Reds.
</br>
<img src="DEP_gower_example3.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_cv"></a>
#### **2.3.2.8 Sample CVs Plots**
In the shiny application, it is named as **Sample CVs**. It is a plot representing distribution of protein level coefficient of variation for each condition. Each plot also contains a vertical line representing median CVs percentage within that condition.
</br>
<img src="DEP_cvplot_example.png" style="margin-left:0px;" width=60% />
</br>


<a name="output_DEP_result_plot_protein_number"></a>
#### **2.3.2.9 Barplot of protein numbers**
In the shiny application, it is named as **Protein Numbers**. It is a barplot representing number of proteins identified and quantified in each sample.
</br>
<img src="DEP_protein_number_example.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_sample_coverage"></a>
#### **2.3.2.10 Sample coverage**
In the shiny application, it is named as **Sample coverage**. It is a plot highlighting overlap between identified proteins across all samples in the experiment.
</br>
<img src="DEP_protein_coverage_example.png" style="margin-left:0px;" width=20% />
</br>


<a name="output_DEP_result_plot_normalization"></a>
#### **2.3.2.11 Normalization plot**
In the shiny application, it is named as **Normalization**. The data is normalized by variance stabilizing transformation (vsn). The normalization can be inspected by checking the distributions of the samples before and after normalization.
</br>
<img src="DEP_normalization_example.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_missing_value_quant"></a>
#### **2.3.2.12 Missing values- Quant**
In the shiny application, it is named as **Missing values- Quant**. the densities and cumulative fractions are plotted for proteins with and without missing values to check whether missing values are biased to lower intense proteins.
</br>
<img src="DEP_missing_values_quant_example.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_missing_value_heatmap"></a>
#### **2.3.2.13 Missing values- Heatmap**
In the shiny application, it is named as **Missing values- Heatmap**. Only proteins with at least one missing value are visualized. It could assist users to identify patterns of missing values (i.e., MAR: Missing values at random or MNAR: Missing values not at random) and guide the choice of imputation method, if desired. If you want to order the columns by yourself, please select your ordered column names through parameter **Custom columns order**, and then uncheck **Cluster columns**
</br>
<img src="DEP_missing_value_heatmap_example.png" style="margin-left:0px;" width=40% />
</br>


<a name="output_DEP_result_plot_imputation"></a>
#### **2.3.2.14 Imputation**
In the shiny application, it is named as **Imputation**. The effect of the imputation on the distributions can be visualized through a desity plot of protein intensity (log2) distrubution for each condition after and before missing value imputation being performed.
</br>
<img src="DEP_imputation_example.png" style="margin-left:0px;" width=40% />
</br>


<hr>
<a name="annotation"></a>
### **3. Gene annotation (Annotation)**
<a name="overview_annotation"></a>
#### **3.1 Overview**
It is usually necessary to do a gene annotation analysis when you get a list of differentially enriched / expressed genes or genes that you are interested in. So you could know more information about these genes so as to understand the relevant biological significance. Here we provide the annotations of gene description, entrezid, the name of [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), the name of [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway, and the name of [pfam](http://pfam.xfam.org/). SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported as input. And now it supports organisms of human and mouse. We support not only pasting a gene list directly, but also importing the result of DEP-LFQ or DEG-RNAseq part to make this application flexible. Finally, we provide the result table of gene annotation to allow user to download. Note that, genes that have no entrezids will not appear in the result table. 


<a name="analytical_steps_annotation"></a>
#### **3.2 Analytical steps**
&emsp;1) Please go to the `Annotation options` on the sidebar
</br>
&emsp;2) You can choose to paste a gene list directly or import the result of DEP-LFQ or DEG-RNAseq part
</br>
&emsp;&emsp;2.1) If you have a gene list that you are interested in, please paste them directly in the box under **Please paste your gene list**
</br>
&emsp;&emsp;2.2) If you want to import the result of DEP-LFQ or DEG-RNAseq part, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part, and then click on the menuitem **Import from** and check **DEP-LFQ or DEG-RNAseq**. It will appears some options for you to choose:

- **UPregu for DEP-LFQ:** up-regulated genes for DEP-LFQ panel
- **DOWNregu for DEP-LFQ:** down-regulated genes for DEP-LFQ panel
- **UPDOWN for DEP-LFQ:** all regulated genes, up- and down- regulated for DEP-LFQ panel
- **UPregu for DEG-RNAseq:** up-regulated genes for DEG-RNAseq panel
- **DOWNregu for DEG-RNAseq:** down-regulated genes for DEG-RNAseq panel
- **UPDOWN for DEG-RNAseq:** all regulated genes, up- and down- regulated for DEG-RNAseq panel

Then, you need to choose the contrast that you are interested in by parameter **Contrast**. `Note that:` [Any significant] represents the genes that are significant in at least one contrast
</br>
&emsp;3) Select the organism
</br>
&emsp;4) Click on the **Analyze** button
</br>
&emsp;5) A table of the result of gene annotation will appear and allow user to download through the **Save table** button on the sidebar

<a name="output_annotation"></a>
#### **3.3 Output**
Users could download the result tabels on the sidebar
- **Result table:** Result of gene annotation. Includes gene name, gene description, entrezid, the name of [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), the name of [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway, and the name of [pfam](http://pfam.xfam.org/)






<hr>
<a name="ORA"></a>
### **4. Over Representation Analysis (ORA)**
<a name="overview_ORA"></a>
#### **4.1 Overview**
Over Representation Analysis (ORA) is a widely used approach to determine whether known biological functions or processes are over-represented (= enriched) in an experimentally-derived gene list, e.g. a list of differentially expressed genes (DEGs). For further information, here is a [link](http://yulab-smu.top/clusterProfiler-book/chapter2.html#over-representation-analysis). We use the most widely used and popular R package [clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) to perform over representation analysis. It needs a gene list (eg: DEGs) as input, and SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported. It supports [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway enrichment. And now it supports organisms of human and mouse. It provides tables to allow to download, a series of figures and parameters that allow user to change. And we support not only pasting a gene list directly, but also importing the result of DEP-LFQ or DEG-RNAseq part to make this application flexible. Except that gene ontology enrichment has a few more parameters such as **Ontology** and **removed redundancy of enriched GO terms** and figures output **Goplot** and **Gograph**, their design and analytical steps are exactly the same. So, in this tutorial, we will just take gene ontology enrichment as an example. Finally, each figure could be downloaded as a proper pdf file.

<a name="analytical_steps_ORA"></a>
#### **4.2 Analytical steps**
`Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter 
</br>
&emsp;1) Please go to the `ORA options` on the sidebar, and choose the type of analysis `gene ontology (GO)` , `KEGG`, or `reactome` pathway enrichment that you want to perform
</br>
&emsp;2) You can choose to paste a gene list directly or import the result of DEP-LFQ or DEG-RNAseq part
</br>
&emsp;&emsp;2.1) If you have a gene list that you are interested in, please paste them directly in the box under **Please paste your gene list**
</br>
&emsp;&emsp;2.2) If you want to import the result of DEP-LFQ or DEG-RNAseq part, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part, and then click on the menuitem **Import from** and check **DEP-LFQ or DEG-RNAseq**. It will appears some options for you to choose:

- **UPregu for DEP-LFQ:** up-regulated genes for DEP-LFQ panel
- **DOWNregu for DEP-LFQ:** down-regulated genes for DEP-LFQ panel
- **UPDOWN for DEP-LFQ:** all regulated genes, up- and down- regulated for DEP-LFQ panel
- **UPregu for DEG-RNAseq:** up-regulated genes for DEG-RNAseq panel
- **DOWNregu for DEG-RNAseq:** down-regulated genes for DEG-RNAseq panel
- **UPDOWN for DEG-RNAseq:** all regulated genes, up- and down- regulated for DEG-RNAseq panel

Then, you need to choose the contrast that you are interested in by parameter **Contrast**. `Note that:` [Any significant] represents the genes that are significant in at least one contrast
</br>
&emsp;3) Select the organism
</br>
&emsp;4) check or uncheck **If with log2 fold change**, it means that whether your gene list contains log2 fold change (total two columns), and it is used only in **Cnetplot** and **Emaplot** to fill the color
</br>
&emsp;5) Click on the **Analyze** 
</br>
&emsp;6) A series of parameters are allowed to set such as **Ontology**, **adj. P value** and so on




<a name="output_ORA"></a>
#### **4.3 Output**
<a name="output_ORA_result_tables"></a>
#### **4.3.1 Result tables**
Users could download the result tabels on the sidebar
- **full_results:** Result of all terms
- **significant_results:** Result of significant terms according to the pvalue, adjusted pvalue and qvalue cutoff. Generally speaking, the default parameters are appropriate


<a name="output_ORA_result_plot"></a>
#### **4.3.2 Visualization of the results**
For every figure, you can change **width** and **height** in order to get a proper pdf file.
<a name="output_ORA_result_plot_bar"></a>
#### **4.3.2.1 Bar plot**
In the shiny application, it is named as **Bar plot**. Bar plot is the most widely used method to visualize enriched terms. Here, it depicts the enrichment scores (p values or adjusted p values) and gene count as bar color and height. You can change the bar color by parameter **colorBy**. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
<img src="ORA_barplot_example1.png" style="margin-left:0px;" width=60% />
</br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together
</br>
<img src="ORA_barplot_example2.png" style="margin-left:0px;" width=60% />
</br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each BP, CC and MF
</br>
<img src="ORA_barplot_example3.png" style="margin-left:0px;" width=60% />
</br>



<a name="output_ORA_result_plot_dot"></a>
#### **4.3.2.2 Dot plot**
In the shiny application, it is named as **Dot plot**. Dot plot is similar to bar plot with the capability to encode another score as dot size. Here, the x-axis represents gene ratio, the color represents p values or adjusted p values by setting **colorBy**, and the dot size represents gene count. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
<img src="ORA_dotplot_example1.png" style="margin-left:0px;" width=60% />
</br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together
</br>
<img src="ORA_dotplot_example2.png" style="margin-left:0px;" width=60% />
</br>

If your selected analysis type of ORA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each BP, CC and MF
</br>
<img src="ORA_dotplot_example3.png" style="margin-left:0px;" width=60% />
</br>


<a name="output_ORA_result_plot_dot_opt"></a>
#### **4.3.2.3 Optimized dot plot**
In the shiny application, it is named as **Dot opt**. It is simplified relative to the dot plot. The color represents p values or adjusted p values by setting **colorBy**, and the dot size represents gene count. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
<img src="ORA_dot_opt_example.png" style="margin-left:0px;" width=80% />
</br>


<a name="output_ORA_result_plot_heat"></a>
#### **4.3.2.4 Heatmap-like functional classification**
In the shiny application, it is named as **Heatplot**. It is similar to cnetplot, while displaying the relationships as a heatmap. The gene-concept network may become too complicated if user want to show a large number significant terms. The heatplot can simplify the result and more easy to identify expression patterns. The rows represent enriched items, the columns represent genes, and the color represents the log2 fold change when parameter **If with log2 fold change** is checked. From this plot, you can see genes of your gene list that the significant items include. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
This is an example that parameter **If with log2 fold change** is unchecked
</br>
</br>
</br>
<img src="ORA_heatplot_example1.png" style="margin-left:0px;" width=80% />
</br>
</br>
</br>

This is an example that parameter **If with log2 fold change** is checked
</br>
</br>
</br>
<img src="ORA_heatplot_example2.png" style="margin-left:0px;" width=80% />
</br>
</br>
</br>


<a name="output_ORA_result_plot_cnet"></a>
#### **4.3.2.5 Gene-Concept Network**
In the shiny application, it is named as **Cnetplot**. Both the barplot and dotplot only displayed most significant enriched terms, while users may want to know which genes are involved in these significant terms. The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network. The yellow point represents items, the colors of the gene points are filled with log2 fold change when parameter **If with log2 fold change** is checked. If a gene belongs to a item, then there will be a line between them. You can change the number of top significant terms to show by setting **ShowCategory**
</br>
This is an example that parameter **If with log2 fold change** is unchecked
</br>
</br>
</br>
<img src="ORA_cnetplot_example1.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>

This is an example that parameter **If with log2 fold change** is checked
</br>
</br>
</br>
<img src="ORA_cnetplot_example2.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>

This is an example that parameter **Circular** is unchecked to change the layout
</br>
</br>
</br>
<img src="ORA_cnetplot_example3.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>


<a name="output_ORA_result_plot_emap"></a>
#### **4.3.2.6 Enrichment Map**
In the shiny application, it is named as **Emaplot**. Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module. You can change the number of top significant terms to show by setting **ShowCategory**
</br>
</br>
</br>
<img src="ORA_emaplot_example.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>

<a name="output_ORA_result_plot_goplot"></a>
#### **4.3.2.7 goplot**
In the shiny application, it is named as **Goplot**. It can accept output of enrichGO and visualized the enriched GO induced graph. It shows the hierarchical relationships of terms
</br>
</br>
</br>
<img src="ORA_goplot_example.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>


<a name="output_ORA_result_plot_gograph"></a>
#### **4.3.2.8 GO graph**
In the shiny application, it is named as **GOgraph**. It is another form of showing hierarchical relationships of terms. `Notice that`, it has its own action button, please click on the **plot** button when you want to plot it
</br>
</br>
</br>
<img src="ORA_gograph_example.png" style="margin-left:0px;" width=80% />






<hr>
<a name="GSEA"></a>
### **5. Gene Set Enrichment Analysis (GSEA)**
<a name="overview_GSEA"></a>
#### **5.1 Overview**
Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states(e.g. phenotypes). For further information, here are some links, [link1](http://www.gsea-msigdb.org/gsea/index.jsp), [link2](http://yulab-smu.top/clusterProfiler-book/chapter2.html#gene-set-enrichment-analysis). And the difference between GSEA and ORA is that:

- ORA mainly focuses on a few genes such as genes that are significantly up-regulated or down-regulated. So it is easy to miss some genes that are not significantly differentially expressed but have important biological significance.
- GSEA is not limited to significantly differentially expressed genes, all genes can be used in GSEA. From the perspective of gene set enrichment, it is theoretically easier to include the impact of slight but coordinated changes on biological pathways, especially when gene sets with not too large differential multiples.

We use the most widely used and popular R package [clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) to perform gene set enrichment analysis. It needs a gene list containing log2 fold change (total two columns) as input, and SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported. It supports [gene ontology](http://geneontology.org/) BP(biological process), CC(cellular component), MF(molecular function), [KEGG](https://www.genome.jp/kegg/) and [reactome](https://reactome.org/) pathway enrichment. And now it supports organisms of human and mouse. It provides tables to allow to download, a series of figures and parameters that allow user to change. And we support not only pasting a gene list containing log2 fold change (total two columns) directly, but also importing the result of DEP-LFQ or DEG-RNAseq part to make this application flexible. Except that gene ontology enrichment has a few more parameters such as **Ontology** and **removed redundancy of enriched GO terms**, their design and analytical steps are exactly the same. So, in this tutorial, we will just take gene ontology enrichment as an example. Finally, each figure could be downloaded as a proper pdf file.


<a name="analytical_steps_GSEA"></a>
#### **5.2 Analytical steps**
`Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter 
</br>
&emsp;1) Please go to the `GSEA options` on the sidebar, and choose the type of analysis `gseGO` , `gseKEGG`, or `gseReactome` pathway enrichment that you want to perform
</br>
&emsp;2) You can choose to paste a gene list containing log2 fold change (total two columns) directly or import the result of DEP-LFQ or DEG-RNAseq part
</br>
&emsp;&emsp;2.1) If you have a gene list containing log2 fold change (total two columns) that you are interested in, please paste them directly in the box under **Please paste your gene list**
</br>
&emsp;&emsp;2.2) If you want to import the result of DEP-LFQ or DEG-RNAseq part, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part, and then click on the menuitem **Import from** and check **DEP-LFQ or DEG-RNAseq**. It will appears some options for you to choose:

- **DEP-LFQ:** import from DEP-LFQ panel
- **DEG-RNAseq:** import from DEG-RNAseq panel

Then, you need to choose the contrast that you are interested in by parameter **Contrast**.
</br>
&emsp;3) Select the organism
</br>
&emsp;4) Click on the **Analyze** 
</br>
&emsp;5) A series of parameters are allowed to set such as **Ontology**, **adj. P value**, **NES**, **Phenotype** and so on
 

<a name="output_GSEA"></a>
#### **5.3 Output**
<a name="output_GSEA_result_tables"></a>
#### **5.3.1 Result tables**
Users could download the result tabels on the sidebar
- **full_results:** Result of all terms
- **significant_results:** Result of significant terms according to the pvalue, adjusted pvalue and NES (the absolute value of normalized enrichment score) cutoff. Generally speaking, the default parameters are appropriate


<a name="output_GSEA_result_plot"></a>
#### **5.3.2 Visualization of the results**
The parameter **Phenotype** means that the phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed. Default is both activated and suppressed
</br>
For every figure, you can change **width** and **height** in order to get a proper pdf file.
<a name="output_GSEA_result_plot_bar"></a>
#### **5.3.2.1 Bar plot**
In the shiny application, it is named as **Bar plot**. Bar plot is the most widely used method to visualize enriched terms. Here, it depicts the enrichment scores (p values or adjusted p values) and NES as bar color and height. You can change the bar color by parameter **colorBy**. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed
</br>
<img src="GSEA_barplot_example1.png" style="margin-left:0px;" width=50% />
</br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to activated
</br>
<img src="GSEA_barplot_example2.png" style="margin-left:0px;" width=50% />
</br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to suppressed
</br>
<img src="GSEA_barplot_example3.png" style="margin-left:0px;" width=50% />
</br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed
</br>
<img src="GSEA_barplot_example4.png" style="margin-left:0px;" width=50% />
</br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each phenotype that you selected based on each BP, CC and MF. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed
</br>
<img src="GSEA_barplot_example5.png" style="margin-left:0px;" width=60% />
</br>



<a name="output_GSEA_result_plot_dot"></a>
#### **5.3.2.2 Dot plot**
In the shiny application, it is named as **Dot plot**. Dot plot is similar to bar plot with the capability to encode another score as dot size. Here, the x-axis represents NES, the color represents p values or adjusted p values by setting **colorBy**, and the dot size represents gene count. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed
</br>
<img src="GSEA_dotplot_example1.png" style="margin-left:0px;" width=50% />
</br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to activated
</br>
<img src="GSEA_dotplot_example2.png" style="margin-left:0px;" width=50% />
</br>

Here is an example that parameter **ShowCategory** is set to 10, and **Phenotype** is set to suppressed
</br>
<img src="GSEA_dotplot_example3.png" style="margin-left:0px;" width=50% />
</br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By default, the plot shows the top **ShowCategory** number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed
</br>
<img src="GSEA_dotplot_example4.png" style="margin-left:0px;" width=50% />
</br>

If your selected analysis type of GSEA is gene ontology(GO) enrichment, you can show terms of BP, CC, MF together by setting **Ontology** to ALL. By checking **If splited by ontology**, you can show the top **ShowCategory** number of each phenotype that you selected based on each BP, CC and MF. Here parameter **ShowCategory** is set to 10, and **Phenotype** is set to both activated and suppressed
</br>
<img src="GSEA_dotplot_example5.png" style="margin-left:0px;" width=60% />
</br>


<a name="output_GSEA_result_plot_heat"></a>
#### **5.3.2.3 Heatmap-like functional classification**
In the shiny application, it is named as **Heatplot**. It is similar to cnetplot, while displaying the relationships as a heatmap. The gene-concept network may become too complicated if user want to show a large number significant terms. The heatplot can simplify the result and more easy to identify expression patterns. The rows represent enriched items, the columns represent genes, and the color represents the log2 fold change. From this plot, you can see genes of your gene list that the significant items include. And you can change the number of top significant terms to show by setting **ShowCategory**
</br>
</br>
</br>
<img src="GSEA_heatplot_example.png" style="margin-left:0px;" width=100% />
</br>
</br>
</br>


<a name="output_GSEA_result_plot_cnet"></a>
#### **5.3.2.4 Gene-Concept Network**
In the shiny application, it is named as **Cnetplot**. Both the barplot and dotplot only displayed most significant enriched terms, while users may want to know which genes are involved in these significant terms. The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network. The yellow point represents items, the colors of the gene points are filled with log2 fold change. If a gene belongs to a item, then there will be a line between them. You can change the number of top significant terms to show by setting **ShowCategory**
</br>
</br>
</br>
<img src="GSEA_cnetplot_example1.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>

This is an example that parameter **Circular** is unchecked to change the layout
</br>
</br>
</br>
<img src="GSEA_cnetplot_example2.png" style="margin-left:0px;" width=60% />
</br>
</br>
</br>


<a name="output_GSEA_result_plot_emap"></a>
#### **5.3.2.5 Enrichment Map**
In the shiny application, it is named as **Emaplot**. Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module. You can change the number of top significant terms to show by setting **ShowCategory**
</br>
</br>
</br>
<img src="GSEA_emaplot_example.png" style="margin-left:0px;" width=50% />
</br>
</br>
</br>


<a name="output_GSEA_result_plot_gseaplot"></a>
#### **5.3.2.6 Running score and preranked list of GSEA result**
In the shiny application, it is named as **Gseaplot**. Running score and preranked list are traditional methods for visualizing GSEA result. It is a graphical visualization of the distribution of the gene set and the enrichment score. You can select the term that you want to show by parameter **Term**. The pvalue, adjusted pvalue and NES are displayed as a table on the plot 
</br>
</br>
</br>
<img src="GSEA_gseaplot_example.png" style="margin-left:0px;" width=50% />
</br>
</br>
</br>








<hr>
<a name="PPI"></a>
### **6. String protein-protein interaction network (PPITools)**
<a name="overview_PPI"></a>
#### **6.1 Overview**
Proteins and their functional interactions form the backbone of the cellular machinery. Their connectivity network needs to be considered for the full understanding of biological phenomena. [String](https://string-db.org/) is a curated database with known and predicted protein-protein interactions and detailed description of their software can be found in their latest [manuscript](https://pubmed.ncbi.nlm.nih.gov/30476243/). We perform protein–protein interactions based on STRING. It needs a gene list as input, and SYMBOL, ENSEMBL, UNIPROT and ALIAS ids are supported. Now it supports organisms of human and mouse. It provides table and network to allow to download, and a series of parameters that allow user to change. We support not only pasting a gene list directly, but also importing the result of DEP-LFQ or DEG-RNAseq part to make this application flexible.


<a name="analytical_steps_PPI"></a>
#### **6.2 Analytical steps**
`Note:` For detailed explanation of each parameter, please see the explanation that appears when the mouse is hovering over the parameter 
</br>
&emsp;1) Please go to the `PPITools options` on the sidebar
</br>
&emsp;2) You can choose to paste a gene list directly or import the result of DEP-LFQ or DEG-RNAseq part
</br>
&emsp;&emsp;2.1) If you have a gene list that you are interested in, please paste them directly in the box under **Please paste your gene list**
</br>
&emsp;&emsp;2.2) If you want to import the result of DEP-LFQ or DEG-RNAseq part, please ensure that you have completed the analysis of DEP-LFQ or DEG-RNAseq part, and then click on the menuitem **Import from** and check **DEP-LFQ or DEG-RNAseq**. It will appears some options for you to choose:

- **UPregu for DEP-LFQ:** up-regulated genes for DEP-LFQ panel
- **DOWNregu for DEP-LFQ:** down-regulated genes for DEP-LFQ panel
- **UPDOWN for DEP-LFQ:** all regulated genes, up- and down- regulated for DEP-LFQ panel
- **UPregu for DEG-RNAseq:** up-regulated genes for DEG-RNAseq panel
- **DOWNregu for DEG-RNAseq:** down-regulated genes for DEG-RNAseq panel
- **UPDOWN for DEG-RNAseq:** all regulated genes, up- and down- regulated for DEG-RNAseq panel

Then, you need to choose the contrast that you are interested in by parameter **Contrast**. `Note that:` [Any significant] represents the genes that are significant in at least one contrast
</br>
&emsp;3) Select the organism
</br>
&emsp;4) Select concerned scores, it means that select which type of evidence will contribute to the prediction of the score under the active interaction sources
</br>
&emsp;5) Set the minimum required interaction score by parameter **scores cutoff**. The minimum required interaction score puts a threshold on the confidence score, such that only interaction above this score are included in the predicted network. Confidence limits are as follows:

- low confidence - 0.15 (or better)
- medium confidence - 0.4
- high confidence - 0.7
- highest confidence - 0.9

&emsp;6) Click on the **String it!** 
</br>
&emsp;7) A series of parameters are allowed to set such as **input node color**, **input line color**, **choose node shape**, **choose layout style of network** and so on




#### 6. [](#)
#### &emsp; 6.1 [](#)
#### &emsp; 6.2 [](#)
#### &emsp; 6.3 [Output](#output_PPI)
#### &emsp;&emsp; 6.3.1 [Result table](#output_PPI_result_tables)
#### &emsp;&emsp; 6.3.2 [Visualization of the results](#output_PPI_result_plot)
#### &emsp;&emsp;&emsp; 6.3.2.1 [Network](#output_PPI_result_plot_network)


 