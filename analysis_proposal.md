## Introduction

This project investigates how the impact of knockout (KO) of the Tet2 gene in conjunction with a high fat western diet (WD) influences hepatocellular carcinoma (HCC) in mice. By using single cell RNA sequencing data, I aim to characterize the specific changes in cell type composition and biological pathways that drive HCC in Tet2 deficient mice. Achieving this requires us to first preprocess the data, conducting clustering, and annotate the cell types using machine learning. Following this, we will conduct downstream analysis such as finding the statistically significant changes in cell type proportions, upregulated biological processes (gene ontology analysis), and also correlating different biological pathways, cell type proportions, and cell type states (such as T cell exhaustion and senescence) against the clinical annotated tumor scores for each of our mouse samples. 


## Datasets


The primary dataset is single cell RNA sequencing data (scRNA-seq) from 16 male 6 month old mice with different combinations tet2 KO and diet provided by our collaborators. There are 4 mice with tet2 KO/WD (TW), 4 are tet2 KO/normal diet (TN), 4 are wildtype/WD (WW) and 4 are wildtype/normal diet (WN). Sequencing was performed using 10x Genomics Chromium droplet technology and Illumina short read sequencing on liver tissue filtered out for hepatocytes. The raw count matrix for the mice contains 54,860 cells and 32,285 genes. The collaborators have also provided hepatocellular carcinoma scores (HCC) ranging from 1-6 in severity based on the size and quantity of lesions on the liver samples.

Additionally, we utilize a scRNA-seq mouse atlas dataset in order to train a cell type classification model. This atlas is from livercellatlas.org and is taken from the livers of mice with non alcoholic fatty liver disease (Guilliams et al. 2022). This atlas comprises 64,244 cells across 24 distinct cell types, including various macrophage subsets (LAMs, KCs, etc) and lymphoid populations (NK, NKT, specialized T cells, etc). Finally, functional enrichment will be done using public gene sets from the Bader Lab. All datasets are currently on the HPC3 cluster which means that the dataset feasibility is already addressed.


## Analyses


Data analysis begins with a quality control (QC) and preprocessing pipeline using Scanpy (Wolf et al. 2018). To exclude multiplets and empty droplets, I will calculate QC metrics using scanpy.pp.calculate_qc_metrics() and filter cells based on extreme outliers in the distribution of expressed features (number of genes) and total counts (by viewing scanpy.pl.violin() plot). Additionally, cells with a high percentage of mitochondrial gene expression ("MT-" prefix) will be removed to exclude low quality cells. Genes expressed in fewer than a minimum threshold of cells will also be removed. 

After QC we will normalize the expression data for cells to a total of 10,000 reads per cell and then logarithmize the data using the normalize_total() and and log1p() functions respectively. We then calculate the highly variable genes using the highly_variable_genes() functions, run principal component analysis (PCA) to do dimensionality reduction on the gene expression data using pca() and finally compute neighbors of cells using the neighbors() function. In order to visualize the clusters we use the umap() function. We can color the UMAP by batch date and we will look to make sure that cells are not clustering by batch date (since this is a sign of a batch effect artifact). 

The next step is to do cell type annotation. Here we will use CellTypist, a linear regression model that uses feature selection and train and test on an 80/20 split using the mentioned atlas dataset (Dominguez Conde et al. 2022). Once the model is trained and the performance is evaluated we will run it on our data to annotate each of the cells. 

After annotating our cell types we can conduct statistical analyses to see if there are any statistically significant differences between the cell type proportions for any of the cell types across the 4 experimental conditions. This pairwise statistical analysis can be done using ggplot in R (Wickham 2016). We can also use gene ontology analysis by calculating the highly variable genes for the cells and then using the AUCell and ssGSEA tools along with the previously mentioned Bader lab gene sets to determine the scores for all biological pathways in all the cells. The last type of analysis we can do is doing more specific cell type labeling to label the specific states such as exhaustion for CD8+ T cells or senescence. We can then do correlational analysis using scipy.stats to see which of these cell states/pathways correlate with increased HCC score. 


## Visualizations

Project visualizations will include both diagnostic (QC) plots and preprocessing plots as well as plots for biological interpretation. Initial diagnostic plots will include violin plots generated via scanpy.pl.violin() to visualize the distribution of mitochondrial RNA percentages and unique gene features per cell. These plots are needed for establishing the filtering thresholds used to remove low quality data. We will also generate an area under the receiver operator curve plot (AUROC) to determine the predictive quality of the cell type classification model (CellTypist). 

For biological interpretation I will use stacked bar plots using ggplot in R to show cell type composition in the experimental conditions as well as lines between conditions to show statistically significant differences. Heatmaps created with scanpy.pl.heatmap() will illustrate differential activity in biological pathways identified through gene ontology scoring. Finally, I will produce scatter plots with regression lines and Pearson/Spearman correlation coefficients to visualize the relationship between specific cellular states (such as T-cell exhaustion) and clinical Hepatocellular Carcinoma (HCC) scores.


## Conclusions

Overall the proposed steps are relatively standard and can always be tweaked as I continue to read the existing literature and discuss with my collaborators. Given my prior experience with single-cell RNA sequencing workflows and the fact that all datasets are already acquired and preprocessed on the HPC3 cluster I believe this project is highly feasible within the current academic quarter. The completed project will consist of a reproducible analysis pipeline and a variety of figures characterizing the transcriptomic landscape of Tet2 deficient HCC. This study aims to provide a robust computational foundation for understanding how high-fat diets contribute to the tumor progression at a single-cell resolution.


## References

Guilliams M, Bonnardel J, Haest B, Daems B, Eguchi C, Gola A, Scott CL. 2022. Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches. Cell 185: 379–396. 

Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19: 15. 

Domínguez Conde C, Xu C, Garcia-Alonso L, Pembroke J, Tan S, Itzkovitz S, Teichmann SA. 2022. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 375: eabl5197.

Aibar S, González-Blas CB, Moerman T, Huynh-Thu VA, Imrichova H, Hulselmans G, Aerts S. 2017. SCENIC: single-cell regulatory network inference and clustering. Nature Methods 14: 1083–1086.

Wickham H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
