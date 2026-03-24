# Introduction

Clonal hematopoiesis of indeterminate potential (CHIP) is a phenomenon closely associated with age where hematopoietic stem cells (HSCs) acquire mutations leading to an increased risk of cardiovascular disorders and cancer (Sikking et al. 2023). This condition is driven by mutations in a few key genes including DMNT3A and TET2 (Marnell et al. 2022). Recently there have been studies that investigate how CHIP and high fat diets in mice contribute to disease but more research must be done to elucidate the mechanism (Wong et al. 2023). In this project we focus on how CHIP in conjunction with a western diet (WD) in mice influences hepatocellular carcinoma (HCC) using single cell sequencing data. In particular we will characterize some of the cellular processes and patterns seen in CHIP related HCC. This project aims to get a better understanding of how this phenomenon occurs in mice to get actionable insights for humans. 


# Data and Preprocessing

The goal of this project is to understand the patterns and processes involved in CHIP related HCC to develop testable wet-lab hypotheses to get a mechanistic understanding. The data for this project is primarily single-cell RNA sequencing data for 16 different male 6 month old mice. The single cell data is a matrix where the rows are unique IDs assigned to each cell and the columns are genes. Each entry is the expression level for a particular gene for one cell. These 16 mice are split into 4 equal groups with two different conditions. The four groups are TN, TW, WN, and WW where the first letter is a T if it is tet2 KO and a W if it is wildtype and the second letter is a N or W based on if the mice was fed a normal or western diet. Additionally, we are provided with HCC scores for each of the 16 samples based on how many tumor nodules are present on each of the mice's livers. We also have variant allele frequencies (VAF) for tet2 knockout (KO) based on ddPCR data for each of the mice samples. The first step is to process the data by removing low quality cells like multiplets and cells with high amounts of mitochondrial DNA. Next we need to determine the cell types for each of the cells in our data. This can be done in a variety of ways but I have chosen to do this by building a machine learning classifier which was trained on an immune cell atlas from mice with nonalcoholic fatty liver disease (NAFLD) which matches our data (Guilliams et al. 2022). After validating our annotations by looking at marker genes we can begin our analysis. Single cell is a powerful modality and can allow for various directions of analysis so after doing some basic analysis and we will build off of the trends and patterns we observe. 


# Questions: 

After completing the preprocessing of the data we can dive into exploring and beginning the analysis. However, before we even begin this process we have already seen a few interesting patterns we have seen in our preliminary, unpublished data so far. The first question we have is why does the proportion of CD8+ T cells increase with HCC severity (HCC scores)? The current hypothesis is that this is an immune response to the tumor and that these cells are getting exhausted leading to an influx of more T cells. The second question is why does the proportion of myeloid cells decrease with HCC severity and what kinds of myeloid cells are these? Since there are many different myeloid cells in the data including native Kupffer cells this is important to understand. Additionally we would like to know: what are the cellular processes involved with HCC? Lastly we have another question we would like to ask to guide this analysis: what role does T cell exhaustion play in this phenomenon? 

# Proposed Analysis: 

From a computational perspective there are a variety of tools we must employ to answer these questions. Single cell data is in matrix format and there are many tools in python and R such as Scanpy and Seurat respectively that are convenient for single cell data analysis (Wolf et al. 2018, Hao et al 2024). As I am more familiar with python I will do much of the data visualization to generate figures using matplotlib and seaborn (Hunter 2007; Waskom 2021). Some of the analyses we will conduct and visualize include running gene set enrichment analysis (GSEA) to visualize the cellular processes with the largest effect that are statistically significant (Hanzelmann et al. 2013). Some other visualizations we can generate include bar plots to show the cell type proportions in each of the conditions: TW, TN, WW, and WN. Although these are some of the main visualizations and questions we want to answer, we will use these to guide future visualization and additional questions are sure to arise as we continue researching. Another important aspect of these visualization scripts and data processing is to make them scalable and reproducible since we will need to run these scripts later on additional data. This project focuses on male 6 month mice but our collaborators are working on generating additional data so it is important to make our scripts modular and document our work rigorously.


# Conclusion and Feasability

Overall the goal of this project is to investigate some of the cellular processes and patterns that relate and contribute to hepatocellular carcinoma in mice with CHIP. By combining this analysis with wet lab analysis from our collaborators we hope to better understand this phenomenon. Although there are many open questions related to this topic, I feel comfortable pursuing this project due to my familiarity with single cell and my access to external help from my PI and collaborators. There is much work left to be done but I have already preprocessed some of the data and done some of the basic analysis and believe I can generate some novel insights by the end of this class and present them in the final report. I also am familiar with the tools and programming and will document this aspect as well. 


# References:

Guilliams M, Bonnardel J, Haest B, Vanderborght B, Wagner C, 
Remmerie A, Bujko A, Martens L, Thoné T, Browaeys R, et al. 
2022. Spatial proteogenomics reveals distinct and evolutionarily 
conserved hepatic macrophage niches. Cell 185: 379–396.

Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, 
Srivastava A, Molla G, Madad S, Fernandez-Granda C, et al. 2024. 
Dictionary learning for integrative, multimodal and scalable 
single-cell analysis. Nat Biotechnol 42: 293–304.

Hänzelmann S, Castelo R, Guinney J. 2013. GSVA: gene set variation
analysis for microarray and RNA-seq data. BMC Bioinformatics 14: 7.

Hunter JD. 2007. Matplotlib: a 2D graphics environment.
Comput Sci Eng 9: 90–95.

Marnell CS, Bick A, Natarajan P. 2022. Clonal hematopoiesis of
indeterminate potential (CHIP): linking somatic mutations,
hematopoiesis, chronic inflammation and cardiovascular disease.
J Mol Cell Cardiol 161: 98–105.

Mitchell E, Spencer Chapman M, Williams N, et al. 2022. Clonal
dynamics of haematopoiesis across the human lifespan. Nature 606:
343–350.

Sikking MA, Stroeks SLVM, Waring OJ, Henkens MTHM, Riksen NP,
Hoischen A, Heymans SRB, Verdonschot JAJ. 2023. Clonal hematopoiesis
of indeterminate potential from a heart failure specialist's point of
view. J Am Heart Assoc 12: e030603.

Waskom ML. 2021. Seaborn: statistical data visualization.
J Open Source Softw 6: 3021.

Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell
gene expression data analysis. Genome Biol 19: 15.

Wong WJ, Emdin C, Bick AG, Zekavat SM, Niroula A, Pirruccello JP,
Dichtel L, Griffin G, Uddin MM, Gibson CJ, Kovalcik V, Lin AE,
McConkey ME, Vromman A, Sellar RS, Kim PG, Agrawal M, Weinstock J,
Long MT, Yu B, et al. 2023. Clonal haematopoiesis and risk of chronic
liver disease. Nature 616: 747–754.

