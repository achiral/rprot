#DEP
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/DEP") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
#install.packages("VIM")
#install.packages("imputeMissings")
#install.packages("mice")
#install.packages("mlbench")
#install.packages("missForest")
#install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DEP")
#BiocManager::install("BaylorEdPsych") #インストールできない
#############################################################
#ライブラリ読み込み
#library(VIM)
#library(BaylorEdPsych) #インストールできない
#library(imputeMissings)
#library(mice)
#library(mlbench)
#library(missForest)
#library(SummarizedExperiment)
library(DEP)
#library(devtools)
#############################################################
#ライブラリセット
#library(EnhancedVolcano)
#library(magrittr)
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
#library(dplyr)
#library(scales) #muted()関数使用のため
#library(rJava)
#library(genefilter) #ヒートマップ
#library(gplots) #ヒートマップ
#library(ComplexHeatmap) #ヒートマップ
#library(RColorBrewer) #色
library(readxl) #エクセル入力(read_excel)
#library(openxlsx) #エクセル入出力(write.xlsx)
#library(gridExtra) #svg出力のため
#library(cowplot)
#############################################################
#Package DEP
#https://www.bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html#installation
#library(DEP)

#Interactive analysis using the DEP Shiny apps
# For LFQ analysis
#run_app("LFQ")
# For TMT analysis
#run_app("TMT")

#Differential analysis
# Loading a package required for data handling
#library("dplyr")

# The data is provided with the package
#dataa <- UbiLength
# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively. 
#data <- filter(data, Reverse != "+", Potential.contaminant != "+") #フィルタ
#dim(data) #The data.frame dimensions:
#colnames(data) #The data.frame column names:
#The “LFQ.intensity” columns will be used for subsequent analysis.
#############################################################
rm(list = ls(all.names = TRUE))
#xlsx入力
data <- read_excel("AMY.xlsx", 3) #シート3入力
ExpDesign <- read_excel("AMY.xlsx", 2) #シート2入力
dim(data) #The data.frame dimensions:
colnames(data) #The data.frame column names:
#The "(SAL|PCP)" columns will be used for subsequent analysis.
#Data preparation
###############################################################################
#デモデータではlog2-transformedデータになっているため、ここで変換しておく？→#make_se関数で処理されるため不要
#log2-transformation
#data2 <- data
#Sample_columns <- grep("(SAL|PCP)", colnames(data2)) # get Sample column numbers
#Sample_columns
#data2[, Sample_columns] <- log(data2[Sample_columns], 2)

#xlsx出力
#smp <- list("data1"=data, "ExpDesign"=ExpDesign, "data2"=data2) #リスト作成
#write.xlsx(smp, "AMY2.xlsx") #シート出力
#############################################################
#rm(list = ls(all.names = TRUE))
#xlsx入力
#data <- read_excel("AMY2.xlsx", 3) #シート3入力
#ExpDesign <- read_excel("AMY2.xlsx", 2) #シート2入力
#dim(data) #The data.frame dimensions:
#colnames(data) #The data.frame column names:
#The "(SAL|PCP)" columns will be used for subsequent analysis.
#Data preparation
###############################################################################
#置換
#split <- str_replace(data$`Peak Name`, pattern="\\|", replacement="\\_")
#分割
#split <- str_split(data$`Peak Name`, pattern = c("\\|","\\_"), simplify = TRUE) #"_"をうまく認識しない
split <- str_split(data$`Peak Name`, pattern = "\\|", simplify = TRUE)
colnames(split) <- c("sp", "Protein.IDs", "GeneName") #列名変更
class(split)
x <- data.frame(split)
#y <- select(.data = x, GeneName)
#文字抽出
Protein.IDs <- str_sub(x$`Protein.IDs`, start = 1, end = 6) #`Peak Name`列の1-6文字目(Protein.IDs)抽出
Gene.names <- str_sub(x$`GeneName`, start = 1, end = -7) #`GeneName`列の1文字目〜-7文字目(GeneName)抽出
Species <- str_sub(x$`GeneName`, start = -5, end = -1) #`GeneName`列の-5〜-1文字目(Species)抽出
#split <- str_split(data$`Peak Name`, pattern = "\\_", simplify = TRUE) #"_"をうまく認識しない
#文字抽出
#Protein.IDs <- str_sub(data$`Peak Name`, start = 4, end = 9) #`Peak Name`列の4-9文字目(Protein.IDs)抽出
#Gene.names <- str_sub(data$`Peak Name`, start = 11, end = -7) #`Peak Name`列の11文字目〜-7文字目(Gene.names)抽出
data <- cbind(data, Protein.IDs) #dataとProtein.IDsを列ベクトル単位で結合 
data <- cbind(data, Gene.names) #dataとGene.namesを列ベクトル単位で結合
data <- cbind(data, Species) #dataとSpeciesを列ベクトル単位で結合

# Are there any duplicated Protein IDs/gene names/Species?
data$Protein.IDs %>% duplicated() %>% any()
data$Gene.names %>% duplicated() %>% any()
data$Species %>% duplicated() %>% any()

# Make a table of duplicated Protein IDs/gene names/Species
data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data %>% group_by(Species) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

#For further analysis these proteins must get unique names. Additionally, some proteins do not have an annotated gene name and for those we will use the Uniprot ID.
# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
# Are there any duplicated names?
data$name %>% duplicated() %>% any()

#Generate a SummarizedExperiment object
# Generate a SummarizedExperiment object using an experimental design
Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
experimental_design <- ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
data_se <- make_se(data_unique, Sample_columns, experimental_design) #columns=データ数
# Generate a SummarizedExperiment object by parsing condition information from the column names
#Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
#data_se_parsed <- make_se_parse(data_unique, Sample_columns) #SWATHはこっちでよいかも
# Let's have a look at the SummarizedExperiment object
data_se
#data_se_parsed

###############################################################################
#ここから---------------------------------------------
#Prerequisites of the SummarizedExperiment object
#The make_se and make_se_parse functions generate a SummarizedExperiment object that has a couple of specifications. The assay data is log2-transformed and its rownames depict the protein names. The rowData contains, amongst others, the ‘name’ and ‘ID’ columns that were generated by make_unique. The colData contains the experimental design and thereby the sample annotation. Thereby the colData includes the ‘label’, ‘condition’ and ‘replicate’ columns as well as a newly generated ‘ID’ column. The log2-transformed assay data and the specified rowData and colData columns are prerequisites for the subsequent analysis steps.
#デモデータではlog2-transformedデータになっているため、ここで変換しておく？
#log2-transformation
#y <- data.frame(data_se)
#Sample_columns <- grep("(SAL|PCP)", colnames(data)) # get Sample column numbers
#Sample_columns
#data[, Sample_columns] <- log(data[Sample_columns], 2)
###############################################################################



###############################################################################
#Filter on missing values
#The dataset contains proteins which are not quantified in all replicates. Some proteins are even only quantified in a single replicate.
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)
#plot_frequency(data_se_parsed)
#This leaves our dataset with missing values, which need to be imputed. However, this should not be done for proteins that contain too many missing values. Therefore, we first filter out proteins that contain too many missing values. This is done by setting the threshold for the allowed number of missing values per condition in the filter_missval function.
# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)
# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
#data_filt2 <- filter_missval(data_se, thr = 1)
# Filter for proteins that are identified in 1 out of 3 replicates of at least one condition
#data_filt3 <- filter_missval(data_se, thr = 2)
# Filter for proteins that are identified in none out of 3 replicates of at least one condition
#data_filt4 <- filter_missval(data_se, thr = 3)
#After filtering, the number of identified proteins per sample can be plotted as well as the overlap in identifications between samples.
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_se)
#plot_numbers(data_se_parsed)
#plot_numbers(data_filt)
#plot_numbers(data_filt2)
#plot_numbers(data_filt3)
#plot_numbers(data_filt4)
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_se)
#plot_coverage(data_filt)
#plot_coverage(data_filt2)
#plot_coverage(data_filt3)
#plot_coverage(data_filt4)
#Normalization
#The data is background corrected and normalized by variance stabilizing transformation (vsn).
# Normalize the data
#data_norm <- normalize_vsn(data_se)
#data_norm <- normalize_vsn(data_filt)
#The normalization can be inspected by checking the distributions of the samples before and after normalization.
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_se) #3つ以上も表示可能、ggplot2+で編集可能
#plot_normalization(data_se, data_norm) #3つ以上も表示可能、ggplot2+で編集可能
#plot_normalization(data_filt, data_norm) #3つ以上も表示可能、ggplot2+で編集可能

###############################################################################
#Impute data for missing values
#The remaining missing values in the dataset need to be imputed. The data can be missing at random (MAR), for example if proteins are quantified in some replicates but not in others. Data can also be missing not at random (MNAR), for example if proteins are not quantified in specific conditions (e.g. in the control samples). MNAR can indicate that proteins are below the detection limit in specific samples, which could be very well the case in proteomics experiments. For these different conditions, different imputation methods have to be used, as described in the MSnbase vignette and more specifically in the impute function descriptions.
#To explore the pattern of missing values in the data, a heatmap is plotted indicating whether values are missing (0) or not (1). Only proteins with at least one missing value are visualized.
# Plot a heatmap of proteins with missing values
#plot_missval(data_filt)
#plot_missval(data_norm)
#This heatmap indicates that missing values are highly biased to specific samples. The example dataset is an affinity enrichment dataset of ubiquitin interactors, which is likely to have proteins which are below the detection limit in specific samples. These can be proteins that are specifically enriched in the ubiquitin purifications, but are not enriched in the controls samples, or vice versa. To check whether missing values are biased to lower intense proteins, the densities and cumulative fractions are plotted for proteins with and without missing values.
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
#plot_detect(data_filt)
#plot_detect(data_norm)
#Indeed the proteins with missing values have on average low intensities. This data (MNAR and close to the detection limit) should be imputed by a left-censored imputation method, such as the quantile regression-based left-censored function (“QRILC”) or random draws from a left-shifted distribution (“MinProb” and “man”). In contrast, MAR data should be imputed with methods such as k-nearest neighbor (“knn”) or maximum likelihood (“MLE”) functions. See the MSnbase vignette and more specifically the impute function description for more information.
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_se, fun = "")
#impute(data_norm, fun = "")
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
#data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man_se <- impute(data_se, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
#data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
# Impute missing data using the k-nearest neighbour approach (for MAR)
#data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
#The effect of the imputation on the distributions can be visualized.
# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp)
plot_imputation(data_se, data_imp_man_se) #3つ以上も表示可能、ggplot2+で編集可能
plot_imputation(data_se, data_norm, data_imp_man_se, data_imp_man) #3つ以上も表示可能、ggplot2+で編集可能
#plot_imputation(data_norm, data_imp_knn)
plot_normalization(data_se, data_imp_man_se) #3つ以上も表示可能、ggplot2+で編集可能
###############################################################################
#Differential enrichment analysis
#Protein-wise linear models combined with empirical Bayes statistics are used for the differential enrichment analysis (or differential expression analysis). The test_diff function introduced here uses limma and automatically generates the contrasts to be tested. For the contrasts generation, the control sample has to be specified. Additionally, the types of contrasts to be produced need to be indicated, allowing the generation of all possible comparisons (“all”) or the generation of contrasts of every sample versus control (“control”). Alternatively, the user can manually specify the contrasts to be tested (type = “manual”), which need to be specified in the argument test.
# Differential enrichment analysis  based on linear models and empherical Bayes statistics
# Test every sample versus control
data_diff <- test_diff(data_imp_man, type = "control", control = "SC0")
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp_man, type = "all")
# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp_man, type = "manual", test = c("SC10_vs_SC0","SC30_vs_SC0","PC0_vs_SC0"))
#Finally, significant proteins are defined by user-defined cutoffs using add_rejections.
# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
dep2 <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5))
dep3 <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
#The results from the previous analysis can be easily visualized by a number of functions. These visualizations assist in the determination of the optimal cutoffs to be used, highlight the most interesting samples and contrasts, and pinpoint differentially enriched/expressed proteins.
#PCA plot
#The PCA plot can be used to get a high-level overview of the data. This can be very useful to observe batch effects, such as clear differences between replicates.
# Plot the first and second principal components
plot_pca(dep2, x = 1, y = 2, n = 500, point_size = 2)

#Correlation matrix
#A correlation matrix can be plotted as a heatmap, to visualize the Pearson correlations between the different samples.
# Plot the Pearson correlation matrix
plot_cor(dep2, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#Heatmap of all significant proteins
#The heatmap representation gives an overview of all significant proteins (rows) in all samples (columns). This allows to see general trends, for example if one sample or replicate is really different compared to the others. Additionally, the clustering of samples (columns) can indicate closer related samples and clustering of proteins (rows) indicates similarly behaving proteins. The proteins can be clustered by k-means clustering (kmeans argument) and the number of clusters can be defined by argument k.
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep2, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
#The heatmap shows a clustering of replicates and indicates that 4Ubi and 6Ubi enrich a similar repertoire of proteins. The k-means clustering of proteins (general clusters of rows) nicely separates protein classes with different binding behaviors.
#Alternatively, a heatmap can be plotted using the contrasts, i.e. the direct sample comparisons, as columns. Here, this emphasises the enrichment of ubiquitin interactors compared to the control sample.
# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep2, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

#Volcano plots of specific contrasts
#Volcano plots can be used to visualize a specific contrast (comparison between two samples). This allows to inspect the enrichment of proteins between the two samples (x axis) and their corresponding adjusted p value (y axis). The add_names argument can be set to FALSE if the protein labels should be omitted, for example if there are too many names.
# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "PC0_vs_SC0", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "PC10_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "PC30_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "SC10_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "SC30_vs_SC0", label_size = 2, add_names = TRUE)

#Barplots of a protein of interest
#It can also be useful to plot the data of a single protein, for example if this protein is of special interest.
# Plot a barplot for SYN1 and ANK3
plot_single(dep, proteins = c("SYN1", "ANK3"))
# Plot a barplot for the protein SYN1 with the data centered
plot_single(dep, proteins = "SYN1", type = "centered")

#Frequency plot of significant proteins and overlap of conditions
#Proteins can be differentially enriched/expressed in multiple comparisons. To visualize the distribution of significant conditions per protein and the overlap between conditions, the plot_cond function can be used.
# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)
plot_cond(dep2)
plot_cond(dep3)

#Results table
#To extract a table containing the essential results, the get_results function can be used.
# Generate a results table
data_results <- get_results(dep2)
# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
#The resulting table contains the following columns:
# Column names of the results table
colnames(data_results)
#Of these columns, the p.val and p.adj columns contain the raw and adjusted p values, respectively, for the contrast as depicted in the column name. The ratio columns contain the average log2 fold changes. The significant columns indicate whether the protein is differentially enriched/expressed, as defined by the chosen cutoffs. The centered columns contain the average log2 fold changes scaled by protein-wise centering.

#Generate a data.frame from the resulting SummarizedExperiment object
#You might want to obtain an ordinary data.frame of the results. For this purpose, the package provides functions to convert SummarizedExperiment objects to data.frames. get_df_wide will generate a wide table, whereas get_df_long will generate a long table.
# Generate a wide data.frame
df_wide <- get_df_wide(dep2)
# Generate a long data.frame
df_long <- get_df_long(dep2)

#Save your results object for reuse
#To facilitate future analysis and/or visualization of your current data, saving your analyzed data is highly recommended. We save the final data object (dep) as well as intermediates of the analysis, i.e. the initial SummarizedExperiment object (data_se), normalized data (data_norm), imputed data (data_imp) and differentially expression analyzed data (data_diff). This allows us to easily change parameters in future analysis.
# Save analyzed data
save(data_se, data_norm, data_imp_man, data_diff_all_contrasts, dep2, file = "data.RData")
#xlsx出力
smp <- list("df_wide"=df_wide,"df_long"=df_long) #リスト作成
write.xlsx(smp, "AMY3.xlsx") #シート出力
# These data can be loaded in future R sessions using this command
load("data.RData")
#############################################################