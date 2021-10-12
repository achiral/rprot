#Imputation
#欠損値データ
#https://techblog.nhn-techorus.com/archives/6309
#https://www.bioconductor.org/packages/release/bioc/html/DEP.html
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/Imputation") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
#install.packages("VIM")
#install.packages("imputeMissings")
#install.packages("mice")
#install.packages("mlbench")
#install.packages("missForest")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DEP")
#BiocManager::install("BaylorEdPsych") #インストールできない
#############################################################
#ライブラリ読み込み
library(VIM)
#library(BaylorEdPsych) #インストールできない
library(imputeMissings)
library(mice)
library(mlbench)
library(missForest)
library("SummarizedExperiment")
library("DEP")
#############################################################
#ライブラリセット
library(EnhancedVolcano)
library(magrittr)
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(dplyr)
library(scales) #muted()関数使用のため
library(scales) #muted()関数使用のため
library(rJava)
library(genefilter) #ヒートマップ
library(gplots) #ヒートマップ
library(ComplexHeatmap) #ヒートマップ
library(RColorBrewer) #色
library(readxl) #エクセル入力(read_excel)
library(openxlsx) #エクセル入出力(write.xlsx)
library(gridExtra) #svg出力のため
library(cowplot)
#############################################################
data(BostonHousing, package = "mlbench") #BostonHousing データセット読み込み
data <- BostonHousing
summary(data) #完全データの要約表示

set.seed(123)
data.mis <- missForest::prodNA(data, noNA = 0.2) #欠測値の割合を20%に設定
summary(data.mis)

vim.aggr <- aggr(data.mis,
                 col = c('white','red'),
                 numbers = TRUE,
                 sortVars = FALSE,
                 prop = TRUE,
                 labels = names(data.mis),
                 cex.axis = .8,
                 gap = 3) #欠測値の割合や欠測パターンの組み合わせを確認
vim.aggr

#library(BaylorEdPsych) #インストールできない
#test.mcar <- LittleMCAR(data.mis) #LittleのMCAR検定#できない
#print(test.mcar$p.value)

#model <- lm(medv ~ ., data = data)
#summary(model)
#model.mis <- lm(medv ~ ., data = data.mis, na.action = na.omit) #エラー
#summary(model.mis) #エラー

#中央値代入
library(imputeMissings)
data.comp.median <- impute(data.mis, method = "median/mode") #中央値を代入
model.median <- lm(medv ~ ., data = data.comp.median) #線形回帰
summary(model.median)

#多重代入法 (Multivariate Imputation by Chained Equations; MICE)
##mを10、methodを予測平均マッチング (predictive mean matching; pmm)：回帰代入とマッチングを組み合わせた手法で回帰した値に近い観測値をランダムに選択し代入
library(mice)
imp <- mice(data.mis,
            m = 10,
            maxit = 50,
            method = "pmm",
            printFlag = FALSE,
            seed = 12345)
summary(imp)
densityplot(imp) #青線が欠損データ、 赤線が m 個の擬似完全データの分布

#m 個の擬似完全データからパラメータを推定し統合
library(magrittr)
combine <- imp %>%
  with(lm(medv ~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + b + lstat)) %>%
  pool()
summary(combine)
#############################################################
#https://www.bioconductor.org/packages/release/bioc/html/DEP.html
# Variables
replicates = 3
bg_proteins = 3000
DE_proteins = 300
log2_mean_bg = 27
log2_sd_bg = 2
log2_mean_DE_control = 25
log2_mean_DE_treatment = 30
log2_sd_DE = 2

# Loading DEP and packages required for data handling
library("DEP")
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")

# Background data (sampled from null distribution)
sim_null <- data_frame(
  name = paste0("bg_", rep(1:bg_proteins, rep((2*replicates), bg_proteins))),
  ID = rep(1:bg_proteins, rep((2*replicates), bg_proteins)),
  var = rep(c("control_1", "control_2", "control_3", 
              "treatment_1", "treatment_2", "treatment_3"), bg_proteins), 
  val = 2^rnorm(2*replicates*bg_proteins, mean = log2_mean_bg, sd = log2_sd_bg))

# Data for DE proteins (sampled from alternative distribution)
sim_diff <- rbind(
  data_frame(
    name = paste0("DE_", rep(1:DE_proteins, rep(replicates, DE_proteins))),
    ID = rep((bg_proteins+1):(bg_proteins+DE_proteins), 
             rep(replicates, DE_proteins)),
    var = rep(c("control_1", "control_2", "control_3"), DE_proteins), 
    val = 2^rnorm(replicates*DE_proteins, 
                  mean = log2_mean_DE_control, sd = log2_sd_DE)),
  data_frame(
    name = paste0("DE_", rep(1:DE_proteins, rep(replicates, DE_proteins))),
    ID = rep((bg_proteins+1):(bg_proteins+DE_proteins), 
             rep(replicates, DE_proteins)),
    var = rep(c("treatment_1", "treatment_2", "treatment_3"), DE_proteins),
    val = 2^rnorm(replicates*DE_proteins, 
                  mean = log2_mean_DE_treatment, sd = log2_sd_DE)))

# Combine null and DE data
sim <- rbind(sim_null, sim_diff) %>% 
  spread(var, val) %>% 
  arrange(ID)

# Generate experimental design
experimental_design <- data_frame(
  label = colnames(sim)[!colnames(sim) %in% c("name", "ID")],
  condition = c(rep("control", replicates), 
                rep("treatment", replicates)),
  replicate = rep(1:replicates, 2))

#Introduce missing values
# Variables
MAR_fraction = 0.05
MNAR_proteins = 100

# Generate a MAR matrix
MAR_matrix <- matrix(
  data = sample(c(TRUE, FALSE), 
                size = 2*replicates*(bg_proteins+DE_proteins), 
                replace = TRUE, 
                prob = c(MAR_fraction, 1-MAR_fraction)), 
  nrow = bg_proteins+DE_proteins, 
  ncol = 2*replicates)

# Introduce missing values at random (MAR)
controls <- grep("control", colnames(sim))
treatments <- grep("treatment", colnames(sim))
sim[, c(controls, treatments)][MAR_matrix] <- 0
sim$MAR <- apply(MAR_matrix, 1, any)

# Introduce missing values not at random (MNAR)
DE_protein_IDs <- grep("DE", sim$name)
sim[DE_protein_IDs[1:MNAR_proteins], controls] <- 0
sim$MNAR <- FALSE
sim$MNAR[DE_protein_IDs[1:MNAR_proteins]] <- TRUE

#Generate a SummarizedExperiment
# Generate a SummarizedExperiment object
sim_unique_names <- make_unique(sim, "name", "ID", delim = ";")
se <- make_se(sim_unique_names, c(controls, treatments), experimental_design)


# Plot a barplot of the protein quantification overlap between samples
plot_frequency(se)

# No filtering
no_filter <- se

# Filter for proteins that are quantified in all replicates of at least one condition
condition_filter <- filter_proteins(se, "condition", thr = 0)

# Filter for proteins that have no missing values
complete_cases <- filter_proteins(se, "complete")

# Filter for proteins that are quantified in at least 2/3 of the samples.
frac_filtered <- filter_proteins(se, "fraction", min = 0.66)

# Function to extract number of proteins
number_prots <- function(se) {
  names <- rownames(get(se))
  data_frame(Dataset = se,
             bg_proteins = sum(grepl("bg", names)),
             DE_proteins = sum(grepl("DE", names)))
}

# Number of bg and DE proteins still included
objects <- c("no_filter", 
             "condition_filter",
             "complete_cases",
             "frac_filtered")

map_df(objects, number_prots)


# Scale and variance stabilize
no_filter <- normalize_vsn(se)
condition_filter <- normalize_vsn(condition_filter)
complete_cases <- normalize_vsn(complete_cases)
frac_filtered <- normalize_vsn(frac_filtered)

# Mean versus Sd plot
meanSdPlot(no_filter)



#Data imputation of missing data
#Explore the pattern of missing values
# Plot a heatmap of proteins with missing values
plot_missval(no_filter)
# Plot intensity distributions and cumulative fraction of proteins 
# with and without missing values
plot_detect(no_filter)

#Imputation options
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(no_filter, fun = "")
# No imputation
no_imputation <- no_filter
# Impute missing data using random draws from a 
# Gaussian distribution centered around a minimal value (for MNAR)
MinProb_imputation <- impute(no_filter, fun = "MinProb", q = 0.01)
# Impute missing data using random draws from a 
# manually defined left-shifted Gaussian distribution (for MNAR)
manual_imputation <- impute(no_filter, fun = "man", shift = 1.8, scale = 0.3)
# Impute missing data using the k-nearest neighbour approach (for MAR)
knn_imputation <- impute(no_filter, fun = "knn", rowmax = 0.9)
# Plot intensity distributions before and after imputation
plot_imputation(no_filter, MinProb_imputation, 
                manual_imputation, knn_imputation)

#Advanced imputation methods
#Mixed imputation on proteins (rows)
# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(no_filter) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()
# Get a logical vector
MNAR <- names(no_filter) %in% proteins_MNAR
# Perform a mixed imputation
mixed_imputation <- impute(
  no_filter, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "zero") # imputation function for MNAR
#Mixed imputation on samples (columns)
# SummarizedExperiment to MSnSet object conversion
sample_specific_imputation <- no_filter
MSnSet <- as(sample_specific_imputation, "MSnSet")
# Impute differently for two sets of samples
MSnSet_imputed1 <- MSnbase::impute(MSnSet[, 1:3], method = "MinProb")
MSnSet_imputed2 <- MSnbase::impute(MSnSet[, 4:6], method = "knn")
# Combine into the SummarizedExperiment object
assay(sample_specific_imputation) <- cbind(
  MSnbase::exprs(MSnSet_imputed1), 
  MSnbase::exprs(MSnSet_imputed2))
# Plot intensity distributions before and after imputation
plot_imputation(no_filter, mixed_imputation, sample_specific_imputation)

#Test for differential expression
# Function that wraps around test_diff, add_rejections and get_results functions
DE_analysis <- function(se) {
  se %>% 
    test_diff(., type = "control", control = "control") %>%
    add_rejections(., alpha = 0.1, lfc = 0) %>% 
    get_results()
}
# DE analysis on no, knn, MinProb and mixed imputation
no_imputation_results <- DE_analysis(no_imputation)
knn_imputation_results <- DE_analysis(knn_imputation)
MinProb_imputation_results <- DE_analysis(MinProb_imputation)
mixed_imputation_results <- DE_analysis(mixed_imputation)

#Number of identified differentially expressed proteins
# Function to extract number of DE proteins
DE_prots <- function(results) {
  data_frame(Dataset = gsub("_results", "", results),
             significant_proteins = get(results) %>% 
               filter(significant) %>% 
               nrow())
}
# Number of significant proteins
objects <- c("no_imputation_results", 
             "knn_imputation_results",
             "MinProb_imputation_results",
             "mixed_imputation_results")
map_df(objects, DE_prots)

#ROC curves
# Function to obtain ROC data
get_ROC_df <- function(results) {
  get(results) %>% 
    select(name, treatment_vs_control_p.val, significant) %>% 
    mutate(
      DE = grepl("DE", name),
      BG = grepl("bg", name)) %>% 
    arrange(treatment_vs_control_p.val) %>% 
    mutate(
      TPR = cumsum(as.numeric(DE)) / 300,
      FPR = cumsum(as.numeric(BG)) / 3000,
      method = results)
}
# Get ROC data for no, knn, MinProb and mixed imputation
ROC_df <- map_df(objects, get_ROC_df)
# Plot ROC curves
ggplot(ROC_df, aes(FPR, TPR, col = method)) +
  geom_line() +
  theme_DEP1() +
  ggtitle("ROC-curve")
# Plot ROC curves zoom
ggplot(ROC_df, aes(FPR, TPR, col = method)) +
  geom_line() +
  theme_DEP1() +
  xlim(0, 0.1) +
  ggtitle("ROC-curve zoom")

#Differences in response
# Function to obtain summary data
get_rejected_proteins <- function(results) {
  get(results) %>% 
    filter(significant) %>% 
    left_join(., select(sim, name, MAR, MNAR), by = "name") %>% 
    mutate(
      DE = grepl("DE", name),
      BG = grepl("bg", name),
      method = results)
}

# Get summary data for no, knn, MinProb and mixed imputation
objects <- c("no_imputation_results", 
             "knn_imputation_results",
             "MinProb_imputation_results",
             "mixed_imputation_results")

summary_df <- map_df(objects, get_rejected_proteins)

# Plot number of DE proteins (True and False)
summary_df %>% 
  group_by(method) %>% 
  summarize(TP = sum(DE), FP = sum(BG)) %>% 
  gather(category, number, -method) %>% 
  mutate(method = gsub("_results", "", method)) %>% 
  ggplot(aes(method, number, fill = category)) +
  geom_col(position = position_dodge()) +
  theme_DEP2() +
  labs(title = "True and False Hits",
       x = "",
       y = "Number of DE proteins",
       fill = "False or True")

# Plot number of DE proteins with missing values
summary_df %>% 
  group_by(method) %>% 
  summarize(MNAR = sum(MNAR), MAR = sum(MAR)) %>% 
  gather(category, number, -method) %>% 
  mutate(method = gsub("_results", "", method)) %>% 
  ggplot(aes(method, number, fill = category)) +
  geom_col(position = position_dodge()) +
  theme_DEP2() +
  labs(title = "Category of Missing Values",
       x = "",
       y = "Number of DE proteins",
       fill = "")

data <- UbiLength
summary(UbiLength)
data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#############################################################
#Package DEP
#https://www.bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html#installation
library("DEP")

#Interactive analysis using the DEP Shiny apps
# For LFQ analysis
run_app("LFQ")

# For TMT analysis
run_app("TMT")

#Differential analysis
# Loading a package required for data handling
library("dplyr")

#############################################################
# The data is provided with the package
data <- UbiLength
# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively. 
data <- filter(data, Reverse != "+", Potential.contaminant != "+") #フィルタ
dim(data) #The data.frame dimensions:
colnames(data) #The data.frame column names:
#The “LFQ.intensity” columns will be used for subsequent analysis.

#Data preparation
# Are there any duplicated gene names?
data$Gene.names %>% duplicated() %>% any()
# Make a table of duplicated gene names
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
#For further analysis these proteins must get unique names. Additionally, some proteins do not have an annotated gene name and for those we will use the Uniprot ID.
# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
# Are there any duplicated names?
data$name %>% duplicated() %>% any()

#############################################################
#Generate a SummarizedExperiment object
# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- UbiLength_ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
data_se <- make_se(data_unique, LFQ_columns, experimental_design) #coumns=データ数(SWATHであればサンプル名になっているのでカラム名の変更が必要)
# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_unique, LFQ_columns) #SWATHはこっちでよいかも
# Let's have a look at the SummarizedExperiment object
data_se
data_se_parsed

#xlsx出力
smp <- list("UbiLength"=data, "UbiLength_ExpDesign"=experimental_design) #リスト作成
write.xlsx(smp, "UbiLength.xlsx") #シート出力

#############################################################
#Prerequisites of the SummarizedExperiment object
#The make_se and make_se_parse functions generate a SummarizedExperiment object that has a couple of specifications. The assay data is log2-transformed and its rownames depict the protein names. The rowData contains, amongst others, the ‘name’ and ‘ID’ columns that were generated by make_unique. The colData contains the experimental design and thereby the sample annotation. Thereby the colData includes the ‘label’, ‘condition’ and ‘replicate’ columns as well as a newly generated ‘ID’ column. The log2-transformed assay data and the specified rowData and colData columns are prerequisites for the subsequent analysis steps.
#デモデータではlog2-transformedデータになっているため、ここで変換しておく？

#Filter on missing values
#The dataset contains proteins which are not quantified in all replicates. Some proteins are even only quantified in a single replicate.
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)
plot_frequency(data_se_parsed)
#This leaves our dataset with missing values, which need to be imputed. However, this should not be done for proteins that contain too many missing values. Therefore, we first filter out proteins that contain too many missing values. This is done by setting the threshold for the allowed number of missing values per condition in the filter_missval function.
# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)
# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)
# Filter for proteins that are identified in 1 out of 3 replicates of at least one condition
data_filt3 <- filter_missval(data_se, thr = 2)
# Filter for proteins that are identified in none out of 3 replicates of at least one condition
data_filt4 <- filter_missval(data_se, thr = 3)
#After filtering, the number of identified proteins per sample can be plotted as well as the overlap in identifications between samples.
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_se)
plot_numbers(data_filt)
plot_numbers(data_filt2)
plot_numbers(data_filt3)
plot_numbers(data_filt4)
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_se)
plot_coverage(data_filt)
plot_coverage(data_filt2)
plot_coverage(data_filt3)
plot_coverage(data_filt4)
#Normalization
#The data is background corrected and normalized by variance stabilizing transformation (vsn).
# Normalize the data
data_norm <- normalize_vsn(data_filt)
#The normalization can be inspected by checking the distributions of the samples before and after normalization.
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) #3つ以上も表示可能、ggplot2+で編集可能

#Impute data for missing values
#The remaining missing values in the dataset need to be imputed. The data can be missing at random (MAR), for example if proteins are quantified in some replicates but not in others. Data can also be missing not at random (MNAR), for example if proteins are not quantified in specific conditions (e.g. in the control samples). MNAR can indicate that proteins are below the detection limit in specific samples, which could be very well the case in proteomics experiments. For these different conditions, different imputation methods have to be used, as described in the MSnbase vignette and more specifically in the impute function descriptions.
#To explore the pattern of missing values in the data, a heatmap is plotted indicating whether values are missing (0) or not (1). Only proteins with at least one missing value are visualized.
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)
#This heatmap indicates that missing values are highly biased to specific samples. The example dataset is an affinity enrichment dataset of ubiquitin interactors, which is likely to have proteins which are below the detection limit in specific samples. These can be proteins that are specifically enriched in the ubiquitin purifications, but are not enriched in the controls samples, or vice versa. To check whether missing values are biased to lower intense proteins, the densities and cumulative fractions are plotted for proteins with and without missing values.
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)
#Indeed the proteins with missing values have on average low intensities. This data (MNAR and close to the detection limit) should be imputed by a left-censored imputation method, such as the quantile regression-based left-censored function (“QRILC”) or random draws from a left-shifted distribution (“MinProb” and “man”). In contrast, MAR data should be imputed with methods such as k-nearest neighbor (“knn”) or maximum likelihood (“MLE”) functions. See the MSnbase vignette and more specifically the impute function description for more information.
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
#The effect of the imputation on the distributions can be visualized.
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)
plot_imputation(data_norm, data_imp_man)
plot_imputation(data_norm, data_imp_knn)

#Differential enrichment analysis
#Protein-wise linear models combined with empirical Bayes statistics are used for the differential enrichment analysis (or differential expression analysis). The test_diff function introduced here uses limma and automatically generates the contrasts to be tested. For the contrasts generation, the control sample has to be specified. Additionally, the types of contrasts to be produced need to be indicated, allowing the generation of all possible comparisons (“all”) or the generation of contrasts of every sample versus control (“control”). Alternatively, the user can manually specify the contrasts to be tested (type = “manual”), which need to be specified in the argument test.
# Differential enrichment analysis  based on linear models and empherical Bayes statistics
# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")
# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))
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
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

#Correlation matrix
#A correlation matrix can be plotted as a heatmap, to visualize the Pearson correlations between the different samples.
# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#Heatmap of all significant proteins
#The heatmap representation gives an overview of all significant proteins (rows) in all samples (columns). This allows to see general trends, for example if one sample or replicate is really different compared to the others. Additionally, the clustering of samples (columns) can indicate closer related samples and clustering of proteins (rows) indicates similarly behaving proteins. The proteins can be clustered by k-means clustering (kmeans argument) and the number of clusters can be defined by argument k.
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
#The heatmap shows a clustering of replicates and indicates that 4Ubi and 6Ubi enrich a similar repertoire of proteins. The k-means clustering of proteins (general clusters of rows) nicely separates protein classes with different binding behaviors.
#Alternatively, a heatmap can be plotted using the contrasts, i.e. the direct sample comparisons, as columns. Here, this emphasises the enrichment of ubiquitin interactors compared to the control sample.
# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

#Volcano plots of specific contrasts
#Volcano plots can be used to visualize a specific contrast (comparison between two samples). This allows to inspect the enrichment of proteins between the two samples (x axis) and their corresponding adjusted p value (y axis). The add_names argument can be set to FALSE if the protein labels should be omitted, for example if there are too many names.
# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "Ubi1_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "Ubi4_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "Ubi6_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep2, contrast = "Ubi1_vs_Ubi4", label_size = 2, add_names = TRUE) #うまくいかない
plot_volcano(dep2, contrast = "Ubi1_vs_Ubi6", label_size = 2, add_names = TRUE) #うまくいかない

#Barplots of a protein of interest
#It can also be useful to plot the data of a single protein, for example if this protein is of special interest.
# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = c("USP15", "IKBKG"))
# Plot a barplot for the protein USP15 with the data centered
plot_single(dep, proteins = "USP15", type = "centered")

#Frequency plot of significant proteins and overlap of conditions
#Proteins can be differentially enriched/expressed in multiple comparisons. To visualize the distribution of significant conditions per protein and the overlap between conditions, the plot_cond function can be used.
# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)
plot_cond(dep2)
plot_cond(dep3)

#Results table
#To extract a table containing the essential results, the get_results function can be used.
# Generate a results table
data_results <- get_results(dep)
# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
#The resulting table contains the following columns:
# Column names of the results table
colnames(data_results)
#Of these columns, the p.val and p.adj columns contain the raw and adjusted p values, respectively, for the contrast as depicted in the column name. The ratio columns contain the average log2 fold changes. The significant columns indicate whether the protein is differentially enriched/expressed, as defined by the chosen cutoffs. The centered columns contain the average log2 fold changes scaled by protein-wise centering.

#Generate a data.frame from the resulting SummarizedExperiment object
#You might want to obtain an ordinary data.frame of the results. For this purpose, the package provides functions to convert SummarizedExperiment objects to data.frames. get_df_wide will generate a wide table, whereas get_df_long will generate a long table.
# Generate a wide data.frame
df_wide <- get_df_wide(dep)
# Generate a long data.frame
df_long <- get_df_long(dep)

#Save your results object for reuse
#To facilitate future analysis and/or visualization of your current data, saving your analyzed data is highly recommended. We save the final data object (dep) as well as intermediates of the analysis, i.e. the initial SummarizedExperiment object (data_se), normalized data (data_norm), imputed data (data_imp) and differentially expression analyzed data (data_diff). This allows us to easily change parameters in future analysis.
# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")
# These data can be loaded in future R sessions using this command
load("data.RData")

#Workflow functions for the entire analysis
#The package contains workflow functions that entail the complete analysis and generate a report.
#LFQ-based DEP analysis
#Differential enrichment analysis of label-free proteomics data can be performed using the LFQ workflow function.
# The data is provided with the package 
data <- UbiLength
experimental_design <- UbiLength_ExpDesign
# The wrapper function performs the full analysis
data_results <- LFQ(data, experimental_design, fun = "MinProb", 
                    type = "control", control = "Ctrl", alpha = 0.05, lfc = 1)
#This wrapper produces a list of objects, which can be used to create a report and/or for further analysis. The report function produces two reports (pdf and html), a results table and a RData object, which are saved in a generated “Report” folder.
# Make a markdown report and save the results
report(data_results)
#The results generated by the LFQ function contain, among others, the results object (data.frame object) and the dep object (SummarizedExperiment object).
# See all objects saved within the results object
names(data_results)
#The results object contains the essential results and the dep object contains the full SummarizedExperiment object. The results table can be explored by selecting the $results object
# Extract the results table
results_table <- data_results$results
# Number of significant proteins
results_table %>% filter(significant) %>% nrow()
#The full data (dep object) can be used for the plotting functions as described in the chapter “Visualization of the results”, for example a heatmap.
# Extract the sign object
full_data <- data_results$dep
# Use the full data to generate a heatmap
plot_heatmap(full_data, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE)

#TMT-based DEP analysis
#Differential enrichment analysis of Tandem-Mass-Tag labeled proteomics data is also supported. Protein results tables from IsobarQuant can be directly analyzed using the TMT wrapper function.
# Need example data
TMTdata <- example_data
Exp_Design <- example_Exp_Design
# The wrapper function performs the full analysis
TMTdata_results <- TMT(TMTdata, expdesign = Exp_Design, fun = "MinProb",
                       type = "control", control = "Control", alpha = 0.05, lfc = 1)
#############################################################










#############################################################
