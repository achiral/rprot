#log2-Impute(MNAR)-Subtract(Median):like a Perseus
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/log2-Impute(MNAR)-Subtract(Median)_like a Perseus") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
##############################################################
#ライブラリ読み込み
library(DEP)
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(readxl) #エクセル入力(read_excel)
library(openxlsx) #エクセル入出力(write.xlsx)
##############################################################
#xlsx入力
rm(list = ls(all.names = TRUE))
data <- read_excel("AMY.xlsx", 3) #シート3入力
ExpDesign <- read_excel("AMY.xlsx", 2) #シート2入力
dim(data) #The data.frame dimensions:
colnames(data) #The data.frame column names:
##############################################################
#分割
split <- str_split(data$`Peak Name`, pattern = "\\|", simplify = TRUE)
colnames(split) <- c("sp", "Protein.IDs", "GeneName") #列名変更
class(split)
x <- data.frame(split)
#文字抽出
Protein.IDs <- str_sub(x$`Protein.IDs`, start = 1, end = 6) #`Peak Name`列の1-6文字目(Protein.IDs)抽出
Gene.names <- str_sub(x$`GeneName`, start = 1, end = -7) #`GeneName`列の1文字目〜-7文字目(GeneName)抽出
Species <- str_sub(x$`GeneName`, start = -5, end = -1) #`GeneName`列の-5〜-1文字目(Species)抽出
#文字抽出
data <- cbind(data, Protein.IDs) #dataとProtein.IDsを列ベクトル単位で結合 
data <- cbind(data, Gene.names) #dataとGene.namesを列ベクトル単位で結合
data <- cbind(data, Species) #dataとSpeciesを列ベクトル単位で結合
#Search Duplication
data$Protein.IDs %>% duplicated() %>% any()
data$Gene.names %>% duplicated() %>% any()
data$Species %>% duplicated() %>% any()
#Duplication table
data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
data %>% group_by(Species) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
#Unique Uniprot ID
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any() # Are there any duplicated names?
#SummarizedExperiment
Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
experimental_design <- ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
###############################################################################
#Log2-transform
data_se <- make_se(data_unique, Sample_columns, experimental_design) #columns=データ数, #Log2-transformation
#Impute:left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_se, fun = "man", shift = 1.8, scale = 0.3) #Perseus
data2 <- data.frame(data_imp_man@assays@data) #Subtract前
#Subtract(Median):Perseus
standardize <- function(z) {
  colmed <- apply(z, 2, median) #Median of Each Sample's Protein Expression level
  colmad <- apply(z, 2, mad)  # median absolute deviation
  rv <- sweep(z, 2, colmed,"-")  #subtracting median expression
  #rv <- sweep(rv, 2, colmad, "/")  # dividing by median absolute deviation
  return(rv)
}
data3 <- data2 #Subtract前
Sample_columns <- grep("(SC|PC)", colnames(data3)) # get Sample column numbers
data3[, Sample_columns] <- standardize(data3[Sample_columns]) #Subtract(Median)
data_imp_man@assays@data < array(data3) #スロットに戻す
data4 <- data.frame(data_imp_man@assays@data) #Subtract後
#############################################################
#txt出力
write.table (data3, file = "data3.txt", sep = "\t",quote = FALSE, row.names = TRUE) #Subtract前
#txt入力
data5 <- read.table("data3.txt",header=T, row.names=NULL) #Subtract前
#統合
data6 <- left_join(data, data5, by = c("Gene.names" = "row.names")) #元データ+Subtract前
#xlsx出力
smp <- list("data4"=data4,"data5"=data5,"data6"=data6) #リスト作成,Subtract後,Subtract前,元データ+Subtract前
write.xlsx(smp, "data4-6.xlsx") #シート出力
#############################################################
#Differential enrichment analysis:THSD
#multicomp.Rにて実施
#############################################################
plot_frequency(data_se)
plot_frequency(data_imp_man)
plot_numbers(data_se)
plot_numbers(data_imp_man)
plot_coverage(data_se)
plot_coverage(data_imp_man)
plot_normalization(data_se, data_imp_man)
plot_imputation(data_se, data_imp_man) 

#data_filt <- filter_missval(data_imp_man, thr = 0)
#plot_normalization(data_se, data_filt, data_imp_man)
#data_norm <- normalize_vsn(data_filt)
#plot_normalization(data_filt, data_norm, data_imp_man)
###############################################################################
###############################################################################
###############################################################################
#Analysis
###############################################################################
###############################################################################
###############################################################################
#Differential enrichment analysis:limma
data_diff <- test_diff(data_imp_man, type = "control", control = "SC0") # Test every sample versus control
data_diff_all_contrasts <- test_diff(data_imp_man, type = "all") # Test all possible comparisons of samples
#data_diff_manual <- test_diff(data_imp_man, type = "manual", test = c("SC10_vs_SC0","SC30_vs_SC0","PC0_vs_SC0","PC10_vs_SC0","PC30_vs_SC0")) # Test manually defined comparisons
#define cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))
dep2 <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1))
#dep3 <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(1.5))
###############################################################################
###############################################################################
#Differential enrichment analysis:THSD
#multicomp.Rにて実施
###############################################################################
#Visualization
#PCA plot
plot_pca(dep, x = 1, y = 2, indicate = "condition", label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE) # Plot the first and second principal components
plot_pca(dep2, x = 1, y = 2, indicate = "condition", label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE) # Plot the first and second principal components

#Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(dep2, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#Heatmap of all significant proteins:proteins (rows) in all samples (columns)
plot_heatmap(dep2, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
#Heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep2, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
#Volcano plots of specific contrasts:samples (x axis) adjusted p value (y axis)
plot_volcano(dep, contrast = "PC0_vs_SC0", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "PC10_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "PC30_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "SC10_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "SC30_vs_SC0", label_size = 2, add_names = TRUE)

#Barplots
plot_single(dep2, proteins = "SYN1") #protein of interest
plot_single(dep2, proteins = "SYN1", type = "centered") #data centered

#Frequency plot of significant proteins and overlap of conditions
plot_cond(dep)
plot_cond(dep2)

#Results table
data_results <- get_results(dep2)
# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
#The resulting table contains the following columns:
# Column names of the results table
colnames(data_results)
df_wide <- get_df_wide(dep2)
df_long <- get_df_long(dep2)

#SSave analyzed data
save(data_se, data_imp_man, data_diff, data_diff_all_contrasts, dep, dep2, file = "data.RData")
#xlsx出力
smp <- list("df_wide"=df_wide,"df_long"=df_long) #リスト作成
write.xlsx(smp, "AMY3.xlsx") #シート出力
# These data can be loaded in future R sessions using this command
load("data.RData")
#############################################################




###############################################################################
###############################################################################
###############################################################################
#DEP:log2-Filter(Missing)-Normalize(VSN)-Impute(MNAR)-Subtract(Median):
###############################################################################
###############################################################################
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/DEP2") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
##############################################################
###############################################################################
#xlsx入力
rm(list = ls(all.names = TRUE))
data <- read_excel("AMY.xlsx", 3) #シート3入力
ExpDesign <- read_excel("AMY.xlsx", 2) #シート2入力
dim(data) #The data.frame dimensions:
colnames(data) #The data.frame column names:
###############################################################################
#分割
split <- str_split(data$`Peak Name`, pattern = "\\|", simplify = TRUE)
colnames(split) <- c("sp", "Protein.IDs", "GeneName") #列名変更
class(split)
x <- data.frame(split)
#文字抽出
Protein.IDs <- str_sub(x$`Protein.IDs`, start = 1, end = 6) #`Peak Name`列の1-6文字目(Protein.IDs)抽出
Gene.names <- str_sub(x$`GeneName`, start = 1, end = -7) #`GeneName`列の1文字目〜-7文字目(GeneName)抽出
Species <- str_sub(x$`GeneName`, start = -5, end = -1) #`GeneName`列の-5〜-1文字目(Species)抽出
#文字抽出
data <- cbind(data, Protein.IDs) #dataとProtein.IDsを列ベクトル単位で結合 
data <- cbind(data, Gene.names) #dataとGene.namesを列ベクトル単位で結合
data <- cbind(data, Species) #dataとSpeciesを列ベクトル単位で結合
#Search Duplication
data$Protein.IDs %>% duplicated() %>% any()
data$Gene.names %>% duplicated() %>% any()
data$Species %>% duplicated() %>% any()
#Duplication table
data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
data %>% group_by(Species) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
#Unique Uniprot ID
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any() # Are there any duplicated names?
#SummarizedExperiment
Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
experimental_design <- ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
###############################################################################
#Log2-transform
data_se <- make_se(data_unique, Sample_columns, experimental_design) #columns=データ数, #Log2-transformation
###############################################################################
#Filter on missing values:Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)
#plot_frequency(data_se)
plot_frequency(data_filt)
#plot_numbers(data_se)
plot_numbers(data_filt)
#Plot a barplot of the protein identification overlap between samples
#plot_coverage(data_se)
plot_coverage(data_filt)
#Normalization:The data is background corrected and normalized by variance stabilizing transformation (vsn).
#data_norm <- normalize_vsn(data_se)
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm) #3つ以上も表示可能、ggplot2+で編集可能
###############################################################################
#Impute:left-shifted Gaussian distribution (for MNAR)
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)
plot_missval(data_norm)
plot_detect(data_filt)
plot_detect(data_norm)
impute(data_norm, fun = "")
#data_imp_man <- impute(data_filt, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
data2 <- data.frame(data_imp_man@assays@data) #Subtract前
plot_imputation(data_filt, data_norm, data_imp_man) #3つ以上も表示可能、ggplot2+で編集可能
plot_normalization(data_filt, data_norm, data_imp_man) #3つ以上も表示可能、ggplot2+で編集可能
###############################################################################
#Subtract(Median):Perseus
standardize <- function(z) {
  colmed <- apply(z, 2, median) #Median of Each Sample's Protein Expression level
  colmad <- apply(z, 2, mad)  # median absolute deviation
  rv <- sweep(z, 2, colmed,"-")  #subtracting median expression
  #rv <- sweep(rv, 2, colmad, "/")  # dividing by median absolute deviation
  return(rv)
}
data3 <- data2 #Subtract前
Sample_columns <- grep("(SC|PC)", colnames(data3)) # get Sample column numbers
data3[, Sample_columns] <- standardize(data3[Sample_columns]) #Subtract(Median)
data_imp_man@assays@data < array(data3) #スロットに戻す
data4 <- data.frame(data_imp_man@assays@data) #Subtract後
plot_normalization(data_filt, data_norm, data_imp_man) #3つ以上も表示可能、ggplot2+で編集可能
#############################################################
#txt出力
write.table (data3, file = "data3.txt", sep = "\t",quote = FALSE, row.names = TRUE) #Subtract前
#txt入力
data5 <- read.table("data3.txt",header=T, row.names=NULL) #Subtract前
#統合
data6 <- left_join(data, data5, by = c("Gene.names" = "row.names")) #元データ+Subtract前
#xlsx出力
smp <- list("data4"=data4,"data5"=data5,"data6"=data6) #リスト作成,Subtract後,Subtract前,元データ+Subtract前
write.xlsx(smp, "data4-6.xlsx") #シート出力
#############################################################
#plot
plot_frequency(data_se) #Log2-transform
plot_frequency(data_filt) #Filter(Missing)
plot_frequency(data_norm) #Normalize(VSN)
plot_frequency(data_imp_man) #Impute(MNAR)-Subtract(Median)
plot_numbers(data_se)
plot_numbers(data_filt)
plot_numbers(data_norm)
plot_numbers(data_imp_man)
plot_coverage(data_se)
plot_coverage(data_filt)
plot_coverage(data_norm)
plot_coverage(data_imp_man)
plot_normalization(data_se, data_filt, data_norm, data_imp_man)
plot_imputation(data_se, data_filt, data_norm, data_imp_man)
###############################################################################
###############################################################################
###############################################################################
#Analysis
###############################################################################
###############################################################################
###############################################################################
#Differential enrichment analysis:limma
# Test every sample versus control
data_diff <- test_diff(data_imp_man, type = "control", control = "SC0")
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp_man, type = "all")
# Test manually defined comparisons
#data_diff_manual <- test_diff(data_imp_man, type = "manual", test = c("SC10_vs_SC0","SC30_vs_SC0","PC0_vs_SC0","PC10_vs_SC0","PC30_vs_SC0"))
#define cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))
dep2 <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1))
#dep3 <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(1.5))
###############################################################################
#Visualization
#PCA plot
plot_pca(dep, x = 1, y = 2, indicate = "condition", label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE) # Plot the first and second principal components
plot_pca(dep2, x = 1, y = 2, indicate = "condition", label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE) # Plot the first and second principal components

#Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(dep2, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#Heatmap of all significant proteins:proteins (rows) in all samples (columns)
plot_heatmap(dep2, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
#Heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep2, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
#Volcano plots of specific contrasts:samples (x axis) adjusted p value (y axis)
plot_volcano(dep, contrast = "PC0_vs_SC0", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "PC10_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "PC30_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "SC10_vs_SC0", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "SC30_vs_SC0", label_size = 2, add_names = TRUE)

#Barplots
plot_single(dep, proteins = c("SYN1", "ANK3")) #protein of interest
plot_single(dep, proteins = "SYN1", type = "centered") #data centered

#Frequency plot of significant proteins and overlap of conditions
plot_cond(dep)
plot_cond(dep2)

#Results table
data_results <- get_results(dep2)
# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
#The resulting table contains the following columns:
# Column names of the results table
colnames(data_results)
df_wide <- get_df_wide(dep2)
df_long <- get_df_long(dep2)

#SSave analyzed data
save(data_se, data_norm, data_imp_man, data_diff, data_diff_all_contrasts, dep, dep2, file = "data.RData")
#xlsx出力
smp <- list("df_wide"=df_wide,"df_long"=df_long) #リスト作成
write.xlsx(smp, "AMY3.xlsx") #シート出力
# These data can be loaded in future R sessions using this command
load("data.RData")
#############################################################