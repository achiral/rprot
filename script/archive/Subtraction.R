#Subtract Median
#https://stackoverflow.com/questions/50815358/how-to-subtract-a-median-using-complex-condition-in-r
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/DEP") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
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
df=structure(list(SKU = c(11202L, 11202L, 11202L, 11202L, 11202L, 
                          11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 
                          11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L
), stuff = c(8.85947691, 9.450108704, 10.0407405, 10.0407405, 
             10.63137229, 11.22200409, 11.22200409, 11.81263588, 12.40326767, 
             12.40326767, 12.40326767, 12.99389947, 13.58453126, 14.17516306, 
             14.76579485, 15.94705844, 17.12832203, 17.71895382, 21.26274458, 
             25.98779894, 63.19760196), action = c(0L, 0L, 0L, 0L, 0L, 0L, 
                                                   0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L), 
acnumber = c(137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 
             137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 
             137L, 137L, 137L), year = c(2018L, 2018L, 2018L, 2018L, 2018L, 
                                         2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 
                                         2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L)), .Names = c("SKU","stuff", "action", "acnumber", "year"), class = "data.frame", row.names = c(NA,-21L))
#############################################################
structure(list(SKU = c(11202L, 11202L, 11202L, 11202L, 11202L, 
                       11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 
                       11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L, 11202L
), stuff = c(8.85947691, 9.450108704, 10.0407405, 10.0407405, 
             10.63137229, 11.22200409, 11.22200409, 11.81263588, 12.40326767, 
             12.40326767, 12.40326767, 12.99389947, 13.58453126, 14.17516306, 
             14.76579485, 15.94705844, 17.12832203, 17.71895382, 21.26274458, 
             25.98779894, 63.19760196), action = c(0L, 0L, 0L, 0L, 0L, 0L, 
                                                   0L, 0L, 0L, 0L, 1L, NA, NA, NA, NA, NA, NA, NA, NA, 1L, 1L), 
acnumber = c(137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 
             137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 137L, 
             137L, 137L, 137L), year = c(2018L, 2018L, 2018L, 2018L, 2018L, 
                                         2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 
                                         2018L, 2018L, 2018L, 2018L, 2018L, 2018L, 2018L)), .Names = c("SKU", 
                                                                                                       "stuff", "action", "acnumber", "year"), class = "data.frame", row.names = c(NA,-21L))
#############################################################
df %>%
  group_by(SKU,acnumber,year) %>%
  summarize(value = 3*(median(stuff[action==1]) - median(stuff[match(1,action)-3:1])),
            stuff=first(stuff),
            action = sum(action)) %>%
  select(SKU,stuff,action,acnumber,year,value)
#############################################################



#My favorite commands Part3: ‘sweep’ function in R
#https://bioinfomagician.wordpress.com/2014/08/12/my-favorite-commands-part3-sweep-function-in-r/
#############################################################
#In the last two posts, I described useful commands either to expand or reduce your data. 
#In R, “apply” function is used to apply a function you specify to the data frame. 
#For example, you can use “mean” function to calculate mean values for each column and row.
data<-matrix(seq(1,12),ncol=4,nrow=3,byrow=TRUE)
data
colMean<-apply(data,2,mean)
colMean
rowMean<-apply(data,1,mean)
rowMean
#You see that to apply mean function by column-wise manner, you use “2” in the second argument 
#For row-wise mean, you use “1”.
#The third argument in “apply” can be a user-defined one, as long as it takes the right data.  For example,
myfunc<-function(x){
  t=1
  for (i in x){t=t*i}
  return(t)
}
#This function multiplies all the numbers in x and return the value.
multCol<-apply(data,2,myfunc)
multCol
#“apply” function is useful to this kind of simple column-wise or row-wise calculation, 
#however it would be cumbersome if you want to add different values for each column/row.


#############################################################
#sweep function in R
#############################################################
#R has a convenient function to apply different values to data in different columns/rows. 
#For example, you want to subtract “3”, “4”,”5″ ,”6″ from each value in the first, 2nd, 3rd and the last column. 
#You can do this by simply applying sweep function.
sweep(data,2,c(3,4,5,6),"-")
#sweep takes four arguments to work. 
#The first one is the data you want to modify. 
#The second one is the margin, if it is 1 then function is applied in row-wise manner. 
#If it is 2, it is applied in column-wise manner. 
#This is the same as “apply” function. 
#The third argument is a vector with the same length as column or row. 
#If the second argument is “1′, then the vector has to have the same length as row. 
#If it is “2”, the length of vector needs to be the same length as column. 
#The last argument is a binary operator, such as  – ,  + , * , / , < and >.

#############################################################
#Standardizing data matrix using sweep
#############################################################
#Gene expression data usually have large ranges and it is sometimes useful to standardize the data to see which genes go up or down more dramatically than others. 
#The mathematics is pretty simple:

#You have normalized expression data first, calculate median expression for each gene and then subtract the value and divide the results by median absolute deviation. 
#Since we are subtracting median expression which differs for each gene, sweep can come into play here. 
#By standardizing, all genes are centered to their median expression values and expression change is also normalized by dividing with median absolute deviation. 
#In this way, gene expression changes are not biased by expression levels.

standardize <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

data <- data.frame(data)
colnames(data) <- c("sample1","sample2","sample3","sample4")
rownames(data) <- c("gene1","gene2","gene3")
standardize(data)

#There is a pretty similar function called “scale” in R, 
#but it is slightly different from this function. 
#If “scale” function is used, you will get:
scale(t(data),scale=T)
#Note that you need to transpose the data because “scale” function works column wise-manner. 
#You see that values are not the same because “scale” function uses standard deviation instead of median absolute deviation. 
#To get identical results, you use a vector with median absolute deviation in the second argument for “scale” function.
scale(t(data),scale=apply(data,1,mad))

#############################################################
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
#For further analysis these proteins must get unique names. 
#Additionally, some proteins do not have an annotated gene name and for those we will use the Uniprot ID.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any() # Are there any duplicated names?
#Generate a SummarizedExperiment object
Sample_columns <- grep("(SAL|PCP)", colnames(data_unique)) # get Sample column numbers
experimental_design <- ExpDesign #ExperimentalDesignSheet(label,condition,replicate)
###############################################################################
#ここまではraw data
###############################################################################
data_se <- make_se(data_unique, Sample_columns, experimental_design) #columns=データ数, #Log2-transformation
data_se
plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_filt)
data2 <- data.frame(data_filt@assays@data) #Subtract前
###############################################################################
#log2-transformation #make_se関数で処理されるため不要
#data2 <- data
#Sample_columns <- grep("(SAL|PCP)", colnames(data2)) # get Sample column numbers
#Sample_columns
#data2[, Sample_columns] <- log(data2[Sample_columns], 2)
###############################################################################
#Subtract
#data_norm <- normalize_vsn(data_filt) #Subtruct?→Medianに変更する???
#standardize <- function(z) {
#  rowmed <- apply(z, 1, median) #Median of Each Protein level
#  rowmad <- apply(z, 1, mad)  # median absolute deviation
#  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  #rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
#  return(rv)
#}

standardize <- function(z) {
  colmed <- apply(z, 2, median) #Median of Each Sample's Protein Expression level
  colmad <- apply(z, 2, mad)  # median absolute deviation
  rv <- sweep(z, 2, colmed,"-")  #subtracting median expression
  #rv <- sweep(rv, 2, colmad, "/")  # dividing by median absolute deviation
  return(rv)
}



data3 <- data2
Sample_columns <- grep("(SC|PC)", colnames(data3)) # get Sample column numbers
Sample_columns
#data3[, Sample_columns] <- standardize(data3[Sample_columns])
data3[, Sample_columns] <- standardize(data3[Sample_columns])
#data_filt@assays@data <- SimpleList(data3) #スロットに戻す
data_filt@assays@data < array(data3) #スロットに戻す
data4 <- data.frame(data_filt@assays@data) #Subtract後
#############################################################
#scale(t(data3[Sample_columns]),scale=T)
#scale(t(data3[Sample_columns]),scale=apply(data3[Sample_columns],1,mad))
#############################################################
#txt出力
write.table (data3, file = "data3.txt", sep = "\t",quote = FALSE, row.names = TRUE)
#txt入力
data5 <- read.table("data3.txt",header=T, row.names=NULL)
#統合
data6 <- left_join(data, data5, by = c("Gene.names" = "row.names"))
#xlsx出力
smp <- list("data4"=data4,"data5"=data5,"data6"=data6) #リスト作成
write.xlsx(smp, "data4-6.xlsx") #シート出力
#############################################################
plot_coverage(data_filt)
#data_norm <- normalize_vsn(data_filt)
#plot_normalization(data_filt, data_norm) #3つ以上も表示可能、ggplot2+で編集可能
plot_normalization(data_filt)
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
impute(data_filt, fun = "")
#impute(data_norm, fun = "")
# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_filt, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
#data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) #Perseusはこれと思われる
#The effect of the imputation on the distributions can be visualized.
plot_imputation(data_filt, data_imp_man)
#plot_imputation(data_norm, data_imp_man)
###############################################################################
#Differential enrichment analysis
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
save(data_se, data_filt, data_imp_man, data_diff_all_contrasts, dep2, file = "data.RData")
#xlsx出力
smp <- list("df_wide"=df_wide,"df_long"=df_long) #リスト作成
write.xlsx(smp, "AMY3.xlsx") #シート出力
# These data can be loaded in future R sessions using this command
load("data.RData")
#############################################################