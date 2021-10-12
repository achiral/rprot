#SWATH2stats
##############################################################
rm(list = ls(all.names = TRUE))
setwd("/Users/user/Dropbox/0_Work/R/SWATH2stats") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SWATH2stats")
BiocManager::install("aLFQ")
BiocManager::install("PECA")
#BiocManager::install("imsbInfer") #うまくいかない
#BiocManager::install("loadTransitonsMSExperiment") #うまくいかない
#BiocManager::install("readinteger_binary")
library(SWATH2stats)
library(MSstats)
library(data.table)
library(aLFQ)
library(PECA)
#library(imsbInfer) #うまくいかない
#library(loadTransitonsMSExperiment) #うまくいかない
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
#source('https://bioconductor.org/biocLite.R') 
#biocLite('SWATH2stats')
#############################################################
#Part 1: Loading and annotation
#############################################################
#Load the SWATH-MS example data from the package, this is a reduced file in order to limit the file size of the package.
library(SWATH2stats)
library(data.table)
data('Spyogenes', package = 'SWATH2stats')
#Alternatively the original file downloaded from the Peptide Atlas can be loaded from the working directory.
#data <- data.frame(fread('rawOpenSwathResults_1pcnt_only.tsv', sep='\t', header=TRUE))
#############################################################
#Extract the study design information from the file names. 
#Alternatively, the study design table can be provided as an external table.
Study_design <- data.frame(Filename = unique(data$align_origfilename))
Study_design$Filename <- gsub(".*strep_align/(.*)_all_peakgroups.*", "\\1", Study_design$Filename)
Study_design$Condition <- gsub("(Strep.*)_Repl.*", "\\1", Study_design$Filename)
Study_design$BioReplicate <- gsub(".*Repl([[:digit:]])_.*", "\\1", Study_design$Filename)
Study_design$Run <- seq(1:nrow(Study_design))
head(Study_design)
#############################################################
#The SWATH-MS data is annotated using the study design table.
data.annotated <- sample_annotation(data, Study_design, column.file = "align_origfilename")
#Remove the decoy peptides for a subsequent inspection of the data.
data.annotated.nodecoy <- subset(data.annotated, decoy==FALSE)
#############################################################
#Part 2: Analyze correlation, variation and signal
#############################################################
#Count the different analytes for the different injections.
count_analytes(data.annotated.nodecoy)
#Plot the correlation of the signal intensity.
correlation <- plot_correlation_between_samples(data.annotated.nodecoy, column.values = 'Intensity')
#Plot the correlation of the delta_rt, which is the deviation of the retention time from the expected retention time.
correlation <- plot_correlation_between_samples(data.annotated.nodecoy, column.values = 'delta_rt')
#Plot the variation of the signal across replicates.
variation <- plot_variation(data.annotated.nodecoy)
variation[[2]]
#Plot the total variation versus variation within replicates.
variation_total <- plot_variation_vs_total(data.annotated.nodecoy)
variation_total[[2]]
#Calculate the summed signal per peptide and protein across samples.
peptide_signal <- write_matrix_peptides(data.annotated.nodecoy) 
protein_signal <- write_matrix_proteins(data.annotated.nodecoy)
head(peptide_signal)
head(protein_signal)
#############################################################
#Part 3: FDR estimation
#############################################################
#Estimate the overall FDR across runs using a target decoy strategy.
par(mfrow = c(1, 3))
fdr_target_decoy <- assess_fdr_overall(data.annotated, n.range = 10, FFT = 0.25, output = 'Rconsole')
#FDR estimation one would need to filter the data with a lower mscore threshold to reach an overall protein FDR of 5%.
mscore4protfdr(data, FFT = 0.25, fdr_target = 0.05)
#############################################################
#Part 4: Filtering
#############################################################
#Filter data for values that pass the 0.001 mscore criteria in at least two replicates of one condition.
data.filtered <- filter_mscore_condition(data.annotated, 0.001, n.replica = 2)
#Select only the 10 peptides showing strongest signal per protein.
data.filtered2 <- filter_on_max_peptides(data.filtered, n_peptides = 10)
#Filter for proteins that are supported by at least two peptides.
data.filtered3 <- filter_on_min_peptides(data.filtered2, n_peptides = 2)
#############################################################
#Part 5: Conversion
#############################################################
#Convert the data into a transition-level format (one row per transition measured).
data.transition <- disaggregate(data.filtered3)
#Convert the data into the format required by MSstats.
MSstats.input <- convert4MSstats(data.transition)
head(MSstats.input)
#Convert the data into the format required by mapDIA.
mapDIA.input <- convert4mapDIA(data.transition) 
head(mapDIA.input)
#Convert the data into the format required by aLFQ.
aLFQ.input <- convert4aLFQ(data.transition)
head(aLFQ.input)
#Session info on the R version and packages used.
sessionInfo()
#############################################################
#############################################################

#############################################################
#2 Loading and annotating the data
#2.2 Loading the data
#############################################################
data('OpenSWATH_data', package='SWATH2stats')
data <- OpenSWATH_data
data('Study_design', package='SWATH2stats') 
head(Study_design)

# set working directory
#setwd("/Users/user/Dropbox/0_Work/R/SWATH2stats") #作業ディレクトリ設定
# Define input data file (e.g. OpenSWATH_output_file.txt)
#file.name <- 'OpenSWATH_output_file.txt'
# File name for annotation file
#annotation.file <- 'Study_design_file.txt'
# load data
#data <- data.frame(fread(file.name, sep='\t', header=TRUE))

# consult the manual page. > help(import_data)
# rename the columns
#data <- import_data(data)

# reduce number of columns
#data <- reduce_OpenSWATH_output(data)
# remove the iRT peptides (or other proteins)
#data <- data[grep('iRT', data$ProteinName, invert=TRUE),]
#############################################################
#2.3 Annotating the data
#############################################################
# list number and different Files present
nlevels(factor(data$filename))
levels(factor(data$filename))
# load the study design table from the indicated file
Study_design <- read.delim2(annotation.file, dec='.', sep='\t', header=TRUE)

# annotate data
data.annotated <- sample_annotation(data, Study_design)
head(unique(data$ProteinName))


# OPTIONAL: for human, shorten Protein Name to remove non-unique information
#(sp|Q9GZL7|WDR12_HUMAN --> Q9GZL7)
data$ProteinName <- gsub('sp\\|([[:alnum:]]+)\\|[[:alnum:]]*_HUMAN','\\1', data$ProteinName)
head(unique(data$ProteinName))
#############################################################
#3 Analyze data
#3.1 Count analytes
#############################################################
count_analytes(data.annotated)
#############################################################
#3.2 Plot correlation between samples
#############################################################
# Plot correlation of intensity
correlation <- plot_correlation_between_samples(data.annotated, column.values = "Intensity")
head(correlation)
# Plot correlation of retention times
correlation <- plot_correlation_between_samples(data.annotated, column.values = "RT")
#############################################################
#3.3 Plot variation
#############################################################
# plot variation of transitions for each condition across replicates
variation <- plot_variation(data.annotated)
head(variation[[2]])
# plot variation of summed signal for Peptides for each condition across replicates
variation <- plot_variation(data.annotated,Comparison = FullPeptideName + Condition ~ BioReplicate,
                            fun.aggregate = sum)
#############################################################
#3.4 Plot variation within replicates versus total variation
#############################################################
# plot variation of transitions for each condition within replicates > # compared to total variation
variation <- plot_variation_vs_total(data.annotated)
head(variation[[2]])
#############################################################
#3.5 Results on protein level
#############################################################
#Writing the overview matrix of summed intensities per protein entry per MS run:
protein_matrix <- write_matrix_proteins(data,filename = "SWATH2stats_overview_matrix_proteinlevel", rm.decoy = FALSE)
#############################################################
#3.6 Results on peptide level
#############################################################
peptide_matrix <- write_matrix_peptides(data,filename = "SWATH2stats_overview_matrix_peptidelevel",rm.decoy = FALSE)
#############################################################
#4 FDR estimation
#4.1 FDR: Overview and visualization
#FDR = (number of decoys * FFT)/(number of targets)
#############################################################
assess_decoy_rate(data)
# count decoys and targets on assay, peptide and protein level
# and report FDR at a range of m_score cutoffs
assess_fdr_overall(data, FFT = 0.7, output = "pdf_csv", plot = TRUE, filename='assess_fdr_overall_testrun')
# The results can be reported back to R for further calculations
overall_fdr_table <- assess_fdr_overall(data, FFT = 0.7,output = "Rconsole")
# create plots from fdr_table 
plot(overall_fdr_table, output = "Rconsole",filename = "FDR_report_overall")
# count decoys and targets on assay, peptide and protein level per run and report FDR at a range of m_score cutoffs
assess_fdr_byrun(data, FFT = 0.7, output = "pdf_csv", plot = TRUE, filename='assess_fdr_byrun_testrun')
# The results can be reported back to R for further calculations
byrun_fdr_cube <- assess_fdr_byrun(data, FFT = 0.7, output = "Rconsole")
# create plots from fdr_table 
plot(byrun_fdr_cube, output = "Rconsole",filename = "FDR_report_overall")
#############################################################
#4.2 Identification of useful m-score cutoffs to satisfy de- sired FDR criteria
#############################################################
# select and return a useful m_score cutoff in order to achieve the desired FDR quality for the entire table
mscore4assayfdr(data, FFT = 0.7, fdr_target=0.01)
#The function mscore4pepfdr reports an m-score cutoff to achieve a desired overall (global) peptide FDR:
# select and return a useful m_score cutoff in order to achieve the desired FDR quality for the entire table
mscore4pepfdr(data, FFT = 0.7, fdr_target=0.02)
# select and return a useful m_score cutoff in order to achieve the desired FDR quality for the entire table
mscore4protfdr(data, FFT = 0.7, fdr_target=0.02)
#############################################################
#5 Filtering the data
#5.1 Filter on m-score
#############################################################
data.filtered.mscore <- filter_mscore(data.annotated, 0.01)
data.filtered.mscore <- filter_mscore_freqobs(data.annotated, 0.01, 0.8,rm.decoy=FALSE)
data.filtered.mscore <- filter_mscore_condition(data.annotated, 0.01, 3)
data.filtered.fdr <- filter_mscore_fdr(data, FFT=0.7,overall_protein_fdr_target = 0.03, upper_overall_peptide_fdr_limit = 0.05)
#############################################################
#5.2 Filter on proteotypic peptides
#############################################################
data <- filter_proteotypic_peptides(data.filtered.mscore) 
data.all <- filter_all_peptides(data.filtered.mscore)
#############################################################
#5.3 Filter on number of peptides per protein
#############################################################
data.filtered.max <- filter_on_max_peptides(data.filtered.mscore, 5)
data.filtered.max.min <- filter_on_min_peptides(data.filtered.max, 2)
#############################################################
#6 Conversion of data for other tools
#6.1 Convert results to transition-level format
#6.1.1 Conversion within R
#############################################################
data.transition <- disaggregate(data)
head(data.transition)
#############################################################
#6.1.2 Conversion using a python script
#############################################################
data.python <- convert4pythonscript(data)
write.table(data.python, file="input.tsv", sep="\t", row.names=FALSE, quote=FALSE)
#The .tsv table can be transformed into a transition level table using a python script (as example featurealigner2msstats withRT.py from msproteomicstools which is available in the scripts folder of the package).
#python ./featurealigner2msstats.py input.csv output.csv
#Afterwards the generated .csv table is loaded again into R.
#data.transition <- data.frame(fread('output.csv',sep=',', header=TRUE))
data.transition <- data.frame(fread('input.tsv',sep='\t', header=TRUE))
#############################################################
#6.2 MSstats
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("MSstats")
#############################################################
#MSstats.input <- convert4MSstats(data.transition) 
data("OpenSWATH_data", package="SWATH2stats")
data("Study_design", package="SWATH2stats")
data <- sample_annotation(OpenSWATH_data, Study_design)
data.filtered.decoy <- filter_mscore(data, 0.01) 
raw <- disaggregate(data.filtered.decoy)
MSstats.input <- convert4MSstats(raw)

quantData <- dataProcess(MSstats.input)
#############################################################
#6.3 aLFQ
#############################################################
#aLFQ.input <- convert4aLFQ(data.transition)
BiocManager::install("aLFQ")
library(aLFQ)
aLFQ.input <- convert4aLFQ(raw)
prots <- ProteinInference(aLFQ.input, peptide_method = 'top',
                            peptide_topx = 3, peptide_strictness = 'loose', 
                            peptide_summary = 'mean', transition_topx = 3, 
                            transition_strictness = 'loose', transition_summary = 'sum', 
                            fasta = NA, model = NA, combine_precursors = FALSE)
#############################################################
#6.4 mapDIA
#############################################################
mapDIA.input <- convert4mapDIA(raw, RT =TRUE)
head(mapDIA.input)
write.table(mapDIA.input, row.names=FALSE, sep='\t')
#############################################################
#6.5 PECA
#############################################################
PECA.input <- convert4PECA(data) 
head(PECA.input)
library(PECA)
group1 <- c("Hela_Control_1", "Hela_Control_2", "Hela_Control_3")
group2 <- c("Hela_Treatment_1", "Hela_Treatment_2", "Hela_Treatment_3")
# PECA_df
results <- PECA_df(PECA.input, group1, group2, id="ProteinName", test = "rots")
#############################################################
#6.6 imsbInfer
#############################################################
data.annotated.full <- sample_annotation(OpenSWATH_data, Study_design)
data.annotated.full <- filter_mscore(data.annotated.full,mscore4assayfdr(data.annotated.full, 0.01))
data.annotated.full$decoy <- 0 ### imsbInfer needs the decoy column
#library(imsbInfer) #うまくいかない
#specLib <- loadTransitonsMSExperiment(data.annotated.full) #うまくいかない




