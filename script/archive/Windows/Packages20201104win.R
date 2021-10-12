#Rパッケージ
sessionInfo()
#############################################################
#パッケージインストール
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos='http://cran.us.r-project.org')
BiocManager::install("airway")
BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq")
BiocManager::install("DESeq2")
BiocManager::install("genefilter") #ヒートマップ
BiocManager::install("gplots") #ヒートマップ
BiocManager::install("ComplexHeatmap") #ヒートマップ
BiocManager::install("mzR")
BiocManager::install("DEP")
BiocManager::install("SWATH2stats")
BiocManager::install("aLFQ")
BiocManager::install("PECA")
BiocManager::install("BaylorEdPsych") #インストールできない
BiocManager::install("imsbInfer") #NA
BiocManager::install("loadTransitonsMSExperiment") #NA
BiocManager::install("readinteger_binary")
BiocManager::install("MSstats")
BiocManager::install("RCyjs") #cyREST,Cytoscape
#BiocManager::install() #更新
install.packages("VIM")
install.packages("imputeMissings")
install.packages("mice")
install.packages("mlbench")
install.packages("missForest")
install.packages("devtools")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggThemeAssist")
install.packages("ggcorrplot") 
install.packages("readxl") #エクセル入力
install.packages("openxlsx") #エクセル入力
install.packages("xlsx") #エクセル出力
#install.packages("xlsx2") #エクセル出力,NA
#install.packages("tablaxlsx") #エクセル表出力,NA
install.packages("XLConnect") #エクセル入出力
install.packages("writexl") #エクセル出力
install.packages("rJava")
install.packages("cowplot")
install.packages("rvg") #svgファイルで出力（https://www.karada-good.net/analyticsr/r-382）
install.packages("rsvg") #svgファイルで出力
install.packages("gridExtra") #svgファイルで出力
install.packages("sets") #集合演算
install.packages("multcomp") #多重比較検定
install.packages("agricolae") #
install.packages(c("gplots", "lme4", "reshape", "reshape2",
                   "ggplot2", "ggrepel", "data.table", "dplyr", "tidyr",
                   "survival", "doSNOW", "snow", "foreach", 'stringr',
                   "randomForest", "minpack.lm"), 
                 repos='http://cran.us.r-project.org')#MSstat
install.packages("BH")
install.packages("sgof") #bh(),Multiple Hypothesis Testing
#############################################################
#パッケージの読み込み
library(VIM)
#library(BaylorEdPsych) #NA
library(imputeMissings)
library(mice)
library(mlbench)
library(missForest)
library(SummarizedExperiment)
library(mzR)
library(DEP) #エラー?
library(SWATH2stats)
library(MSstats) #エラー
library(data.table)
library(aLFQ)
library(PECA)
#library(imsbInfer) #NA
#library(loadTransitonsMSExperiment) #NA

library(EnhancedVolcano)
library(airway)
library(magrittr)
library(DESeq)
library(DESeq2)
library(tidyverse) #ggplot2,dplyr使用
library(dplyr)
library(stringr)
library(scales) #muted()関数使用のため
library(scales) #muted()関数使用のため
library(rJava)
library(ggcorrplot) #グラフ作成のため
library(genefilter) #ヒートマップ
library(gplots) #ヒートマップ
library(ComplexHeatmap) #ヒートマップ
library(RColorBrewer) #色
library(cowplot)

library(readxl) #エクセル入力(read_excel)
library(openxlsx) #エクセル入出力(write.xlsx)
#library(XLConnect) #エクセル入出力,NA(JAVA8-11)
library(xlsx) #エクセル出力
#library(xlsx2) #エクセル出力,NA
#library(tablaxlsx) #エクセル表出力,NA
library(writexl) #エクセル出力
library(gridExtra) #svg出力
library(rvg) #svg出力
library(rsvg) #svg出力

library(sets) #集合演算
library(multcomp) #多重比較検定
library(BH) #FDR
library(sgof) #bh(),Multiple Hypothesis Testing
#############################################################
sessionInfo() 
#############################################################
