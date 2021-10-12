#ANOVA, TukeyHSD
#https://github.com/gmiaslab/MathematicaBioinformatics/blob/master/sandberg-sampledata.txt
#https://stat.ethz.ch/pipermail/bioconductor/2013-May/052446.html
#sandberg-sampledata.txt
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/ANOVA") #作業ディレクトリ設定
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
sdata<-read.table("sandberg-sampledata.txt", header=T, row.names=1)
strain <- gl(2,12,24, label=c("129","bl6")) #カテゴリー2,繰返し12,要素24,ラベル2種
region <- gl(6,2,24, label=c("ag", "cb", "cx", "ec", "hp", "mb")) #部位6,繰返し2,要素24(12x2リピート),ラベル6種
# define ANOVA function ##  ANOVA
aof <- function(x) { 
  m <- data.frame(strain,region, x); 
  anova(aov(x ~ strain + region + strain*region, m))
}
# apply analysis to the data and get the pvalues.
anovaresults <- apply(sdata, 1, aof)
pvalues <- data.frame(lapply(anovaresults, function(x) { x["Pr(>F)"][1:3,] }))
pvalues2 <- data.frame(t(pvalues))
colnames(pvalues2) <- c("p_strain","p_region","p_interaction")
write.table(pvalues2, file="anova-results.txt", quote=F, sep='\t')

# Get the genes with good region effect pvalues. 
reg.hi.p <- sort(t(data.frame(pvalues[2, pvalues[2,] < 0.0001 & pvalues[3,] > 0.1])))
reg.hi.pdata <- sdata[row.names(reg.hi.p),]
cbtempl <- c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0) #cerebellum
cor.test(t(sdata[14,]),cbtempl)
template.match <- function(x, template) {
  k <- cor.test(x,template)
  k$p.value
}

#############################################################
#上手く行かない
#############################################################
cbtempl.results <- apply(reg.hi.pdata, 1, template.match, cbtempl) 
write.table(cbtempl.results, file="cbtempl-results.txt", quote=F, sep='\t')

t.test(t(reg.hi.pdata[14,region!="cb"]), t(reg.hi.pdata[14,region=="cb"]), var.equal=T)

regionp <- sort(t(pvalues[2,])) # the p values must be sorted in increasing order! 
fdr.result <- bh(regionp, 0.05) # reports that 192 genes are selected
bhthresh <- cbind(regionp, fdr.result) # put the bh results in our table. 
write.table(bhthresh, "bhthresh.txt", sep='\t', quote=F) # print to a file.


anovaresults2 <- data.frame(anovaresults)
mcps <- TukeyHSD(anovaresults, ordered=TRUE) 
Error in UseMethod("TukeyHSD") :  no applicable method for 'TukeyHSD' applied to an object of class "list"

