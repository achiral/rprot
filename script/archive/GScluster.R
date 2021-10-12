#GScluster
#https://github.com/unistbig/GScluster
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/GScluster")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#############################################################
#Open R program and type following commands in R console.
#install.packages('devtools') 
#library(remotes)
#install.packages("openxlsx")
#install.packages("XLConnect", type="source") # Mac
#library(devtools) 
#devtools::install_github("daattali/shinyjs")
#devtools::install_github('jhk0530/shinyCyJS')
#install.packages('shinyCyJS') #https://github.com/jhk0530/shinyCyJS
#install_github('unistbig/GScluster', force = T) 

#パッケージの読み込み
library(EnhancedVolcano)
#library(airway)
library(magrittr)
#library("DESeq2")
library(tidyverse) #ggplot2,dplyr使用
#tidyverse(ggplot2,dplyr),gcookbook読み込み
library(dplyr)
#library(stringr)
library(scales) #muted()関数使用のため
#library(ggcorrplot) #グラフ作成のため
#library(cowplot)
library(rJava)
library(openxlsx) #エクセル入出力（write.xlsx）
#library(XLConnect) #エクセル入出力,動作しない
#library(xlsx) #エクセル出力
#library(xlsx2) #エクセル出力,動作しない
library(readxl) #エクセル入力
#library(tablaxlsx) #エクセル表出力
library(gridExtra) #svg出力
#library(XLConnect) #動作しない
#library(rvg) #svg出力
#library(rsvg) #svg出力
#library(sets) #集合演算
library(GScluster) 
library(shinyCyJS) #https://stackoverflow.com/questions/46982593/r-shiny-error-there-is-no-package-called-shinyjs
#############################################################
#Example Run
#GScluster() 
#Read gene-set analysis result table.
#GSAresult=read.delim('https://github.com/unistbig/GScluster/raw/master/sample_geneset.txt', stringsAsFactors=FALSE)
#Read gene score table.
#GeneScores=read.delim('https://github.com/unistbig/GScluster/raw/master/sample_genescore.txt', header=F)
#Run GScluster
#GScluster(GSAresult = GSAresult, GeneScores = GeneScores, Species = 'H', alpha = 1, GsQCutoff = 0.25, GQCutoff = 0.05)
#エクセル入力
#read_excel(   path,  # Excelのファイルパス名
#              sheet = 1,  # シート名かシート番号
#              col_names = TRUE,  # ヘッダの有無
#              col_types = NULL,  # 列ごとのデータ型を定義。難しいです。
#              na = "",  # 欠損値をどの文字に置き換えるか 
#              skip = 0) # 読み飛ばす行数
dat1 <- read_excel("ClueGOResultTable_PC10_PFC_GSc.xls", 1) #シート1の読み込み
dat2 <- read_excel("ClueGOResultTable_PC10_PFC_GSc.xls", 2) #シート2の読み込み
dat3 <- read_excel("ClueGOResultTable_PC10_PFC_GSc.xls", 3) #シート3の読み込み
dat4 <- read_excel("ClueGOResultTable_PC10_PFC_GSc.xls", 4) #シート4の読み込み
dat5 <- read_excel("ClueGOResultTable_PC10_PFC_GSc.xls", 5) #シート5の読み込み
dat6 <- read_excel("ClueGOResultTable_PC10_PFC_GSc.xls", 6) #シート6の読み込み
datGS <- read.csv("Genescore.csv", stringsAsFactors=FALSE, header=FALSE) #ヘッダ無
datGS <- datGS[-1,] #1行目削除
options(digits=2) #change digit2桁表示指定
#列抽出
t(names(dat1)) #行列の名前確認
dat1 <- dat1[,c(2,12,7)] #列の順番を入れ替えて抽出
dat2 <- dat2[,c(2,12,7)] #列の順番を入れ替えて抽出
dat3 <- dat3[,c(2,12,7)] #列の順番を入れ替えて抽出
dat4 <- dat4[,c(2,12,7)] #列の順番を入れ替えて抽出
dat5 <- dat5[,c(2,12,7)] #列の順番を入れ替えて抽出
dat6 <- dat6[,c(2,12,7)] #列の順番を入れ替えて抽出
colnames(dat1) <- c("Geneset", "Gene", "qvalue") #列名変更
colnames(dat2) <- c("Geneset", "Gene", "qvalue") #列名変更
colnames(dat3) <- c("Geneset", "Gene", "qvalue") #列名変更
colnames(dat4) <- c("Geneset", "Gene", "qvalue") #列名変更
colnames(dat5) <- c("Geneset", "Gene", "qvalue") #列名変更
colnames(dat6) <- c("Geneset", "Gene", "qvalue") #列名変更
dat1["Geneset"] <- lapply(dat1["Geneset"], gsub, pattern=" ", replacement = "_") #スペースを_に置換
dat2["Geneset"] <- lapply(dat2["Geneset"], gsub, pattern=" ", replacement = "_") #スペースを_に置換
dat3["Geneset"] <- lapply(dat3["Geneset"], gsub, pattern=" ", replacement = "_") #スペースを_に置換
dat4["Geneset"] <- lapply(dat4["Geneset"], gsub, pattern=" ", replacement = "_") #スペースを_に置換
dat5["Geneset"] <- lapply(dat5["Geneset"], gsub, pattern=" ", replacement = "_") #スペースを_に置換
dat6["Geneset"] <- lapply(dat6["Geneset"], gsub, pattern=" ", replacement = "_") #スペースを_に置換
dat1$Gene <- toupper(dat1$Gene) #大文字変換
dat2$Gene <- toupper(dat2$Gene) #大文字変換
dat3$Gene <- toupper(dat3$Gene) #大文字変換
dat4$Gene <- toupper(dat4$Gene) #大文字変換
dat5$Gene <- toupper(dat5$Gene) #大文字変換
dat6$Gene <- toupper(dat6$Gene) #大文字変換
dat1 <- dat1 %>% filter(qvalue < 0.05) #q<0.05
dat2 <- dat2 %>% filter(qvalue < 0.05) #q<0.05
dat3 <- dat3 %>% filter(qvalue < 0.05) #q<0.05
dat4 <- dat4 %>% filter(qvalue < 0.05) #q<0.05
dat5 <- dat5 %>% filter(qvalue < 0.05) #q<0.05
dat6 <- dat6 %>% filter(qvalue < 0.05) #q<0.05
write.table(dat1,"dat1.txt", sep="\t", row.names=FALSE)
write.table(dat2,"dat2.txt", sep="\t", row.names=FALSE)
write.table(dat3,"dat3.txt", sep="\t", row.names=FALSE)
write.table(dat4,"dat4.txt", sep="\t", row.names=FALSE)
write.table(dat5,"dat5.txt", sep="\t", row.names=FALSE)
write.table(dat6,"dat6.txt", sep="\t", row.names=FALSE)
write.table(datGS,"datGS.txt", sep="\t", row.names=FALSE)
#Read gene-set analysis result table.
#GSAresult=read.delim('geneset_AML.txt', stringsAsFactors=FALSE)
GSAresult1=read.delim('dat1.txt', stringsAsFactors=FALSE)
GSAresult2=read.delim('dat2.txt', stringsAsFactors=FALSE)
GSAresult3=read.delim('dat3.txt', stringsAsFactors=FALSE)
GSAresult4=read.delim('dat4.txt', stringsAsFactors=FALSE)
GSAresult5=read.delim('dat5.txt', stringsAsFactors=FALSE)
GSAresult6=read.delim('dat6.txt', stringsAsFactors=FALSE)
#Read gene score table.
#GeneScores=read.delim('genescore_AML.txt', header=F)
datGS <- read.delim('datGS.txt', header=F)
GeneScores1 <- datGS[,c(1,2)]
GeneScores2 <- datGS[,c(1,3)]
GeneScores3 <- datGS[,c(1,4)]
GeneScores4 <- datGS[,c(1,2)]
GeneScores5 <- datGS[,c(1,3)]
GeneScores6 <- datGS[,c(1,4)]
#Run GScluster
GScluster(GSAresult = GSAresult1, 
          GeneScores = GeneScores1, 
          Species = 'H', 
          alpha = 1, 
          GsQCutoff = 0.25, 
          GQCutoff = 0.05)
GScluster(GSAresult = GSAresult2, 
          GeneScores = GeneScores2, 
          Species = 'M', 
          alpha = 1, 
          GsQCutoff = 0.25, 
          GQCutoff = 0.05)

