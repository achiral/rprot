#ClueGO_05_heatmap_datapreparation.R
#ヒートマップ(データ整形)
#############################################################
setwd("/Users/user/Dropbox/0_Work/R/SWATH") #作業ディレクトリ設定
getwd() #作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
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
#library(rvg) #svg出力
#library(rsvg) #svg出力
#library(sets) #集合演算
#############################################################
#csv入力
dat1 <- read.csv("MCL1 default  node.csv", header=TRUE) #csv入力
dat2 <- read.csv("MCL2 default  node.csv", header=TRUE) #csv入力
dat3 <- read.csv("MCL3 default  node.csv", header=TRUE) #csv入力
dat4 <- read.csv("MCL4 default  node.csv", header=TRUE) #csv入力
dat5 <- read.csv("MCL5 default  node.csv", header=TRUE) #csv入力
options(digits=2) #change digit2桁表示指定
#xlsx入力
dat6 <- read_excel("heat.xlsx", 1) #シート1の読み込み
#dat6 <- read.delim("heat.txt", header=TRUE, row.names=1) #行名を読み込ませる
#列抽出
t(names(dat1)) #行列の名前確認
#dat1[,19]
dat1 <- dat1 %>% select(19) #列抽出
dat2 <- dat2 %>% select(19) #列抽出
dat3 <- dat3 %>% select(19) #列抽出
dat4 <- dat4 %>% select(19) #列抽出
dat5 <- dat5 %>% select(19) #列抽出
#列名変更
colnames(dat1) <- "GeneName" #列名変更
colnames(dat2) <- "GeneName" #列名変更
colnames(dat3) <- "GeneName" #列名変更
colnames(dat4) <- "GeneName" #列名変更
colnames(dat5) <- "GeneName" #列名変更
#置換
dat1$GeneName <- gsub("Eef1b2", "Eef1b", dat1$GeneName) #置換
dat5$GeneName <- gsub("Eno1b", "Eno1", dat5$GeneName) #置換
dat5$GeneName <- gsub("Gpi1", "Gpi", dat5$GeneName) #置換
#############################################################
#結合
dat1 <- left_join(dat1, dat6, by = c("GeneName" = "GeneName")) #結合
dat2 <- left_join(dat2, dat6, by = c("GeneName" = "GeneName")) #結合
dat3 <- left_join(dat3, dat6, by = c("GeneName" = "GeneName")) #結合
dat4 <- left_join(dat4, dat6, by = c("GeneName" = "GeneName")) #結合
dat5 <- left_join(dat5, dat6, by = c("GeneName" = "GeneName")) #結合
#############################################################
#txt出力
write.table(dat1, file="dat1.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat2, file="dat2.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat3, file="dat3.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat4, file="dat4.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat5, file="dat5.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat6, file="dat6.txt", sep="\t", row.names=F) #行名出力しない
#xlsx出力
#smp <- list("dat1"=dat1, "dat2"=dat2, "dat3"=dat3, "dat4"=dat4, "dat5"=dat5, "dat6"=dat6) #dat1-6リスト作成
#write.xlsx(smp, "heatCL.xlsx") #シート出力
#############################################################
#txt入力
dat1 <- read.table("dat1.txt", header=TRUE, sep="\t") #行名入力しない, row.names=1
dat2 <- read.delim("dat2.txt", header=TRUE) #行名入力しない
dat3 <- read.delim("dat3.txt", header=TRUE) #行名入力しない
dat4 <- read.delim("dat4.txt", header=TRUE) #行名入力しない
dat5 <- read.delim("dat5.txt", header=TRUE) #行名入力しない
dat6 <- read.delim("dat6.txt", header=TRUE) #行名入力しない
#xlsx入力
dat7 <- read_excel("FC_P.xlsx", 1) #シート1の読み込み
#############################################################
#Log2FC値符号変換
t(names(dat7)) #行列の名前確認
m <- mutate(dat7, Log2FC_SC10vsSC0 = Log2FC_SC0vsSC10 * -1)
m <- mutate(m, Log2FC_PC0vsSC0 = Log2FC_SC0vsPC0 * -1)
m <- mutate(m, Log2FC_PC10vsSC0 = Log2FC_SC0vsPC10 * -1)
m <- mutate(m, Log2FC_PC0vsSC10 = Log2FC_SC10vsPC0 * -1)
m <- mutate(m, Log2FC_PC10vsSC10 = Log2FC_SC10vsPC10 * -1)
dat7 <- mutate(m, Log2FC_PC10vsPC0 = Log2FC_PC0vsPC10 * -1)

#列名変更
names(dat7)[which(names(dat7)=="Group")] <- "Protein name"

#シート結合
dat1 <- left_join(dat1, dat7, by="GeneName")#識別子IDとしてdat1にdat7の情報を統合
dat2 <- left_join(dat2, dat7, by="GeneName")#識別子IDとしてdat2にdat7の情報を統合
dat3 <- left_join(dat3, dat7, by="GeneName")#識別子IDとしてdat3にdat7の情報を統合
dat4 <- left_join(dat4, dat7, by="GeneName")#識別子IDとしてdat4にdat7の情報を統合
dat5 <- left_join(dat5, dat7, by="GeneName")#識別子IDとしてdat5にdat7の情報を統合
dat6 <- left_join(dat6, dat7, by="GeneName")#識別子IDとしてdat6にdat7の情報を統合

#列名変更
names(dat1)[which(names(dat1)=="GeneName")] <- "Gene symbol"
names(dat2)[which(names(dat2)=="GeneName")] <- "Gene symbol"
names(dat3)[which(names(dat3)=="GeneName")] <- "Gene symbol"
names(dat4)[which(names(dat4)=="GeneName")] <- "Gene symbol"
names(dat5)[which(names(dat5)=="GeneName")] <- "Gene symbol"
names(dat6)[which(names(dat6)=="GeneName")] <- "Gene symbol"

#列抽出
t(names(dat1)) #行列の名前確認
dat1 <- dat1[,c(1,46,47,19:21,54,48,55,50,56,52,57,53,58,51,59,49)] #列抽出+並び替え
dat2 <- dat2[,c(1,46,47,19:21,54,48,55,50,56,52,57,53,58,51,59,49)] #列抽出+並び替え
dat3 <- dat3[,c(1,46,47,19:21,54,48,55,50,56,52,57,53,58,51,59,49)] #列抽出+並び替え
dat4 <- dat4[,c(1,46,47,19:21,54,48,55,50,56,52,57,53,58,51,59,49)] #列抽出+並び替え
dat5 <- dat5[,c(1,46,47,19:21,54,48,55,50,56,52,57,53,58,51,59,49)] #列抽出+並び替え
dat6 <- dat6[,c(1,46,47,19:21,54,48,55,50,56,52,57,53,58,51,59,49)] #列抽出+並び替え

#dat1 <- dat1[,c(1,34,35,19:21,24,22,23,28,26,27,32,30,31,36:41)] #列抽出+並び替え
#dat2 <- dat2[,c(1,34,35,19:21,24,22,23,28,26,27,32,30,31,36:41)] #列抽出+並び替え
#dat3 <- dat3[,c(1,34,35,19:21,24,22,23,28,26,27,32,30,31,36:41)] #列抽出+並び替え
#dat4 <- dat4[,c(1,34,35,19:21,24,22,23,28,26,27,32,30,31,36:41)] #列抽出+並び替え
#dat5 <- dat5[,c(1,34,35,19:21,24,22,23,28,26,27,32,30,31,36:41)] #列抽出+並び替え
#dat6 <- dat6[,c(1,34,35,19:21,24,22,23,28,26,27,32,30,31,36:41)] #列抽出+並び替え

#txt出力
write.table(dat1, file="dat1fcp.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat2, file="dat2fcp.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat3, file="dat3fcp.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat4, file="dat4fcp.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat5, file="dat5fcp.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat6, file="dat6fcp.txt", sep="\t", row.names=F) #行名出力しない