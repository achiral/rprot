#MergeFunction
#臨床情報シートに新規情報を追加
#SWATHデータの統合
#merge(dataA,dataB,by=”キー変数”,incomparables=NA)
#############################################################
#作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Nicotine")
setwd("/Users/user/Dropbox/0_Work/R/SWATH")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#############################################################
#Excelファイル読み込み
#install.packages("readxl")
#y
library(readxl)
data1 <- read_excel("Info.xlsx", 1) #シート1の読み込み
data2 <- read_excel("Gx.xlsx", sheet = 1) #シート1の読み込み
#dataシートのマージ
#data3 <- merge(data1, data2, by="ID",incomparables=NA) #識別子IDとしてdata1とdata2に共通する情報を統合,NA対象外(full_joinと同じ)
#dplyrパッケージによる共通列マージ(RCookbook2p171)
library(dplyr)
data4 <- left_join(data1, data2, by="ID")#識別子IDとしてdata1にdata2の情報を統合
view(data4)
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(data4, sheetname = "sheet1", file = "out_file.xlsx") 
#############################################################
#Excelファイル読み込み
#install.packages("readxl")
#y
library(readxl)
data1 <- read_excel("PC10_matrix10merge.xlsx", 1) #シート1の読み込み
#data2 <- read.table("PC10_matrix23.txt", header=T, sep="\t") #ヘッダあり
data2 <- read_excel("PC10_matrix23.xlsx", 1) #シート1の読み込み
data3 <- read_excel("PC10_matrix25.xlsx", 1) #シート1の読み込み
data4 <- read_excel("PC10_matrix27.xlsx", 1) #シート1の読み込み
data5 <- read_excel("PC10_matrix29.xlsx", 1) #シート1の読み込み
#dataシートのマージ
#merge <- merge(data1, data2, by="ID",incomparables=NA) #識別子IDとしてdata1とdata2に共通する情報を統合,NA対象外(full_joinと同じ)
#dplyrパッケージによる共通列マージ(RCookbook2p171)
#install.packages("dplyr")
library(dplyr)
data1m <- left_join(data1, data2, by="ID")#識別子IDとしてdata1にdata2の情報を統合
data1m <- left_join(data1m, data3, by="ID")#識別子IDとしてdata1にdata2の情報を統合
data1m <- left_join(data1m, data4, by="ID")#識別子IDとしてdata1にdata2の情報を統合
data1m <- left_join(data1m, data5, by="ID")#識別子IDとしてdata1にdata2の情報を統合
View(data1m)
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(data1m, sheetname = "sheet1", file = "out_file.xlsx") 