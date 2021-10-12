#setwd("/Users/user/Dropbox/0_Work/R/AMY") #作業ディレクトリ設定
setwd("~/Google Drive/マイドライブ/0_Work/R/AMY")
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
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
#dplyrパッケージによる共通列マージ(RCookbook2p171)
#library(dplyr)
data1 <- read_excel("ANOVA.xlsx", 1) #シート1の読み込み
data2 <- read_excel("ANOVA.xlsx", sheet = 2) #シート2の読み込み
data3 <- read_excel("ANOVA.xlsx", sheet = 3) #シート3の読み込み
data4 <- read_excel("ANOVA.xlsx", sheet = 4) #シート4の読み込み
data5 <- read_excel("ANOVA.xlsx", sheet = 5) #シート5の読み込み
data6 <- read_excel("ANOVA.xlsx", sheet = 6) #シート6の読み込み
#dataシートのマージ
data9 <- left_join(data1, data2, by="GN") #識別子IDとしてdata1とdata2の情報を統合
data9 <- left_join(data9, data3, by="GN") #
data9 <- left_join(data9, data4, by="GN") #
data9 <- left_join(data9, data5, by="GN") #
data9 <- left_join(data9, data6, by="GN") #
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
#library(openxlsx)
write.xlsx(data9, sheetname = "sheet1", file = "out_file1.xlsx") 

