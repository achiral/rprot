#Fisher
#Fisher's exact test
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html
#fisher.test(x, y = NULL, workspace = 200000, hybrid = FALSE,
#         hybridPars = c(expect = 5, percent = 80, Emin = 1),
#         control = list(), or = 1, alternative = "two.sided",
#         conf.int = TRUE, conf.level = 0.95,
#         simulate.p.value = FALSE, B = 2000)
#############################################################
setwd("/Users/user/Dropbox/0_Work/R/Stat") #作業ディレクトリ設定
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
#Excelファイル読み込み
data1 <- read_excel("F.xlsx", 1) #シート1の読み込み
data2 <- read_excel("F.xlsx", sheet = 2) #シート2の読み込み
options(digits=2) #change digit2桁表示指定
#str(data1) #データ構造確認
#str(data2) #データ構造確認
#############################################################
fisher.test(data1, y = NULL, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = FALSE, B = 2000)
fisher.test(data2)
#############################################################
summary(data1)
summary(data2)
#############################################################