#microarray_volcanoplot.R
#############################################################
#setwd("/Users/user/Dropbox/0_Work/R/Directory_DEF") #作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/Directory_PCP") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
library(EnhancedVolcano)
library(magrittr)
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(dplyr)
library(scales) #muted()関数使用のため
library(readxl)
#############################################################
#Excelファイル読み込み
data1 <- read_excel("Stat.xlsx", 1) #シート1の読み込み
data2 <- read_excel("Stat.xlsx", 2) #シート2の読み込み
data3 <- read_excel("Stat.xlsx", 3) #シート3の読み込み
data4 <- read_excel("Stat.xlsx", 4) #シート4の読み込み
data5 <- read_excel("Stat.xlsx", 5) #シート5の読み込み
#############################################################
#ggplot2
plot1 <- ggplot(data1, aes(x = log2FC, y = p, colour = p)) +
  geom_point(size = 0.5)+
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_reverse()
plot1 +
  scale_colour_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 0.05
  ) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_reverse()

plot1_mod <- ggplot(data1_mod, aes(x = log2FC, y = log10p, colour = log10p)) +
  geom_point(size = 0.5)+
  scale_x_continuous(limits = c(-2, 2))
plot1_mod +
  scale_colour_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 2
  ) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_reverse()

#log10pを計算して保存
data1_mod <- data1 %>%
  mutate(log10p = -log10(p))
data1_mod

#p<0.05, log2FC>1で分類
data1_mod2 <- data1_mod %>%
  mutate(p0.05 = if_else(p < 0.05, "lower", "higher"),
         log2FC1 = if_else(log2FC > 1, "higher", "lower")
  )
data1_mod2

data1_mod3 <- data1_mod2 %>%
  filter(p0.05 == "lower") %>%
  mutate(p0.05log2FC1 = if_else(log2FC > 1, "p<0.05, FC>2", "p>=0.05, FC=<2")
  )
data1_mod3
str(data1_mod3) #データ構造の確認

data1_mod4 <- data1_mod3 %>%
  mutate(type_of_gene = as.factor(type_of_gene),
         p0.05 = as.factor(p0.05),
         log2FC1 = as.factor(log2FC1),
         p0.05log2FC1 = as.factor(p0.05log2FC1)
  ) #ファクタに変換
data1_mod4
str(data1_mod4) #データ構造の確認

plot1_mod2 <- ggplot(data1_mod3, aes(x = log2FC, y = log10p, colour = p0.05log2FC1, fill = p0.05log2FC1)) +
  geom_point(size = 1, shape = 21, alpha = 0.8) +
  scale_fill_manual(
    values = c("red", "black"),
    guide = guide_legend(override.aes = list(shape = 21))) +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(limits = c(-6.5, 6.5))
plot1_mod2

plot1_mod3 <- ggplot(data1_mod3, aes(x = log2FC, y = log10p, colour = log10p)) +
  geom_point(size = 3, shape = 20, alpha = 0.8) +
  scale_colour_gradient2(low = "black", 
                         mid = "black",
                         high = "red",
                         midpoint = 1.30102999566398) +
  scale_x_continuous(limits = c(-6.5, 6.5))
plot1_mod3
#############################################################
#EnhancedVolcano############################################
#EnhancedVolcano(data1_mod4, #データを指定
#                lab = data1_mod4$Gene, #遺伝子名を指定、rownames(データ名)のみ？
#                x = 'log2FC', #x軸の数値に使う列の指定
#                y = 'p', #y軸の数値に使う列の指定
#                xlim = c(-4, 6.5)) #x軸の表示範囲を設定、なくても良い

vol1 <- EnhancedVolcano(data1, #データを指定
                lab = NA, #遺伝子名を指定、rownames(データ名)のみ？
                x = 'log2FC', #x軸の数値に使う列の指定
                y = 'p', #y軸の数値に使う列の指定
                xlim = c(-2, 2), #x軸の表示範囲を設定、なくても良い
                ylim = c(0, 5), #y軸の表示範囲を設定、なくても良い
                pCutoff = 5.0*10e-2,#p値のカットオフ値を設定
                FCcutoff = 1.3,#発現量（a倍など）のカットオフ値を設定
                ) 
vol1


vol2 <- EnhancedVolcano(data2, #データを指定
                lab = NA, #遺伝子名を指定、rownames(データ名)のみ？
                x = 'log2FC', #x軸の数値に使う列の指定
                y = 'p', #y軸の数値に使う列の指定
                xlim = c(-3.8, 3.8), #x軸の表示範囲を設定、なくても良い
                ylim = c(0, 5), #y軸の表示範囲を設定、なくても良い
                pCutoff = 5.0*10e-2,#p値のカットオフ値を設定
                FCcutoff = 1.3, #発現量（a倍など）のカットオフ値を設定
                ) 
vol2

vol3 <- EnhancedVolcano(data3, #データを指定
                lab = NA, #遺伝子名を指定、rownames(データ名)のみ？
                x = 'log2FC', #x軸の数値に使う列の指定
                y = 'p', #y軸の数値に使う列の指定
                xlim = c(-3.1, 3.1), #x軸の表示範囲を設定、なくても良い
                ylim = c(0, 5), #y軸の表示範囲を設定、なくても良い
                pCutoff = 5.0*10e-2,#p値のカットオフ値を設定
                FCcutoff = 1.3, #発現量（a倍など）のカットオフ値を設定
                ) 
vol3

vol4 <- EnhancedVolcano(data4, #データを指定
                lab = NA, #遺伝子名を指定、rownames(データ名)のみ？
                x = 'log2FC', #x軸の数値に使う列の指定
                y = 'p', #y軸の数値に使う列の指定
                xlim = c(-3, 3), #x軸の表示範囲を設定、なくても良い
                ylim = c(0, 5), #y軸の表示範囲を設定、なくても良い
                pCutoff = 5.0*10e-2,#p値のカットオフ値を設定
                FCcutoff = 1.3, #発現量（a倍など）のカットオフ値を設定
                ) 
vol4

vol5 <- EnhancedVolcano(data5, #データを指定
                lab = NA, #遺伝子名を指定、rownames(データ名)のみ？
                x = 'log2FC', #x軸の数値に使う列の指定
                y = 'p', #y軸の数値に使う列の指定
                xlim = c(-3, 3), #x軸の表示範囲を設定、なくても良い
                ylim = c(0, 5), #y軸の表示範囲を設定、なくても良い
                pCutoff = 5.0*10e-2,#p値のカットオフ値を設定
                FCcutoff = 1.3, #発現量（a倍など）のカットオフ値を設定
                ) 
vol5

#P359##############################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("vol1.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol1) #プロット作成
dev.off() #PDF出力
pdf("vol2.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol2) #プロット作成
dev.off() #PDF出力
pdf("vol3.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol3) #プロット作成
dev.off() #PDF出力
pdf("vol4.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol4) #プロット作成
dev.off() #PDF出力
pdf("vol5.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol5) #プロット作成
dev.off() #PDF出力
##############################################




