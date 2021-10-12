#Correlation
#############################################################
#作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Directory_DEF")
#setwd("/Users/user/Dropbox/0_Work/R/Directory_PCP")
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/Directory_PCP")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#############################################################
#Download the package from Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("airway")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("DESeq")
#BiocManager::install("DESeq2")
#install.packages("tidyverse")
#install.packages("ggThemeAssist")

#Load the package into R session
#Quick start
#library(EnhancedVolcano)
#library(airway)
library(magrittr)
#library("DESeq2")
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(dplyr)
library(scales) #muted()関数使用のため

#Excelファイル読み込み
#install.packages("readxl")
library(readxl)
data6 <- read_excel("Stat.xlsx", 6) #シート6の読み込み
data7 <- read_excel("Stat.xlsx", sheet = 7) #シート7の読み込み
#data3 <- read_excel("Stat.xlsx", sheet = "datafile3") #シート名で読み込み
#data4 <- read_excel("Stat.xlsx", sheet = "datafile4", col_names = TRUE, col_types = c("blank", "text", "date", "text", "numeric"))  #最初の列を削除(blank)し、後続3列の型を指定



#log2FCで分類(DEF)
#data6_mod1 <- data6 %>%
#  mutate(UpDown = if_else(log2FCDEF2h < 0 & log2FCDEF24h < 0, "DD", if_else(log2FCDEF2h >= 0 & log2FCDEF24h >= 0, "UU", "UD")),
#         )
#write.table (data6_mod1, file = "/Users/user/Dropbox/0_Work/R/Directory_DEF/output.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, append = FALSE)
#str(data6_mod1) #データ構造の確認

#log2FCで分類(PCP)
data6_mod1 <- data6 %>%
  mutate(UpDown = if_else(log2FCPCP2h < 0 & log2FCPCP24h < 0, "DD", if_else(log2FCPCP2h >= 0 & log2FCPCP24h >= 0, "UU", "UD")),
  )
#write.table (data6_mod1, file = "/Users/user/Dropbox/0_Work/R/Directory_PCP/output.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, append = FALSE)
#write.table (data6_mod1, file = "/Users/user/Dropbox/0_Work/R/Directory_PCP/output.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, append = FALSE)
write.table (data6_mod1, file = "/Users/user/Dropbox/My Mac (MacBook-Pro.local)/Desktop/Directory_PCP/output.txt", sep = "\t",
quote = FALSE, row.names = FALSE, append = FALSE)
str(data6_mod1) #データ構造の確認

#protein_codingでフィルタ:エクセル置換で-を_に修正済み################################
data6_mod1$ToG <- recode(data6_mod1$type_of_gene, protein_coding = "1", ncRNA ="2", pseudo = "3", rRNA = "4", unknown = "5") #数字に置き換えて新しい列を追加
data6_mod1$ToG <- as.numeric(as.character(data6_mod1$ToG)) #一旦文字列にしてから数値化
str(data6_mod1) #データ構造の確認

data6_mod2 <- data6_mod1 %>%
  filter(data6_mod1$ToG == 1) #protein_coding(1)でフィルタ
data6_mod2

#log2FCで分類(DEF)
#data7_mod1 <- data7 %>%
#  mutate(UpDown = if_else(log2FCDEF2h < 0 & log2FCDEF24h < 0, "DD", if_else(log2FCDEF2h >= 0 & log2FCDEF24h >= 0, "UU", "UD")),
#  )
#write.table (data7_mod1, file = "/Users/user/Dropbox/0_Work/R/Directory_DEF/output.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, append = FALSE)
#str(data7_mod1) #データ構造の確認

#log2FCで分類(PCP)
data7_mod1 <- data7 %>%
  mutate(UpDown = if_else(log2FCPCP2h < 0 & log2FCPCP24h < 0, "DD", if_else(log2FCPCP2h >= 0 & log2FCPCP24h >= 0, "UU", "UD")),
  )
write.table (data7_mod1, file = "/Users/user/Dropbox/My Mac (MacBook-Pro.local)/Desktop/Directory_PCP/output.txt", sep = "\t",
             quote = FALSE, row.names = FALSE, append = FALSE)
str(data7_mod1) #データ構造の確認

#protein_codingでフィルタ:エクセル置換で-を_に修正済み################################
data7_mod1$ToG <- recode(data7_mod1$type_of_gene, protein_coding = "1", ncRNA ="2", pseudo = "3", rRNA = "4", unknown = "5") #数字に置き換えて新しい列を追加
data7_mod1$ToG <- as.numeric(as.character(data7_mod1$ToG)) #一旦文字列にしてから数値化
str(data7_mod1) #データ構造の確認

data7_mod2 <- data7_mod1 %>%
  filter(data7_mod1$ToG == 1) #protein_coding(1)でフィルタ
data7_mod2

################################################################################################
#data6_mod2 <- data6_mod1 %>%
#  mutate(UpDown = as.factor(UpDown)) #ファクタに変換
#data6_mod2 <- read.table("output.txt", header = TRUE, sep = "\t") #テキスト読み込み
################################################################################################

#作図################################################################################################
#DEF
#plot6 <- ggplot(data6_mod2, aes(x = log2FCDEF24h, y = log2FCDEF2h, colour = UpDown,fill = UpDown)) +
#  geom_point(shape = 21, size = 2.5, alpha = 0.5) + #基本プロット
#  scale_colour_manual(values = c("blue", "black", "red"), guide = FALSE) +
#  scale_fill_manual(values = c("blue", "black", "red"), guide = FALSE) #凡例なし

#PCP
plot6 <- ggplot(data6_mod2, aes(x = log2FCPCP24h, y = log2FCPCP2h, colour = UpDown,fill = UpDown)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.5) + #基本プロット
  scale_colour_manual(values = c("blue", "black", "red"), guide = FALSE) +
  scale_fill_manual(values = c("blue", "black", "red"), guide = FALSE) #凡例なし

#Correlation curve
corr6 <- plot6 +
  xlim(-4.5,4.5) + ylim(-4.5,4.5) + #x軸y軸範囲指定
  theme_classic() + #クラシックに上書きしていく
  theme(
    axis.line = element_line(colour = "black", size = 1), #xy軸の指定
    panel.grid.major = element_line(colour = "grey", size = 0.5),
    panel.grid.minor = element_line(colour = "grey", linetype = "dashed", size = 0.5)
  )
corr6

#DEF
#plot7 <- ggplot(data7_mod2, aes(x = log2FCDEF24h, y = log2FCDEF2h, colour = UpDown,fill = UpDown)) +
#  geom_point(shape = 21, size = 2.5, alpha = 0.5) + #基本プロット
#  scale_colour_manual(values = c("blue", "black", "red"), guide = FALSE) +
#  scale_fill_manual(values = c("blue", "black", "red"), guide = FALSE) #凡例なし

#PCP
plot7 <- ggplot(data7_mod2, aes(x = log2FCPCP24h, y = log2FCPCP2h, colour = UpDown,fill = UpDown)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.5) + #基本プロット
  scale_colour_manual(values = c("blue", "black", "red"), guide = FALSE) +
  scale_fill_manual(values = c("blue", "black", "red"), guide = FALSE) #凡例なし

#Correlation curve
corr7 <- plot7 +
  xlim(-3,3) + ylim(-3,3) + #x軸y軸範囲指定
  theme_classic() + #クラシックに上書きしていく
  theme(
    axis.line = element_line(colour = "black", size = 1), #xy軸の指定
    panel.grid.major = element_line(colour = "grey", size = 0.5),
    panel.grid.minor = element_line(colour = "grey", linetype = "dashed", size = 0.5)
  )
corr7

#P359#########################################################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("corr6.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(corr6) #プロット作成
dev.off() #PDF出力
pdf("corr7.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(corr7) #プロット作成
dev.off() #PDF出力
##############################################################################################################
