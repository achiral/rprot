#01_volcanoplot_ClueGO.R
#ClueGO出力データからvolcano plot描写
#重複削除
#https://a-habakiri.hateblo.jp/entry/2016/11/29/215013
#############################################################
setwd("/Users/user/Dropbox/0_Work/R/SWATH") #作業ディレクトリ設定
getwd() #作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
library(EnhancedVolcano)
library(magrittr)
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(dplyr)
library(scales) #muted()関数使用のため
library(gridExtra) #svg出力のため
library(rJava)
library(readxl) #エクセル読み込み
library(openxlsx) #JAVA不使用で大きなデータも読み込める
library(cowplot)
#############################################################
#Excelファイル読み込み
data1 <- read_excel("ClueGOResultTable_PC10_PFC.xls", 1) #シート1の読み込み
data2 <- read_excel("ClueGOResultTable_PC10_PFC.xls", 2) #シート2の読み込み
data3 <- read_excel("ClueGOResultTable_PC10_PFC.xls", 3) #シート3の読み込み
data4 <- read_excel("ClueGOResultTable_PC10_PFC.xls", 4) #シート4の読み込み
data5 <- read_excel("ClueGOResultTable_PC10_PFC.xls", 5) #シート5の読み込み
data6 <- read_excel("ClueGOResultTable_PC10_PFC.xls", 6) #シート6の読み込み
#############################################################
#行番号挿入,重複削除
data1_mod1 <- data1 %>% mutate(No = row_number()) %>% #行番号挿入
  distinct(GOTerm,.keep_all=TRUE) #重複削除
#write.xlsx(data1_mod1, "out1.xlsx", sheetname="sheet1") #出力
#str(data1_mod1)
data2_mod1 <- data2 %>% mutate(No = row_number()) %>% #行番号挿入
  distinct(GOTerm,.keep_all=TRUE) #重複削除
str(data2_mod1)
data3_mod1 <- data3 %>% mutate(No = row_number()) %>% #行番号挿入
  distinct(GOTerm,.keep_all=TRUE) #重複削除
data4_mod1 <- data4 %>% mutate(No = row_number()) %>% #行番号挿入
  distinct(GOTerm,.keep_all=TRUE) #重複削除
data5_mod1 <- data5 %>% mutate(No = row_number()) %>% #行番号挿入
  distinct(GOTerm,.keep_all=TRUE) #重複削除
data6_mod1 <- data6 %>% mutate(No = row_number()) %>% #行番号挿入
  distinct(GOTerm,.keep_all=TRUE) #重複削除

#p値→-log10p値変換
data1_mod2 <- data1_mod1 %>%
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))
data2_mod2 <- data2_mod1 %>%
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))
data3_mod2 <- data3_mod1 %>%
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))
data4_mod2 <- data4_mod1 %>%
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))
data5_mod2 <- data5_mod1 %>%
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))
data6_mod2 <- data6_mod1 %>%
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))

#Inf(無限大)の置換,ClueGOでは-Log10FDRが45くらいで頭打ち?
ifelse((is.infinite(data1_mod2$log10p) & (data1_mod2$log10p > 0)), 45, data1_mod2$log10p) #Infの場合に50に置換したい
replace(data1_mod2$log10p, which(is.infinite(data1_mod2$log10p) & (data1_mod2$log10p > 0)), 50) #Infの場合に50に置換
data1_mod2$log10p[is.infinite(data1_mod2$log10p)& (data1_mod2$log10p > 0)] <- 50 #置換の実行
#data1_mod2$log10p #置換結果の表示
ifelse((is.infinite(data2_mod2$log10p) & (data2_mod2$log10p > 0)), 45, data2_mod2$log10p) #Infの場合に50に置換したい
replace(data2_mod2$log10p, which(is.infinite(data2_mod2$log10p) & (data2_mod2$log10p > 0)), 50) #Infの場合に50に置換
data2_mod2$log10p[is.infinite(data2_mod2$log10p)& (data2_mod2$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data3_mod2$log10p) & (data3_mod2$log10p > 0)), 45, data3_mod2$log10p) #Infの場合に50に置換したい
replace(data3_mod2$log10p, which(is.infinite(data3_mod2$log10p) & (data3_mod2$log10p > 0)), 50) #Infの場合に50に置換
data3_mod2$log10p[is.infinite(data3_mod2$log10p)& (data3_mod2$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data4_mod2$log10p) & (data4_mod2$log10p > 0)), 45, data4_mod2$log10p) #Infの場合に50に置換したい
replace(data4_mod2$log10p, which(is.infinite(data4_mod2$log10p) & (data4_mod2$log10p > 0)), 50) #Infの場合に50に置換
data4_mod2$log10p[is.infinite(data4_mod2$log10p)& (data4_mod2$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data5_mod2$log10p) & (data5_mod2$log10p > 0)), 45, data5_mod2$log10p) #Infの場合に50に置換したい
replace(data5_mod2$log10p, which(is.infinite(data5_mod2$log10p) & (data5_mod2$log10p > 0)), 50) #Infの場合に50に置換
data5_mod2$log10p[is.infinite(data5_mod2$log10p)& (data5_mod2$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data6_mod2$log10p) & (data6_mod2$log10p > 0)), 45, data6_mod2$log10p) #Infの場合に50に置換したい
replace(data6_mod2$log10p, which(is.infinite(data6_mod2$log10p) & (data6_mod2$log10p > 0)), 50) #Infの場合に50に置換
data6_mod2$log10p[is.infinite(data6_mod2$log10p)& (data6_mod2$log10p > 0)] <- 50 #置換の実行

#top50データの抽出
data1_top <- data1_mod2 %>% #data1_mod2のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data1_mod3 <- left_join(data1_mod2, data1_top, by = "GOID") #GOID識別子としてmerge
data2_top <- data2_mod2 %>% #data2_mod2のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data2_mod3 <- left_join(data2_mod2, data2_top, by = "GOID") #GOID識別子としてmerge
data3_top <- data3_mod2 %>% #data3_mod2のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data3_mod3 <- left_join(data3_mod2, data3_top, by = "GOID") #GOID識別子としてmerge
data4_top <- data4_mod2 %>% #data4_mod2のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data4_mod3 <- left_join(data4_mod2, data4_top, by = "GOID") #GOID識別子としてmerge
data5_top <- data5_mod2 %>% #data5_mod2のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data5_mod3 <- left_join(data5_mod2, data5_top, by = "GOID") #GOID識別子としてmerge
data6_top <- data6_mod2 %>% #data6_mod2のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data6_mod3 <- left_join(data6_mod2, data6_top, by = "GOID") #GOID識別子としてmerge

#NA(欠損値)の置換
ifelse(is.na(data1_mod3$log10p.y), 0, data1_mod3$log10p.y)
replace(data1_mod3$log10p.y, which(is.na(data1_mod3$log10p.y)), 0)
data1_mod3$log10p.y[is.na(data1_mod3$log10p.y)] <- 0
ifelse(is.na(data2_mod3$log10p.y), 0, data2_mod3$log10p.y)
replace(data2_mod3$log10p.y, which(is.na(data2_mod3$log10p.y)), 0)
data2_mod3$log10p.y[is.na(data2_mod3$log10p.y)] <- 0
ifelse(is.na(data3_mod3$log10p.y), 0, data3_mod3$log10p.y)
replace(data3_mod3$log10p.y, which(is.na(data3_mod3$log10p.y)), 0)
data3_mod3$log10p.y[is.na(data3_mod3$log10p.y)] <- 0
ifelse(is.na(data4_mod3$log10p.y), 0, data4_mod3$log10p.y)
replace(data4_mod3$log10p.y, which(is.na(data4_mod3$log10p.y)), 0)
data4_mod3$log10p.y[is.na(data4_mod3$log10p.y)] <- 0
ifelse(is.na(data5_mod3$log10p.y), 0, data5_mod3$log10p.y)
replace(data5_mod3$log10p.y, which(is.na(data5_mod3$log10p.y)), 0)
data5_mod3$log10p.y[is.na(data5_mod3$log10p.y)] <- 0
ifelse(is.na(data6_mod3$log10p.y), 0, data6_mod3$log10p.y)
replace(data6_mod3$log10p.y, which(is.na(data6_mod3$log10p.y)), 0)
data6_mod3$log10p.y[is.na(data6_mod3$log10p.y)] <- 0

#top50データをラベル
data1_mod4 <- mutate(data1_mod3, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data2_mod4 <- mutate(data2_mod3, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data3_mod4 <- mutate(data2_mod3, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data4_mod4 <- mutate(data2_mod3, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data5_mod4 <- mutate(data2_mod3, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data6_mod4 <- mutate(data2_mod3, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))

#ggplot2#############################################################
#volcano plot
plot1 <- ggplot(data1_mod4, aes(x = `Nr. Genes.x`, y = log10p.x, #xy値設定,`% Associated Genes.x`はやめる
                                shape = top50,
                                colour = top50, 
                                fill = `Nr. Genes.x`#カラー設定
)) +
  geom_point(size = data1_mod4$`Nr. Genes.x` * 0.02,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
plot1
vol1 <- plot1 +
  theme_bw()+
  scale_shape_manual(values = c(21,21))+  #プロット形状設定
  scale_colour_manual(values = c("black","red"))+  #プロット形状設定
  scale_fill_gradient2(
    low = "dark grey",
    mid = "black",
    high  = "dark red",
    midpoint = 1
  )+
  scale_x_continuous(limits = c(0, 200)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 50)) #x軸の範囲設定
vol1

plot2 <- ggplot(data2_mod4, aes(x = `Nr. Genes.x`, y = log10p.x, #xy値設定,`% Associated Genes.x`はやめる
                                shape = top50,
                                colour = top50, 
                                fill = `Nr. Genes.x`#カラー設定
)) +
  geom_point(size = data2_mod4$`Nr. Genes.x` * 0.2,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
vol2 <- plot2 +
  theme_bw()+
  scale_shape_manual(values = c(21,21))+  #プロット形状設定
  scale_colour_manual(values = c("black","red"))+  #プロット形状設定
  scale_fill_gradient2(
    low = "dark grey",
    mid = "black",
    high  = "dark red",
    midpoint = 1
  )+
  scale_x_continuous(limits = c(0, 30)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 15)) #x軸の範囲設定
vol2



plot3 <- ggplot(data3_mod4, aes(x = `Nr. Genes.x`, y = log10p.x, #xy値設定,`% Associated Genes.x`はやめる
                                shape = top50,
                                colour = top50, 
                                fill = `Nr. Genes.x`#カラー設定
)) +
  geom_point(size = data3_mod4$`Nr. Genes.x` * 0.2,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
plot3
vol3 <- plot3 +
  theme_bw()+
  scale_shape_manual(values = c(21,21))+  #プロット形状設定
  scale_colour_manual(values = c("black","red"))+  #プロット形状設定
  scale_fill_gradient2(
    low = "dark grey",
    mid = "black",
    high  = "dark red",
    midpoint = 1
  )+
  scale_x_continuous(limits = c(0, 30)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 15)) #x軸の範囲設定
vol3

plot4 <- ggplot(data4_mod4, aes(x = `Nr. Genes.x`, y = log10p.x, #xy値設定,`% Associated Genes.x`はやめる
                                shape = top50,
                                colour = top50, 
                                fill = `Nr. Genes.x`#カラー設定
)) +
  geom_point(size = data4_mod4$`Nr. Genes.x` * 0.2,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
vol4 <- plot4 +
  theme_bw()+
  scale_shape_manual(values = c(21,21))+  #プロット形状設定
  scale_colour_manual(values = c("black","red"))+  #プロット形状設定
  scale_fill_gradient2(
    low = "dark grey",
    mid = "black",
    high  = "dark red",
    midpoint = 1
  )+
  scale_x_continuous(limits = c(0, 30)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 15)) #x軸の範囲設定
vol4


plot5 <- ggplot(data5_mod4, aes(x = `Nr. Genes.x`, y = log10p.x, #xy値設定,`% Associated Genes.x`はやめる
                                shape = top50,
                                colour = top50, 
                                fill = `Nr. Genes.x`#カラー設定
)) +
  geom_point(size = data5_mod4$`Nr. Genes.x` * 0.2,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
plot5
vol5 <- plot5 +
  theme_bw()+
  scale_shape_manual(values = c(21,21))+  #プロット形状設定
  scale_colour_manual(values = c("black","red"))+  #プロット形状設定
  scale_fill_gradient2(
    low = "dark grey",
    mid = "black",
    high  = "dark red",
    midpoint = 1
  )+
  scale_x_continuous(limits = c(0, 30)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 15)) #x軸の範囲設定
vol5

plot6 <- ggplot(data6_mod4, aes(x = `Nr. Genes.x`, y = log10p.x, #xy値設定,`% Associated Genes.x`はやめる
                                shape = top50,
                                colour = top50, 
                                fill = `Nr. Genes.x`#カラー設定
)) +
  geom_point(size = data6_mod4$`Nr. Genes.x` * 0.2,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
vol6 <- plot6 +
  theme_bw()+
  scale_shape_manual(values = c(21,21))+  #プロット形状設定
  scale_colour_manual(values = c("black","red"))+  #プロット形状設定
  scale_fill_gradient2(
    low = "dark grey",
    mid = "black",
    high  = "dark red",
    midpoint = 1
  )+
  scale_x_continuous(limits = c(0, 30)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 15)) #x軸の範囲設定
vol6
#####################################################################

#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("ClueGOvol1.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol1) #プロット作成
dev.off() #PDF出力
pdf("ClueGOvol2.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol2) #プロット作成
dev.off() #PDF出力
pdf("ClueGOvol3.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol3) #プロット作成
dev.off() #PDF出力
pdf("ClueGOvol4.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol4) #プロット作成
dev.off() #PDF出力
pdf("ClueGOvol5.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol5) #プロット作成
dev.off() #PDF出力
pdf("ClueGOvol6.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol6) #プロット作成
dev.off() #PDF出力
############################################################################################
#プロット内容を作業フォルダにsvgで出力:dsvg() #作動しない
svg(file="ClueGOvol1.svg") #ファイル名指定
print(vol1) #プロット作成
dev.off() #svg出力
svg(file="ClueGOvol2.svg") #ファイル名指定
print(vol2) #プロット作成
dev.off() #svg出力
svg(file="ClueGOvol3.svg") #ファイル名指定
print(vol3) #プロット作成
dev.off() #svg出力
svg(file="ClueGOvol4.svg") #ファイル名指定
print(vol4) #プロット作成
dev.off() #svg出力
svg(file="ClueGOvol5.svg") #ファイル名指定
print(vol5) #プロット作成
dev.off() #svg出力
svg(file="ClueGOvol6.svg") #ファイル名指定
print(vol6) #プロット作成
dev.off() #svg出力
