#WebGestalt_01_volcanoplot.R
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/SWATH")
getwd() #作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
library(EnhancedVolcano)
library(magrittr)
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(dplyr)
library(scales) #muted()関数使用のため
library(gridExtra) #svg出力のため
#############################################################
#Excelファイル読み込み
#install.packages("readxl")
library(readxl)
data1 <- read_excel("PC10_PFC_pathway.xlsx", 2) #シート2の読み込み
data2 <- read_excel("PC10_PFC_pathway.xlsx", sheet = 3) #シート3の読み込み
data3 <- read_excel("PC10_PFC_pathway.xlsx", sheet = 4) #シート4の読み込み
data4 <- read_excel("PC10_PFC_pathway.xlsx", sheet = 5) #シート5の読み込み
data5 <- read_excel("PC10_PFC_pathway.xlsx", sheet = 6) #シート6の読み込み
data6 <- read_excel("PC10_PFC_pathway.xlsx", sheet = 7) #シート7の読み込み
#data4 <- read_excel("Stat.xlsx", sheet = "datafile3") #シート名で読み込み
#data5 <- read_excel("Stat.xlsx", sheet = "datafile4", col_names = TRUE, col_types = c("blank", "text", "date", "text", "numeric"))  #最初の列を削除(blank)し、後続3列の型を指定

#Log2ER列の作成
data1_mod <- mutate(data1, Log2ER = log2(data1$enrichmentRatio))
data2_mod <- mutate(data2, Log2ER = log2(data2$enrichmentRatio))
data3_mod <- mutate(data3, Log2ER = log2(data3$enrichmentRatio))
data4_mod <- mutate(data4, Log2ER = log2(data4$enrichmentRatio))
data5_mod <- mutate(data5, Log2ER = log2(data5$enrichmentRatio))
data6_mod <- mutate(data6, Log2ER = log2(data6$enrichmentRatio))

#-Log10FDRの作成
data1_mod2 <- mutate(data1_mod, Log10FDR = -log10(data1$FDR))
data2_mod2 <- mutate(data2_mod, Log10FDR = -log10(data2$FDR))
data3_mod2 <- mutate(data3_mod, Log10FDR = -log10(data3$FDR))
data4_mod2 <- mutate(data4_mod, Log10FDR = -log10(data4$FDR))
data5_mod2 <- mutate(data5_mod, Log10FDR = -log10(data5$FDR))
data6_mod2 <- mutate(data6_mod, Log10FDR = -log10(data6$FDR))

#Inf(無限大)の置換,WebGestalt2019では-Log10FDRが15.65で頭打ち
ifelse((is.infinite(data1_mod2$Log10FDR) & (data1_mod2$Log10FDR > 0)), 15.65, data1_mod2$Log10FDR) #Infの場合に15.65に置換したい
replace(data1_mod2$Log10FDR, which(is.infinite(data1_mod2$Log10FDR) & (data1_mod2$Log10FDR > 0)), 15.65) #Infの場合に15.65に置換
data1_mod2$Log10FDR[is.infinite(data1_mod2$Log10FDR)& (data1_mod2$Log10FDR > 0)] <- 15.65 #置換の実行
data1_mod2$Log10FDR #置換結果の表示
ifelse((is.infinite(data2_mod2$Log10FDR) & (data2_mod2$Log10FDR > 0)), 15.65, data2_mod2$Log10FDR) #Infの場合に15.65に置換したい
replace(data2_mod2$Log10FDR, which(is.infinite(data2_mod2$Log10FDR) & (data2_mod2$Log10FDR > 0)), 15.65) #Infの場合に15.65に置換
data2_mod2$Log10FDR[is.infinite(data2_mod2$Log10FDR)& (data2_mod2$Log10FDR > 0)] <- 15.65 #置換の実行
data2_mod2$Log10FDR #置換結果の表示
ifelse((is.infinite(data3_mod2$Log10FDR) & (data3_mod2$Log10FDR > 0)), 15.65, data3_mod2$Log10FDR) #Infの場合に15.65に置換したい
replace(data3_mod2$Log10FDR, which(is.infinite(data3_mod2$Log10FDR) & (data3_mod2$Log10FDR > 0)), 15.65) #Infの場合に15.65に置換
data3_mod2$Log10FDR[is.infinite(data3_mod2$Log10FDR)& (data3_mod2$Log10FDR > 0)] <- 15.65 #置換の実行
data3_mod2$Log10FDR #置換結果の表示
ifelse((is.infinite(data4_mod2$Log10FDR) & (data4_mod2$Log10FDR > 0)), 15.65, data4_mod2$Log10FDR) #Infの場合に15.65に置換したい
replace(data4_mod2$Log10FDR, which(is.infinite(data4_mod2$Log10FDR) & (data4_mod2$Log10FDR > 0)), 15.65) #Infの場合に15.65に置換
data4_mod2$Log10FDR[is.infinite(data4_mod2$Log10FDR)& (data4_mod2$Log10FDR > 0)] <- 15.65 #置換の実行
data4_mod2$Log10FDR #置換結果の表示
ifelse((is.infinite(data5_mod2$Log10FDR) & (data5_mod2$Log10FDR > 0)), 15.65, data5_mod2$Log10FDR) #Infの場合に15.65に置換したい
replace(data5_mod2$Log10FDR, which(is.infinite(data5_mod2$Log10FDR) & (data5_mod2$Log10FDR > 0)), 15.65) #Infの場合に15.65に置換
data5_mod2$Log10FDR[is.infinite(data5_mod2$Log10FDR)& (data5_mod2$Log10FDR > 0)] <- 15.65 #置換の実行
data5_mod2$Log10FDR #置換結果の表示
ifelse((is.infinite(data6_mod2$Log10FDR) & (data6_mod2$Log10FDR > 0)), 15.65, data6_mod2$Log10FDR) #Infの場合に15.65に置換したい
replace(data6_mod2$Log10FDR, which(is.infinite(data6_mod2$Log10FDR) & (data6_mod2$Log10FDR > 0)), 15.65) #Infの場合に15.65に置換
data6_mod2$Log10FDR[is.infinite(data6_mod2$Log10FDR)& (data6_mod2$Log10FDR > 0)] <- 15.65 #置換の実行
data6_mod2$Log10FDR #置換結果の表示

#top50データの抽出
data4_top50 <- data4_mod2 %>% #data4_mod2のデータから
  arrange(desc(Log10FDR)) %>% #Log10FDRで降順ソート
  top_n(n = 50, wt = Log10FDR) #Log10FDR上位50のデータを保存
data4_top50
data4_mod3 <- left_join(data4_mod2, data4_top50, by = "geneSet")

data5_top50 <- data5_mod2 %>% #data5_mod2のデータから
  arrange(desc(Log10FDR)) %>% #Log10FDRで降順ソート
  top_n(n = 50, wt = Log10FDR) #Log10FDR上位50のデータを保存
data5_top50
data5_mod3 <- left_join(data5_mod2, data5_top50, by = "geneSet")

data6_top50 <- data6_mod2 %>% #data6_mod2のデータから
  arrange(desc(Log10FDR)) %>% #Log10FDRで降順ソート
  top_n(n = 50, wt = Log10FDR) #Log10FDR上位50のデータを保存
data6_top50
data6_mod3 <- left_join(data6_mod2, data6_top50, by = "geneSet")

#NA(欠損値)の置換
ifelse(is.na(data4_mod3$Log10FDR.y), 0, data4_mod3$Log10FDR.y)
replace(data4_mod3$Log10FDR.y, which(is.na(data4_mod3$Log10FDR.y)), 0)
data4_mod3$Log10FDR.y[is.na(data4_mod3$Log10FDR.y)] <- 0

ifelse(is.na(data5_mod3$Log10FDR.y), 0, data5_mod3$Log10FDR.y)
replace(data5_mod3$Log10FDR.y, which(is.na(data5_mod3$Log10FDR.y)), 0)
data5_mod3$Log10FDR.y[is.na(data5_mod3$Log10FDR.y)] <- 0

ifelse(is.na(data6_mod3$Log10FDR.y), 0, data6_mod3$Log10FDR.y)
replace(data6_mod3$Log10FDR.y, which(is.na(data6_mod3$Log10FDR.y)), 0)
data6_mod3$Log10FDR.y[is.na(data6_mod3$Log10FDR.y)] <- 0

#top50データをラベル
data4_mod4 <- mutate(data4_mod3, top50 = if_else(Log10FDR.y > 0, true = TRUE, false = FALSE))
data5_mod4 <- mutate(data5_mod3, top50 = if_else(Log10FDR.y > 0, true = TRUE, false = FALSE))
data6_mod4 <- mutate(data6_mod3, top50 = if_else(Log10FDR.y > 0, true = TRUE, false = FALSE))



#ggplot2#############################################################
plot1 <- ggplot(data1_mod2, aes(x = Log2ER, y = Log10FDR, #xy値設定
                                colour = overlap, #カラー設定
)) +
  geom_point(size = data1_mod2$overlap * 0.1,  #プロットサイズ設定
             #shape = 21,  #プロット形状設定
             alpha = 0.8) + #プロット透明度設定
  scale_x_continuous(limits = c(0, 5)) #x軸の範囲設定
#+ scale_y_reverse() #y軸の逆転設定
vol1 <- plot1 +
  theme_bw()+
  scale_colour_gradient2(
    low = "white",
    mid = "dark grey",
    high  = "black",
    midpoint = 1
  )#+
vol1
#  scale_x_continuous(limits = c(0, 5))
#  + scale_y_reverse()


plot2 <- ggplot(data2_mod2, aes(x = Log2ER, y = Log10FDR, #xy値設定
                                colour = overlap, #カラー設定
)) +
  geom_point(size = data2_mod2$overlap * 0.1,  #プロットサイズ設定
             #shape = 21,  #プロット形状設定
             alpha = 0.8) + #プロット透明度設定
  scale_x_continuous(limits = c(0, 5)) #x軸の範囲設定
vol2 <- plot2 +
  theme_bw()+
  scale_colour_gradient2(
    low = "white",
    mid = "dark grey",
    high  = "black",
    midpoint = 1
  )
vol2


plot3 <- ggplot(data3_mod2, aes(x = Log2ER, y = Log10FDR, #xy値設定
                                colour = overlap, #カラー設定
)) +
  geom_point(size = data3_mod2$overlap * 0.1,  #プロットサイズ設定
             #shape = 21,  #プロット形状設定
             alpha = 0.8) + #プロット透明度設定
  scale_x_continuous(limits = c(0, 5)) #x軸の範囲設定
vol3 <- plot3 +
  theme_bw()+
  scale_colour_gradient2(
    low = "white",
    mid = "dark grey",
    high  = "black",
    midpoint = 1
  )
vol3



plot4 <- ggplot(data4_mod4, aes(x = Log2ER.x, y = Log10FDR.x, #xy値設定
                                shape = top50,
                                colour = top50, 
                                fill = overlap.x#カラー設定
)) +
  geom_point(size = data4_mod4$overlap.x * 0.1,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(-4, 5)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 16)) #x軸の範囲設定
vol4


plot5 <- ggplot(data5_mod4, aes(x = Log2ER.x, y = Log10FDR.x, #xy値設定
                                shape = top50,
                                colour = top50, 
                                fill = overlap.x#カラー設定
)) +
  geom_point(size = data5_mod4$overlap.x * 0.1,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
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
  scale_x_continuous(limits = c(-4, 5)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 16)) #x軸の範囲設定
vol5



plot6 <- ggplot(data6_mod4, aes(x = Log2ER.x, y = Log10FDR.x, #xy値設定
                                shape = top50,
                                colour = top50, 
                                fill = overlap.x#カラー設定
)) +
  geom_point(size = data6_mod4$overlap.x * 0.1,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(-4, 5)) +#x軸の範囲設定
  scale_y_continuous(limits = c(0, 16)) #x軸の範囲設定
vol6
#####################################################################

#P359################################################################
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
pdf("vol6.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
print(vol6) #プロット作成
dev.off() #PDF出力
############################################################################################
#プロット内容を作業フォルダにsvgで出力:dsvg() #作動しない
#dsvg(file = "vol6.svg", #ファイル名指定
#     width = 6, height = 6, #プロットサイズをインチで指定
#     bg = "white", #背景色
#     pointsize = 12, 
#     standalone = TRUE, #ファイルのみでプロットを表示
#     fontname_serif = "Times New Roman", #文字種類を指定
#     fontname_sans = "Calibri", 
#     fontname_mono = "Courier New", 
#     fontname_symbol = "Symbol")
#print(vol6) #プロット作成
#dev.off() #svg出力

svg(file="vol1.svg") #ファイル名指定
print(vol1) #プロット作成
dev.off() #svg出力
svg(file="vol2.svg") #ファイル名指定
print(vol2) #プロット作成
dev.off() #svg出力
svg(file="vol3.svg") #ファイル名指定
print(vol3) #プロット作成
dev.off() #svg出力
svg(file="vol4.svg") #ファイル名指定
print(vol4) #プロット作成
dev.off() #svg出力
svg(file="vol5.svg") #ファイル名指定
print(vol5) #プロット作成
dev.off() #svg出力
svg(file="vol6.svg") #ファイル名指定
print(vol6) #プロット作成
dev.off() #svg出力

#docxまたはpptxにプロットをsvgで埋め込み出力
#複雑なプロットはファイルエラーが出るので注意
#docxファイル:write_docxコマンド

#library(rvg) 
#write_docx(file = "my_plot.docx", code = vol(1:6)) #作動しない
#pptxファイル:write_pptxコマンド
#write_pptx(file = "my_plot.pptx", code = vol(1:6)) #作動しない
############################################################################################









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