#ClueGO_02_GO_KEGG_list_preparation1.R
#縦(行)にデータを追加
#https://bioscryptome.t-ohashi.info/r/dataframe-bind/#データフレームに行を追加する縦に結合する
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/SWATH")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
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
data2 <- read_excel("ClueGOResultTable_PC10_PFC.xls", sheet = 2) #シート2の読み込み
data3 <- read_excel("ClueGOResultTable_PC10_PFC.xls", sheet = 3) #シート3の読み込み
data4 <- read_excel("ClueGOResultTable_PC10_PFC.xls", sheet = 4) #シート4の読み込み
data5 <- read_excel("ClueGOResultTable_PC10_PFC.xls", sheet = 5) #シート5の読み込み
data6 <- read_excel("ClueGOResultTable_PC10_PFC.xls", sheet = 6) #シート6の読み込み
#data3 <- read_excel("Stat.xlsx", sheet = "datafile3") #シート名で読み込み
#data4 <- read_excel("Stat.xlsx", sheet = "datafile4", col_names = TRUE, col_types = c("blank", "text", "date", "text", "numeric"))  #最初の列を削除(blank)し、後続3列の型を指定
#str(data1) #データ構造確認

#行番号挿入,シート番号挿入,重複削除,-log10p値変換
data1_mod1 <- data1 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "sheet1")%>% #シート番号挿入
  distinct(GOTerm,.keep_all=TRUE)%>% #重複削除
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))#p値→-log10p値変換
#write.xlsx(data1_mod1, "out1.xlsx", sheetname="sheet1") #出力
#str(data1_mod1)
data2_mod1 <- data2 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "sheet2")%>% #シート番号挿入
  distinct(GOTerm,.keep_all=TRUE)%>% #重複削除
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))#p値→-log10p値変換
data3_mod1 <- data3 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "sheet3")%>% #シート番号挿入
  distinct(GOTerm,.keep_all=TRUE)%>% #重複削除
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))#p値→-log10p値変換
data4_mod1 <- data4 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "sheet4")%>% #シート番号挿入
  distinct(GOTerm,.keep_all=TRUE)%>% #重複削除
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))#p値→-log10p値変換
data5_mod1 <- data5 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "sheet5")%>% #シート番号挿入
  distinct(GOTerm,.keep_all=TRUE)%>% #重複削除
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))#p値→-log10p値変換
data6_mod1 <- data6 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "sheet6")%>% #シート番号挿入
  distinct(GOTerm,.keep_all=TRUE)%>% #重複削除
  mutate(log10p = -log10(`Term PValue Corrected with Bonferroni step down`))#p値→-log10p値変換

#Inf(無限大)の置換,ClueGOでは-Log10FDRが50くらいで頭打ち?
ifelse((is.infinite(data1_mod1$log10p) & (data1_mod1$log10p > 0)), 50, data1_mod1$log10p) #Infの場合に50に置換したい
replace(data1_mod1$log10p, which(is.infinite(data1_mod1$log10p) & (data1_mod1$log10p > 0)), 50) #Infの場合に50に置換
data1_mod1$log10p[is.infinite(data1_mod1$log10p) & (data1_mod1$log10p > 0)] <- 50 #置換の実行
#data1_mod1$log10p #置換結果の表示
ifelse((is.infinite(data2_mod1$log10p) & (data2_mod1$log10p > 0)), 50, data2_mod1$log10p) #Infの場合に50に置換したい
replace(data2_mod1$log10p, which(is.infinite(data2_mod1$log10p) & (data2_mod1$log10p > 0)), 50) #Infの場合に50に置換
data2_mod1$log10p[is.infinite(data2_mod1$log10p) & (data2_mod1$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data3_mod1$log10p) & (data3_mod1$log10p > 0)), 50, data3_mod1$log10p) #Infの場合に50に置換したい
replace(data3_mod1$log10p, which(is.infinite(data3_mod1$log10p) & (data3_mod1$log10p > 0)), 50) #Infの場合に50に置換
data3_mod1$log10p[is.infinite(data3_mod1$log10p) & (data3_mod1$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data4_mod1$log10p) & (data4_mod1$log10p > 0)), 50, data4_mod1$log10p) #Infの場合に50に置換したい
replace(data4_mod1$log10p, which(is.infinite(data4_mod1$log10p) & (data4_mod1$log10p > 0)), 50) #Infの場合に50に置換
data4_mod1$log10p[is.infinite(data4_mod1$log10p) & (data4_mod1$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data5_mod1$log10p) & (data5_mod1$log10p > 0)), 50, data5_mod1$log10p) #Infの場合に50に置換したい
replace(data5_mod1$log10p, which(is.infinite(data5_mod1$log10p) & (data5_mod1$log10p > 0)), 50) #Infの場合に50に置換
data5_mod1$log10p[is.infinite(data5_mod1$log10p) & (data5_mod1$log10p > 0)] <- 50 #置換の実行
ifelse((is.infinite(data6_mod1$log10p) & (data6_mod1$log10p > 0)), 50, data6_mod1$log10p) #Infの場合に50に置換したい
replace(data6_mod1$log10p, which(is.infinite(data6_mod1$log10p) & (data6_mod1$log10p > 0)), 50) #Infの場合に50に置換
data6_mod1$log10p[is.infinite(data6_mod1$log10p) & (data6_mod1$log10p > 0)] <- 50 #置換の実行

#top50データの抽出
data1_top <- data1_mod1 %>% #data1_mod1のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data1_mod2 <- left_join(data1_mod1, data1_top, by = "GOID") #GOID識別子としてmerge
data2_top <- data2_mod1 %>% #data1_mod1のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data2_mod2 <- left_join(data2_mod1, data2_top, by = "GOID") #GOID識別子としてmerge
data3_top <- data3_mod1 %>% #data1_mod1のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data3_mod2 <- left_join(data3_mod1, data3_top, by = "GOID") #GOID識別子としてmerge
data4_top <- data4_mod1 %>% #data1_mod1のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data4_mod2 <- left_join(data4_mod1, data4_top, by = "GOID") #GOID識別子としてmerge
data5_top <- data5_mod1 %>% #data1_mod1のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data5_mod2 <- left_join(data5_mod1, data5_top, by = "GOID") #GOID識別子としてmerge
data6_top <- data6_mod1 %>% #data1_mod1のデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
data6_mod2 <- left_join(data6_mod1, data6_top, by = "GOID") #GOID識別子としてmerge

#NA(欠損値)の置換
ifelse(is.na(data1_mod2$log10p.y), 0, data1_mod2$log10p.y)
replace(data1_mod2$log10p.y, which(is.na(data1_mod2$log10p.y)), 0)
data1_mod2$log10p.y[is.na(data1_mod2$log10p.y)] <- 0
ifelse(is.na(data2_mod2$log10p.y), 0, data2_mod2$log10p.y)
replace(data2_mod2$log10p.y, which(is.na(data2_mod2$log10p.y)), 0)
data2_mod2$log10p.y[is.na(data2_mod2$log10p.y)] <- 0
ifelse(is.na(data3_mod2$log10p.y), 0, data3_mod2$log10p.y)
replace(data3_mod2$log10p.y, which(is.na(data3_mod2$log10p.y)), 0)
data3_mod2$log10p.y[is.na(data3_mod2$log10p.y)] <- 0
ifelse(is.na(data4_mod2$log10p.y), 0, data4_mod2$log10p.y)
replace(data4_mod2$log10p.y, which(is.na(data4_mod2$log10p.y)), 0)
data4_mod2$log10p.y[is.na(data4_mod2$log10p.y)] <- 0
ifelse(is.na(data5_mod2$log10p.y), 0, data5_mod2$log10p.y)
replace(data5_mod2$log10p.y, which(is.na(data5_mod2$log10p.y)), 0)
data5_mod2$log10p.y[is.na(data5_mod2$log10p.y)] <- 0
ifelse(is.na(data6_mod2$log10p.y), 0, data6_mod2$log10p.y)
replace(data6_mod2$log10p.y, which(is.na(data6_mod2$log10p.y)), 0)
data6_mod2$log10p.y[is.na(data6_mod2$log10p.y)] <- 0

#top50データをラベル
data1_mod3 <- mutate(data1_mod2, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data2_mod3 <- mutate(data2_mod2, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data3_mod3 <- mutate(data3_mod2, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data4_mod3 <- mutate(data4_mod2, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data5_mod3 <- mutate(data5_mod2, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
data6_mod3 <- mutate(data6_mod2, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))

#ここまではvolcanoplotと同じ？#############################################################################
#data1,3,5(GO)とdata2,4,6(KEGG)を1つのデータフレーム(GO)に結合
GO_tmp <- rbind(data1_mod3, data3_mod3) #data1と3を結合
GO <- rbind(GO_tmp, data5_mod3) #data1+3に5を結合
KEGG_tmp <- rbind(data2_mod3, data4_mod3) #data2と4を結合
KEGG <- rbind(KEGG_tmp, data6_mod3) #data2+4に6を結合

#重複削除,top20データの抽出
GO_mod1 <- GO %>%
  distinct(GOTerm.x,.keep_all=TRUE) #重複削除
GO_top <- GO_mod1 %>%
  arrange(desc(log10p.x)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p.x) #log10p上位20のデータを保存
KEGG_mod1 <- KEGG %>%
  distinct(GOTerm.x,.keep_all=TRUE) #重複削除
KEGG_top <- KEGG_mod1 %>%
  arrange(desc(log10p.x)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p.x) #log10p上位20のデータを保存

#top20の絞り込み
data1_top2 <- left_join(GO_top, data1_mod1, by = "GOID") #GOID識別子としてmerge
data3_top2 <- left_join(GO_top, data3_mod1, by = "GOID") #GOID識別子としてmerge
data5_top2 <- left_join(GO_top, data5_mod1, by = "GOID") #GOID識別子としてmerge
data2_top2 <- left_join(KEGG_top, data2_mod1, by = "GOID") #GOID識別子としてmerge
data4_top2 <- left_join(KEGG_top, data4_mod1, by = "GOID") #GOID識別子としてmerge
data6_top2 <- left_join(KEGG_top, data6_mod1, by = "GOID") #GOID識別子としてmerge

#グループ挿入
data1_top3 <- data1_top2 %>% mutate(group2 = "PCP") #グループ挿入
data2_top3 <- data2_top2 %>% mutate(group2 = "PCP") #グループ挿入
data3_top3 <- data3_top2 %>% mutate(group2 = "CLZ") #グループ挿入  
data4_top3 <- data4_top2 %>% mutate(group2 = "CLZ") #グループ挿入
data5_top3 <- data5_top2 %>% mutate(group2 = "Interaction") #グループ挿入
data6_top3 <- data6_top2 %>% mutate(group2 = "Interaction") #グループ挿入

#data1,3,5(GO)とdata2,4,6(KEGG)を1つのデータフレーム(GO)に再結合
GO_tmp2 <- rbind(data1_top3, data3_top3) #data1と3を結合
GO_top2 <- rbind(GO_tmp2, data5_top3) #data1+3に5を結合
KEGG_tmp2 <- rbind(data2_top3, data4_top3) #data2と4を結合
KEGG_top2 <- rbind(KEGG_tmp2, data6_top3) #data2+4に6を結合

#NA(欠損値)の置換
ifelse(is.na(GO_top2$log10p), 0, GO_top2$log10p)
replace(GO_top2$log10p, which(is.na(GO_top2$log10p)), 0)
GO_top2$log10p[is.na(GO_top2$log10p)] <- 0
ifelse(is.na(KEGG_top2$log10p), 0, KEGG_top2$log10p)
replace(KEGG_top2$log10p, which(is.na(KEGG_top2$log10p)), 0)
KEGG_top2$log10p[is.na(KEGG_top2$log10p)] <- 0

#NA(欠損値)の置換
ifelse(is.na(GO_top2$`Nr. Genes`), 0, GO_top2$`Nr. Genes`)
replace(GO_top2$`Nr. Genes`, which(is.na(GO_top2$`Nr. Genes`)), 0)
GO_top2$`Nr. Genes`[is.na(GO_top2$`Nr. Genes`)] <- 0
ifelse(is.na(KEGG_top2$`Nr. Genes`), 0, KEGG_top2$`Nr. Genes`)
replace(KEGG_top2$`Nr. Genes`, which(is.na(KEGG_top2$`Nr. Genes`)), 0)
KEGG_top2$`Nr. Genes`[is.na(KEGG_top2$`Nr. Genes`)] <- 0

#GOTerm.xの順序の入れ替えのために
#データをファクタに変換する
#library(dplyr)
GO_top2 <- GO_top2 %>%
  mutate(GOTerm.x = as.factor(GOTerm.x)) #GOTerm.xをファクタに変換
levels(GO_top2$GOTerm.x)
str(GO_top2)
KEGG_top2 <- KEGG_top2 %>%
  mutate(GOTerm.x = as.factor(GOTerm.x)) #GOTerm.xをファクタに変換
levels(KEGG_top2$GOTerm.x)
str(KEGG_top2)



#ggplot2#############################################################
plot1 <- ggplot()+
  theme_gray(base_family = "Arial")+ #フォント指定"HiraKakuPro-W3"など
  geom_point(data = GO_top2, aes(x = log10p, 
                                 y = reorder(x = GOTerm.x, X = log10p), 
                                  #position = "jitter", #geom_jitter
                                  color = group2,
                                  alpha = 0.9,
                                  size = `Nr. Genes`)
  )
plot1

plot2 <- ggplot()+
  theme_gray(base_family = "Arial")+ #フォント指定"HiraKakuPro-W3"など
  geom_point(data = KEGG_top2, aes(x = log10p, 
                                 y = reorder(x = GOTerm.x, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group2,
                                 alpha = 0.9,
                                 size = `Nr. Genes`)
  )
plot2
#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
#上手く行かない
#pdf("ClueGOplotGO.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
#print(plot1) #プロット作成
#dev.off() #PDF出力
#pdf("ClueGOplotKEGG.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅と高さはインチ指定2.54で除してcm変換
#print(plot2) #プロット作成
#dev.off() #PDF出力
############################################################################################
#プロット内容を作業フォルダにsvgで出力:dsvg() #作動しない
svg(file="ClueGOplotGO.svg") #ファイル名指定
print(plot1) #プロット作成
dev.off() #svg出力
svg(file="ClueGOplotKEGG.svg") #ファイル名指定
print(plot2) #プロット作成
dev.off() #svg出力

