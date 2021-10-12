#stringApp_01_GO_KEGG_volcanoplot.R
#stringApp1.5.1_Cytoscape3.8.0→stringApp1.6.0_Cytoscape3.8.1_Java11.0.6
#stringApp出力データから
#1_P,C,IntにおけるGO_KEGG_top50ラベル
#2_P,C,Intにおけるvolcano plot描写(top50赤色)
#3_P,C,Int全てについてGO_KEGG_top20抽出
#4_GO_KEGG termのplot描写
#5_GO_KEGG_top20に含まれるDEP list作成 →別script
#############################################################
#setwd("/Users/user/Dropbox/0_Work/R/SWATH") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY/2_Cytoscape") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP/2_Cytoscape") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc/2_Cytoscape") #作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC/2_Cytoscape") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR/2_Cytoscape") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
options(digits=2) #change digit2桁表示指定
##############################################################
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
#CSV入力
dat1 <- read.csv("String_twANOVA_Pq005_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
dat2 <- read.csv("String_twANOVA_Cq005_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
dat3 <- read.csv("String_twANOVA_PxCq005_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(dat1)
#############################################################
#行番号挿入,重複削除
dat1 <- dat1 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "PCP") %>% #シート番号挿入
  distinct(description,.keep_all=TRUE) %>% #重複削除
  mutate(log10p = -log10(`FDR.value`)) #p値→-log10p値変換
dat2 <- dat2 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "CLZ") %>% #シート番号挿入
  distinct(description,.keep_all=TRUE) %>% #重複削除
  mutate(log10p = -log10(`FDR.value`)) #p値→-log10p値変換
dat3 <- dat3 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "Interaction") %>% #シート番号挿入
  distinct(description,.keep_all=TRUE) %>% #重複削除
  mutate(log10p = -log10(`FDR.value`)) #p値→-log10p値変換
#############################################################
#Inf(無限大)の置換:無いので処理しない
#ifelse((is.infinite(dat1$log10p) & (dat1$log10p > 0)), 100, dat1$log10p) #Infの場合に100に置換したい
#replace(dat1$log10p, which(is.infinite(dat1$log10p) & (dat1$log10p > 0)), 100) #Infの場合に100に置換
#dat1$log10p[is.infinite(dat1$log10p)& (dat1$log10p > 0)] <- 100 #置換の実行
#dat1$log10p #置換結果の表示
#############################################################
#i列目の条件で抽出
#i <- 3 #3列目の条件
#varname <- colnames(dat1)[i] #条件抽出する列名varname
#x <- dat1 %>% dplyr::filter_(paste(varname, "==", 1)) #数字のみ抽出可能
#############################################################
#GO抽出
dat1_f <- dat1 %>% #dat1データの
  dplyr::filter(category=="GO Component" | #category列がGO Componentまたは
                  category=="GO Process" | #category列がGO Processまたは
                  category=="GO Function") #category列がGO Function
dat2_f <- dat2 %>% #dat2データの
  dplyr::filter(category=="GO Component" | #category列がGO Componentまたは
                  category=="GO Process" | #category列がGO Processまたは
                  category=="GO Function") #category列がGO Function
dat3_f <- dat3 %>% #dat3データの
  dplyr::filter(category=="GO Component" | #category列がGO Componentまたは
                  category=="GO Process" | #category列がGO Processまたは
                  category=="GO Function") #category列がGO Function
#KEGG抽出
dat4_f <- dat1 %>% #dat1データの
  dplyr::filter(category=="KEGG Pathways") #category列がKEGG Pathways
dat5_f <- dat2 %>% #dat2データの
  dplyr::filter(category=="KEGG Pathways") #category列がKEGG Pathways
dat6_f <- dat3 %>% #dat3データの
  dplyr::filter(category=="KEGG Pathways") #category列がKEGG Pathways
#############################################################
#top50データの抽出
dat1_top <- dat1_f %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat1_f <- left_join(dat1_f, dat1_top, by = "term.name") #term.name識別子としてmerge
dat2_top <- dat2_f %>% #P_GOのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat2_f <- left_join(dat2_f, dat2_top, by = "term.name") #term.name識別子としてmerge
dat3_top <- dat3_f %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat3_f <- left_join(dat3_f, dat3_top, by = "term.name") #term.name識別子としてmerge
dat4_top <- dat4_f %>% #P_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat4_f <- left_join(dat4_f, dat4_top, by = "term.name") #term.name識別子としてmerge
dat5_top <- dat5_f %>% #C_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat5_f <- left_join(dat5_f, dat5_top, by = "term.name") #term.name識別子としてmerge
dat6_top <- dat6_f %>% #PxC_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat6_f <- left_join(dat6_f, dat6_top, by = "term.name") #term.name識別子としてmerge
#############################################################
#NA(欠損値)の置換
ifelse(is.na(dat1_f$log10p.y), 0, dat1_f$log10p.y)
replace(dat1_f$log10p.y, which(is.na(dat1_f$log10p.y)), 0)
dat1_f$log10p.y[is.na(dat1_f$log10p.y)] <- 0
ifelse(is.na(dat2_f$log10p.y), 0, dat2_f$log10p.y)
replace(dat2_f$log10p.y, which(is.na(dat2_f$log10p.y)), 0)
dat2_f$log10p.y[is.na(dat2_f$log10p.y)] <- 0
ifelse(is.na(dat3_f$log10p.y), 0, dat3_f$log10p.y)
replace(dat3_f$log10p.y, which(is.na(dat3_f$log10p.y)), 0)
dat3_f$log10p.y[is.na(dat3_f$log10p.y)] <- 0
ifelse(is.na(dat4_f$log10p.y), 0, dat4_f$log10p.y)
replace(dat4_f$log10p.y, which(is.na(dat4_f$log10p.y)), 0)
dat4_f$log10p.y[is.na(dat4_f$log10p.y)] <- 0
ifelse(is.na(dat5_f$log10p.y), 0, dat5_f$log10p.y)
replace(dat5_f$log10p.y, which(is.na(dat5_f$log10p.y)), 0)
dat5_f$log10p.y[is.na(dat5_f$log10p.y)] <- 0
ifelse(is.na(dat6_f$log10p.y), 0, dat6_f$log10p.y)
replace(dat6_f$log10p.y, which(is.na(dat6_f$log10p.y)), 0)
dat6_f$log10p.y[is.na(dat6_f$log10p.y)] <- 0
#############################################################
#top50データをラベル
dat1_f <- mutate(dat1_f, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
dat2_f <- mutate(dat2_f, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
dat3_f <- mutate(dat3_f, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
dat4_f <- mutate(dat4_f, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
dat5_f <- mutate(dat5_f, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
dat6_f <- mutate(dat6_f, top50 = if_else(log10p.y > 0, true = TRUE, false = FALSE))
#############################################################
#dat1,2,3(GO)とdat4,5,6(KEGG)を1つのデータフレーム(GO)に結合
GO_tmp <- rbind(dat1_f, dat2_f) #dat1と2を結合
GO <- rbind(GO_tmp, dat3_f) #data1+2に3を結合
KEGG_tmp <- rbind(dat4_f, dat5_f) #dat4と5を結合
KEGG <- rbind(KEGG_tmp, dat6_f) #data4+5に6を結合
#############################################################
#重複削除,top20データの抽出
GO_mod1 <- GO %>%
  distinct(description.x,.keep_all=TRUE) #重複削除
GO_top <- GO_mod1 %>%
  arrange(desc(log10p.x)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p.x) #log10p上位20のデータを保存
KEGG_mod1 <- KEGG %>%
  distinct(description.x,.keep_all=TRUE) #重複削除
KEGG_top <- KEGG_mod1 %>%
  arrange(desc(log10p.x)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p.x) #log10p上位20のデータを保存
#############################################################
#top20の絞り込み
dat1_20 <- left_join(GO_top, dat1_f, by = "term.name") #term.name識別子としてmerge
dat2_20 <- left_join(GO_top, dat2_f, by = "term.name") #term.name識別子としてmerge
dat3_20 <- left_join(GO_top, dat3_f, by = "term.name") #term.name識別子としてmerge
dat4_20 <- left_join(KEGG_top, dat4_f, by = "term.name") #term.name識別子としてmerge
dat5_20 <- left_join(KEGG_top, dat5_f, by = "term.name") #term.name識別子としてmerge
dat6_20 <- left_join(KEGG_top, dat6_f, by = "term.name") #term.name識別子としてmerge
#############################################################
#グループ挿入
dat1_20 <- dat1_20 %>% mutate(group2 = "PCP") #グループ挿入
dat2_20 <- dat2_20 %>% mutate(group2 = "CLZ") #グループ挿入
dat3_20 <- dat3_20 %>% mutate(group2 = "Interaction") #グループ挿入
dat4_20 <- dat4_20 %>% mutate(group2 = "PCP") #グループ挿入
dat5_20 <- dat5_20 %>% mutate(group2 = "CLZ") #グループ挿入
dat6_20 <- dat6_20 %>% mutate(group2 = "Interaction") #グループ挿入
#############################################################
#dat1,2,3(GO)とdat4,5,6(KEGG)を1つのデータフレーム(GO)に再結合
GO_tmp <- rbind(dat1_20, dat2_20) #dat1と2を結合
GO_top <- rbind(GO_tmp, dat3_20) #dat1+2に3を結合
KEGG_tmp <- rbind(dat4_20, dat5_20) #dat4と5を結合
KEGG_top <- rbind(KEGG_tmp, dat6_20) #dat4+5に6を結合
#############################################################
#列名の置換
names(GO_top)[which(names(GO_top)=="log10p.y.y" ) ] <- "log10p"
names(KEGG_top)[which(names(KEGG_top)=="log10p.y.y" ) ] <- "log10p"
names(GO_top)[which(names(GO_top)=="X..genes.y.y" ) ] <- "Overlap"
names(KEGG_top)[which(names(KEGG_top)=="X..genes.y.y" ) ] <- "Overlap"
#############################################################
#NA(欠損値)の置換
ifelse(is.na(GO_top$log10p), 0, GO_top$log10p)
replace(GO_top$log10p, which(is.na(GO_top$log10p)), 0)
GO_top$log10p[is.na(GO_top$log10p)] <- 0
ifelse(is.na(KEGG_top$log10p), 0, KEGG_top$log10p)
replace(KEGG_top$log10p, which(is.na(KEGG_top$log10p)), 0)
KEGG_top$log10p[is.na(KEGG_top$log10p)] <- 0
#############################################################
#NA(欠損値)の置換
ifelse(is.na(GO_top$`Overlap`), 0, GO_top$`Overlap`)
replace(GO_top$`Overlap`, which(is.na(GO_top$`Overlap`)), 0)
GO_top$`Overlap`[is.na(GO_top$`Overlap`)] <- 0
ifelse(is.na(KEGG_top$`Overlap`), 0, KEGG_top$`Overlap`)
replace(KEGG_top$`Overlap`, which(is.na(KEGG_top$`Overlap`)), 0)
KEGG_top$`Overlap`[is.na(KEGG_top$`Overlap`)] <- 0
#############################################################
#description.x.xの順序の入れ替えのために
#データをファクタに変換する
#library(dplyr)
GO_top <- GO_top %>%
  mutate(description.x.x = as.factor(description.x.x)) #description.x.xをファクタに変換
KEGG_top <- KEGG_top %>%
  mutate(description.x.x = as.factor(description.x.x)) #description.x.xをファクタに変換
#############################################################
#列の作成と結合
GO_top <- GO_top %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category.x.x` == "GO Process" ~ "BP",      #category.x.xのGO ProcessをBP
    `category.x.x` == "GO Component" ~ "CC",    #category.x.xのGO ComponentをCC
    `category.x.x` == "GO Function" ~ "MP",     #category.x.xのGO FunctionをMP
    `category.x.x` == "KEGG Pathways" ~ "KEGG", #category.x.xのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category.x.x`)
  ))%>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description.x.x", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description.x.x_DBsourseとして結合
KEGG_top <- KEGG_top %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category.x.x` == "GO Process" ~ "BP",      #category.x.xのGO ProcessをBP
    `category.x.x` == "GO Component" ~ "CC",    #category.x.xのGO ComponentをCC
    `category.x.x` == "GO Function" ~ "MP",     #category.x.xのGO FunctionをMP
    `category.x.x` == "KEGG Pathways" ~ "KEGG", #category.x.xのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category.x.x`)
  ))%>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description.x.x", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description.x.x_DBsourseとして結合
#############################################################
#ggplot2
plot1 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = GO_top, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group2,
                                 alpha = 0.9,
                                 size = `Overlap`))
plot1

plot2 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = KEGG_top, aes(x = log10p,
                                  y = reorder(x = term, X = log10p),
                                   #position = "jitter", #geom_jitter
                                   color = group2,
                                   alpha = 0.9,
                                   size = `Overlap`))
plot2
############################################################################################
#svg出力
svg(file="stringAppGO.svg") #ファイル名指定
print(plot1) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG.svg") #ファイル名指定
print(plot2) #プロット作成
dev.off() #svg出力
#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("stringAppGO.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot1) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot2) #プロット作成
dev.off() #PDF出力
#############################################################
#ggplot2
#volcano plot
plot1 <- ggplot(dat1_f, aes(x = `X..genes.x`, y = log10p.x, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `X..genes.x` #カラー設定
)) +
  geom_point(size = dat1_f$`X..genes.x` * 0.02,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
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
  scale_x_continuous(limits = c(0, 250)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 55))    #y軸の範囲設定
vol1

plot2 <- ggplot(dat2_f, aes(x = `X..genes.x`, y = log10p.x, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `X..genes.x` #カラー設定
)) +
  geom_point(size = dat2_f$`X..genes.x` * 0.02,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 350)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 80))    #y軸の範囲設定
vol2

plot3 <- ggplot(dat3_f, aes(x = `X..genes.x`, y = log10p.x, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `X..genes.x` #カラー設定
)) +
  geom_point(size = dat3_f$`X..genes.x` * 0.02,  #プロットサイズ設定
             alpha = 0.8)  #プロット透明度設定
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
  scale_x_continuous(limits = c(0, 200)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 80))    #y軸の範囲設定
vol3

plot4 <- ggplot(dat4_f, aes(x = `X..genes.x`, y = log10p.x, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `X..genes.x` #カラー設定
)) +
  geom_point(size = dat4_f$`X..genes.x` * 0.1,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 30)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 10))    #y軸の範囲設定
vol4

plot5 <- ggplot(dat5_f, aes(x = `X..genes.x`, y = log10p.x, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `X..genes.x` #カラー設定
)) +
  geom_point(size = dat5_f$`X..genes.x` * 0.1,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 70)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 20))    #y軸の範囲設定
vol5

plot6 <- ggplot(dat6_f, aes(x = `X..genes.x`, y = log10p.x, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `X..genes.x` #カラー設定
)) +
  geom_point(size = dat6_f$`X..genes.x` * 0.1,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 35)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 20))    #y軸の範囲設定
vol6
#############################################################
#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("stringApp_vol1.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol1) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol2.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol2) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol3.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol3) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol4.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol4) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol5.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol5) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol6.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol6) #プロット作成
dev.off() #PDF出力
#############################################################
#xslx出力
smp <- list("GO"=GO_top, "KEGG"=KEGG_top) #GO,KEGGのtop20リスト作成
write.xlsx(smp, "GO_KEGG.xlsx") #GOシート,KEGGシート出力
#エクセルでRの特殊記号"|"を","に置換し、GO_KEGGrep.xlsxとして保存！！！
q()
