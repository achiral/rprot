#stringApp_01_GO_KEGG_volcanoplot.R
#stringApp1.5.1_Cytoscape3.8.0→stringApp1.6.0_Cytoscape3.8.1_Java11.0.6
#stringApp出力データから
#1_P,C,IntにおけるGO_KEGG_top50ラベル
#2_P,C,Int全てについてGO_KEGG_top50抽出
#3_GO_KEGG termのplot描写
#4_P,C,Intにおけるvolcano plot描写(top50赤色)
#5_GO_KEGG_top50に含まれるDEP list作成 →別script
#############################################################
rm(list = ls(all.names = TRUE))
#setwd("/Users/user/Dropbox/0_Work/R/SWATH") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR/stringApp") #作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC/stringApp") #作業ディレクトリ設定
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
#dat1 <- read.csv("Untitled_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(dat1)
#############################################################
#行番号挿入,重複削除,列名置換
dat1 <- dat1 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "PCP") %>% #シート番号挿入
  distinct(description,.keep_all=TRUE) %>% #重複削除
  mutate(log10p = -log10(`FDR.value`)) #p値→-log10p値変換
names(dat1)[which(names(dat1)=="X..genes" ) ] <- "Overlap" #列名置換
dat2 <- dat2 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "CLZ") %>% #シート番号挿入
  distinct(description,.keep_all=TRUE) %>% #重複削除
  mutate(log10p = -log10(`FDR.value`)) #p値→-log10p値変換
names(dat2)[which(names(dat2)=="X..genes" ) ] <- "Overlap" #列名置換
dat3 <- dat3 %>% mutate(No = row_number()) %>% #行番号挿入
  mutate(group = "Interaction") %>% #シート番号挿入
  distinct(description,.keep_all=TRUE) %>% #重複削除
  mutate(log10p = -log10(`FDR.value`)) #p値→-log10p値変換
names(dat3)[which(names(dat3)=="X..genes" ) ] <- "Overlap" #列名置換
#############################################################
#列の作成と結合
dat1 <- dat1 %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category` == "GO Process" ~ "BP",      #categoryのGO ProcessをBP
    `category` == "GO Component" ~ "CC",    #categoryのGO ComponentをCC
    `category` == "GO Function" ~ "MF",     #categoryのGO FunctionをMF
    `category` == "KEGG Pathways" ~ "KEGG", #categoryのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category`)
  )) %>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description_DBsourseとして結合
dat2 <- dat2 %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category` == "GO Process" ~ "BP",      #categoryのGO ProcessをBP
    `category` == "GO Component" ~ "CC",    #categoryのGO ComponentをCC
    `category` == "GO Function" ~ "MF",     #categoryのGO FunctionをMF
    `category` == "KEGG Pathways" ~ "KEGG", #categoryのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category`)
  )) %>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description_DBsourseとして結合
dat3 <- dat3 %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category` == "GO Process" ~ "BP",      #categoryのGO ProcessをBP
    `category` == "GO Component" ~ "CC",    #categoryのGO ComponentをCC
    `category` == "GO Function" ~ "MF",     #categoryのGO FunctionをMF
    `category` == "KEGG Pathways" ~ "KEGG", #categoryのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category`)
  )) %>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description_DBsourseとして結合
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
dat1_f_CC <- dat1 %>% dplyr::filter(category=="GO Component") 
dat1_f_BP <- dat1 %>% dplyr::filter(category=="GO Process") 
dat1_f_MF <- dat1 %>% dplyr::filter(category=="GO Function")



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
dat4_f_KEGG <- dat1 %>% dplyr::filter(category=="KEGG Pathways")
dat5_f <- dat2 %>% #dat2データの
  dplyr::filter(category=="KEGG Pathways") #category列がKEGG Pathways
dat6_f <- dat3 %>% #dat3データの
  dplyr::filter(category=="KEGG Pathways") #category列がKEGG Pathways
#############################################################
#top10データの抽出
dat1_f_BP10 <- dat1_f_BP %>% arrange(desc(log10p)) %>% top_n(n = 10, wt = log10p)
dat1_f_CC10 <- dat1_f_CC %>% arrange(desc(log10p)) %>% top_n(n = 10, wt = log10p)
dat1_f_MF10 <- dat1_f_MF %>% arrange(desc(log10p)) %>% top_n(n = 10, wt = log10p)
dat4_f_KEGG10 <- dat4_f_KEGG %>% arrange(desc(log10p)) %>% top_n(n = 10, wt = log10p)
#top20データの抽出
dat1_ft <- dat1_f %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p) #log10p上位20のデータを保存
dat2_ft <- dat2_f %>% #P_GOのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p) #log10p上位20のデータを保存
dat3_ft <- dat3_f %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p) #log10p上位20のデータを保存
dat4_ft <- dat4_f %>% #P_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p) #log10p上位20のデータを保存
dat5_ft <- dat5_f %>% #C_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p) #log10p上位20のデータを保存
dat6_ft <- dat6_f %>% #PxC_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 20, wt = log10p) #log10p上位20のデータを保存
#############################################################
#top50データの抽出
dat1_ft <- dat1_f %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat2_ft <- dat2_f %>% #P_GOのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat3_ft <- dat3_f %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat4_ft <- dat4_f %>% #P_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat5_ft <- dat5_f %>% #C_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
dat6_ft <- dat6_f %>% #PxC_KEGGのデータから
  arrange(desc(log10p)) %>% #log10pで降順ソート
  top_n(n = 50, wt = log10p) #log10p上位50のデータを保存
#############################################################
#NA(欠損値)の置換
#ifelse(is.na(dat1_ft$log10p), 0, dat1_ft$log10p)
#replace(dat1_ft$log10p, which(is.na(dat1_ft$log10p)), 0)
#dat1_ft$log10p[is.na(dat1_ft$log10p)] <- 0
#############################################################
#top10データをラベル
dat1_f_BP <- mutate(dat1_f_BP, top10 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat1_f_CC <- mutate(dat1_f_CC, top10 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat1_f_MF <- mutate(dat1_f_MF, top10 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat4_f_KEGG <- mutate(dat4_f_KEGG, top10 = if_else(log10p > 0, true = TRUE, false = FALSE))
#top20データをラベル
dat1_f <- mutate(dat1_ft, top20 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat2_f <- mutate(dat2_ft, top20 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat3_f <- mutate(dat3_ft, top20 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat4_f <- mutate(dat4_ft, top20 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat5_f <- mutate(dat5_ft, top20 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat6_f <- mutate(dat6_ft, top20 = if_else(log10p > 0, true = TRUE, false = FALSE))
#############################################################
#top50データをラベル
dat1_f <- mutate(dat1_ft, top50 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat2_f <- mutate(dat2_ft, top50 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat3_f <- mutate(dat3_ft, top50 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat4_f <- mutate(dat4_ft, top50 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat5_f <- mutate(dat5_ft, top50 = if_else(log10p > 0, true = TRUE, false = FALSE))
dat6_f <- mutate(dat6_ft, top50 = if_else(log10p > 0, true = TRUE, false = FALSE))
#############################################################
#dat1,2,3(GO)とdat4,5,6(KEGG)を1つのデータフレーム(GO)に結合
#Top10
GO <- rbind(dat1_f_BP10, dat1_f_CC10, dat1_f_MF10) #data1,2,3結合
KEGG <- dat4_f_KEGG10 #data4,5,6結合
#Top20
GO <- rbind(dat1_f, dat2_f, dat3_f) #data1,2,3結合
KEGG <- rbind(dat4_f, dat5_f, dat6_f) #data4,5,6結合
#############################################################
#重複削除,top10データの抽出
GO_t <- GO %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = 10, wt = log10p) #log10p上位
KEGG_t <- KEGG %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = 10, wt = log10p) #log10p上位
#############################################################
#重複削除,top20データの抽出
GO_t <- GO %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = 20, wt = log10p) #log10p上位
KEGG_t <- KEGG %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = 20, wt = log10p) #log10p上位
#############################################################
#重複削除,top50データの抽出
GO_t <- GO %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = 50, wt = log10p) #log10p上位
KEGG_t <- KEGG %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = 50, wt = log10p) #log10p上位
#############################################################
#topの絞り込み
t(colnames(GO_t))
GO_term <- data.frame(GO_t[,c(3,5,12)])
t(colnames(KEGG_t))
KEGG_term <- data.frame(KEGG_t[,c(3,5,12)])
#Annotation作成
t(colnames(dat1_f))
dat1_f <- dat1_f[,c(-3,-5)]
dat2_f <- dat2_f[,c(-3,-5)]
dat3_f <- dat3_f[,c(-3,-5)]
dat4_f <- dat4_f[,c(-3,-5)]
dat5_f <- dat5_f[,c(-3,-5)]
dat6_f <- dat6_f[,c(-3,-5)]
#データにAnnotationを結合
dat1_t <- left_join(GO_term, dat1_f, by = "term.name") #term.name識別子としてmerge
dat2_t <- left_join(GO_term, dat2_f, by = "term.name") #term.name識別子としてmerge
dat3_t <- left_join(GO_term, dat3_f, by = "term.name") #term.name識別子としてmerge
dat4_t <- left_join(KEGG_term, dat4_f, by = "term.name") #term.name識別子としてmerge
dat5_t <- left_join(KEGG_term, dat5_f, by = "term.name") #term.name識別子としてmerge
dat6_t <- left_join(KEGG_term, dat6_f, by = "term.name") #term.name識別子としてmerge
#############################################################
#グループ挿入
dat1_t$group <- "PCP" #グループ挿入
dat2_t$group <- "CLZ" #グループ挿入
dat3_t$group <- "Interaction" #グループ挿入
dat4_t$group <- "PCP" #グループ挿入
dat5_t$group <- "CLZ" #グループ挿入
dat6_t$group <- "Interaction" #グループ挿入
#############################################################
#dat1,2,3(GO)とdat4,5,6(KEGG)を1つのデータフレーム(GO)に再結合
GO_t <- rbind(dat1_t, dat2_t, dat3_t) #dat1,2,3結合
KEGG_t <- rbind(dat4_t, dat5_t, dat6_t) #dat4,5,6結合
#############################################################
#NA(欠損値)の置換
ifelse(is.na(GO_t$log10p), 0, GO_t$log10p)
replace(GO_t$log10p, which(is.na(GO_t$log10p)), 0)
GO_t$log10p[is.na(GO_t$log10p)] <- 0
ifelse(is.na(KEGG_t$log10p), 0, KEGG_t$log10p)
replace(KEGG_t$log10p, which(is.na(KEGG_t$log10p)), 0)
KEGG_t$log10p[is.na(KEGG_t$log10p)] <- 0
#############################################################
#NA(欠損値)の置換
ifelse(is.na(GO_t$`Overlap`), 0, GO_t$`Overlap`)
replace(GO_t$`Overlap`, which(is.na(GO_t$`Overlap`)), 0)
GO_t$`Overlap`[is.na(GO_t$`Overlap`)] <- 0
ifelse(is.na(KEGG_t$`Overlap`), 0, KEGG_t$`Overlap`)
replace(KEGG_t$`Overlap`, which(is.na(KEGG_t$`Overlap`)), 0)
KEGG_t$`Overlap`[is.na(KEGG_t$`Overlap`)] <- 0
#############################################################
#descriptionの順序の入れ替えのために
#データをファクタに変換する
#library(dplyr)
dat1_ft <- dat1_ft %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
dat2_ft <- dat2_ft %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
dat3_ft <- dat3_ft %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
dat4_ft <- dat4_ft %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
dat5_ft <- dat5_ft %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
dat6_ft <- dat6_ft %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
GO_t <- GO_t %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
KEGG_t <- KEGG_t %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
#############################################################
#列の作成と結合
GO_t <- GO_t %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category` == "GO Process" ~ "BP",      #categoryのGO ProcessをBP
    `category` == "GO Component" ~ "CC",    #categoryのGO ComponentをCC
    `category` == "GO Function" ~ "MF",     #categoryのGO FunctionをMF
    `category` == "KEGG Pathways" ~ "KEGG", #categoryのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category`)
  )) %>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description_DBsourseとして結合
KEGG_t <- KEGG_t %>% as.data.frame() %>%
  mutate(DBsourse = case_when(                  #DBsourse列作成
    `category` == "GO Process" ~ "BP",      #categoryのGO ProcessをBP
    `category` == "GO Component" ~ "CC",    #categoryのGO ComponentをCC
    `category` == "GO Function" ~ "MF",     #categoryのGO FunctionをMF
    `category` == "KEGG Pathways" ~ "KEGG", #categoryのKEGG PathwaysをKEGG
    TRUE ~ as.character(`category`)
  )) %>%
  mutate(term = paste(                                #term列作成
    !!!rlang::syms(c("description", "DBsourse")), #シンボル化：map(as.symbol)と同じ
    sep="_"))                                         #description_DBsourseとして結合
#############################################################
#ggplot2,GO
plot1 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat1_f_BP10, aes(x = log10p, 
                              y = reorder(x = term, X = log10p), 
                              #position = "jitter", #geom_jitter
                              color = group,
                              alpha = 0.9,
                              size = `Overlap`))
plot2 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat1_f_CC10, aes(x = log10p, 
                                     y = reorder(x = term, X = log10p), 
                                     #position = "jitter", #geom_jitter
                                     color = group,
                                     alpha = 0.9,
                                     size = `Overlap`))
plot3 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat1_f_MF10, aes(x = log10p, 
                                     y = reorder(x = term, X = log10p), 
                                     #position = "jitter", #geom_jitter
                                     color = group,
                                     alpha = 0.9,
                                     size = `Overlap`))






#ggplot2,GO
plot1 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = GO_t, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group,
                                 alpha = 0.9,
                                 size = `Overlap`))
plot2 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat1_ft, aes(x = log10p, 
                                y = reorder(x = term, X = log10p), 
                                #position = "jitter", #geom_jitter
                                color = group,
                                alpha = 0.9,
                                size = `Overlap`))
plot3 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat2_ft, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group,
                                 alpha = 0.9,
                                 size = `Overlap`))
plot4 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat3_ft, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group,
                                 alpha = 0.9, #プロット透明度設定
                                 size = `Overlap`))
#############################################################
#ggplot2,KEGG
plot5 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat4_f_KEGG10, aes(x = log10p,
                                y = reorder(x = term, X = log10p),
                                #position = "jitter", #geom_jitter
                                color = group,
                                alpha = 0.9,
                                size = `Overlap`))


#ggplot2,KEGG
plot5 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = KEGG_t, aes(x = log10p,
                                y = reorder(x = term, X = log10p),
                                #position = "jitter", #geom_jitter
                                color = group,
                                alpha = 0.9,
                                size = `Overlap`))
plot6 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat4_ft, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group,
                                 alpha = 0.9,
                                 size = `Overlap`))
plot7 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat5_ft, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group,
                                 alpha = 0.9,
                                 size = `Overlap`))
plot8 <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = dat6_ft, aes(x = log10p, 
                                 y = reorder(x = term, X = log10p), 
                                 #position = "jitter", #geom_jitter
                                 color = group,
                                 alpha = 0.9,
                                 size = `Overlap`))
############################################################################################
#svg出力
svg(file="BP_top10.svg") #ファイル名指定
print(plot1) #プロット作成
dev.off() #svg出力
svg(file="CC_top10.svg") #ファイル名指定
print(plot2) #プロット作成
dev.off() #svg出力
svg(file="MF_top10.svg") #ファイル名指定
print(plot3) #プロット作成
dev.off() #svg出力
svg(file="KEGG_top10.svg") #ファイル名指定
print(plot5) #プロット作成
dev.off() #svg出力





#svg出力
svg(file="stringAppGO_top20.svg") #ファイル名指定
print(plot1) #プロット作成
dev.off() #svg出力
svg(file="stringAppGO_P_top20.svg") #ファイル名指定
print(plot2) #プロット作成
dev.off() #svg出力
svg(file="stringAppGO_C_top20.svg") #ファイル名指定
print(plot3) #プロット作成
dev.off() #svg出力
svg(file="stringAppGO_PxC_top20.svg") #ファイル名指定
print(plot4) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_top20.svg") #ファイル名指定
print(plot5) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_P_top20.svg") #ファイル名指定
print(plot6) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_C_top20.svg") #ファイル名指定
print(plot7) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_PxC_top20.svg") #ファイル名指定
print(plot8) #プロット作成
dev.off() #svg出力
#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("BP_top10.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot1) #プロット作成
dev.off() #PDF出力
pdf("CC_top10.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot2) #プロット作成
dev.off() #PDF出力
pdf("MF_top10.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot3) #プロット作成
dev.off() #PDF出力
pdf("KEGG_top10.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot5) #プロット作成
dev.off() #PDF出力








#レシピ14.1 PDFベクタファイルへの出力
pdf("stringAppGO_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot1) #プロット作成
dev.off() #PDF出力
pdf("stringAppGO_P_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot2) #プロット作成
dev.off() #PDF出力
pdf("stringAppGO_C_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot3) #プロット作成
dev.off() #PDF出力
pdf("stringAppGO_PxC_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot4) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot5) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_P_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot6) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_C_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot7) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_PxC_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot8) #プロット作成
dev.off() #PDF出力
#############################################################
#ggplot2
#volcano plot
plot1 <- ggplot(dat1_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top20,
                            colour = top20, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat1_f$`Overlap` * 0.05,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 150)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 35))    #y軸の範囲設定
vol1 #グラフをみて範囲を調整
plot2 <- ggplot(dat2_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top20,
                            colour = top20, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat2_f$`Overlap` * 0.02,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 200)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 60))    #y軸の範囲設定
vol2 #グラフをみて範囲を調整
plot3 <- ggplot(dat3_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top20,
                            colour = top20, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat3_f$`Overlap` * 0.01,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 175)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 80))    #y軸の範囲設定
vol3 #グラフをみて範囲を調整
plot4 <- ggplot(dat4_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top20,
                            colour = top20, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat4_f$`Overlap` * 0.2,  #プロットサイズ設定
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
vol4 #グラフをみて範囲を調整
plot5 <- ggplot(dat5_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top20,
                            colour = top20, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat5_f$`Overlap` * 0.2,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 30)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 20))    #y軸の範囲設定
vol5 #グラフをみて範囲を調整
plot6 <- ggplot(dat6_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top20,
                            colour = top20, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat6_f$`Overlap` * 0.07,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 30)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 20))    #y軸の範囲設定
vol6 #グラフをみて範囲を調整
#############################################################
#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("stringApp_vol1_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol1) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol2_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol2) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol3_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol3) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol4_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol4) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol5_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol5) #プロット作成
dev.off() #PDF出力
pdf("stringApp_vol6_top20.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(vol6) #プロット作成
dev.off() #PDF出力
#############################################################
#xslx出力
smp <- list("GO"=GO_t, "GO_P"=dat1_ft, "GO_C"=dat2_ft, "GO_PxC"=dat3_ft, "KEGG"=KEGG_t, "KEGG_P"=dat4_ft, "KEGG_C"=dat5_ft, "KEGG_PxC"=dat6_ft) #GO,KEGGのtop20リスト作成
write.xlsx(smp, "GO_KEGG_top20.xlsx") #GOシート,KEGGシート出力
#############################################################
#エクセルでRの特殊記号"|"を","に置換し、GO_KEGGr.xlsxとして保存！！！









#############################################################
#svg出力
svg(file="stringAppGO.svg") #ファイル名指定
print(plot1) #プロット作成
dev.off() #svg出力
svg(file="stringAppGO_P.svg") #ファイル名指定
print(plot2) #プロット作成
dev.off() #svg出力
svg(file="stringAppGO_C.svg") #ファイル名指定
print(plot3) #プロット作成
dev.off() #svg出力
svg(file="stringAppGO_PxC.svg") #ファイル名指定
print(plot4) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG.svg") #ファイル名指定
print(plot5) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_P.svg") #ファイル名指定
print(plot6) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_C.svg") #ファイル名指定
print(plot7) #プロット作成
dev.off() #svg出力
svg(file="stringAppKEGG_PxC.svg") #ファイル名指定
print(plot8) #プロット作成
dev.off() #svg出力
#P359################################################################
#レシピ14.1 PDFベクタファイルへの出力
pdf("stringAppGO.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot1) #プロット作成
dev.off() #PDF出力
pdf("stringAppGO_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot2) #プロット作成
dev.off() #PDF出力
pdf("stringAppGO_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot3) #プロット作成
dev.off() #PDF出力
pdf("stringAppGO_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot4) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot5) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot6) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot7) #プロット作成
dev.off() #PDF出力
pdf("stringAppKEGG_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plot8) #プロット作成
dev.off() #PDF出力
#############################################################
#ggplot2
#volcano plot
plot1 <- ggplot(dat1_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat1_f$`Overlap` * 0.05,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 150)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 35))    #y軸の範囲設定
vol1 #グラフをみて範囲を調整
plot2 <- ggplot(dat2_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat2_f$`Overlap` * 0.02,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 150)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 35))    #y軸の範囲設定
vol2 #グラフをみて範囲を調整
plot3 <- ggplot(dat3_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat3_f$`Overlap` * 0.01,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 350)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 70))    #y軸の範囲設定
vol3 #グラフをみて範囲を調整
plot4 <- ggplot(dat4_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat4_f$`Overlap` * 0.2,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 20)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 10))    #y軸の範囲設定
vol4 #グラフをみて範囲を調整
plot5 <- ggplot(dat5_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat5_f$`Overlap` * 0.2,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 30)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 15))    #y軸の範囲設定
vol5 #グラフをみて範囲を調整
plot6 <- ggplot(dat6_f, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top50,
                            colour = top50, 
                            fill = `Overlap` #カラー設定
)) +
  geom_point(size = dat6_f$`Overlap` * 0.07,  #プロットサイズ設定
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
  scale_x_continuous(limits = c(0, 60)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 40))    #y軸の範囲設定
vol6 #グラフをみて範囲を調整
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
smp <- list("GO"=GO_t, "GO_P"=dat1_ft, "GO_C"=dat2_ft, "GO_PxC"=dat3_ft, "KEGG"=KEGG_t, "KEGG_P"=dat4_ft, "KEGG_C"=dat5_ft, "KEGG_PxC"=dat6_ft) #GO,KEGGのtop50リスト作成
write.xlsx(smp, "GO_KEGG.xlsx") #GOシート,KEGGシート出力
#############################################################
#エクセルでRの特殊記号"|"を","に置換し、GO_KEGGr.xlsxとして保存！！！

