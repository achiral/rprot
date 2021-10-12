#stringApp_01_GO_KEGG_volcanoplot.R
#stringApp1.5.1_Cytoscape3.8.0→stringApp1.6.0_Cytoscape3.8.1_Java11.0.6
#stringApp出力データから
#1_P,C,IntにおけるGO_KEGG_top50ラベル
#2_P,C,Int全てについてGO_KEGG_top50抽出
#3_GO_KEGG termのplot描写
#4_P,C,Intにおけるvolcano plot描写(top50赤色)
#5_GO_KEGG_top50に含まれるDEP list作成 →別script
################################################################################
rm(list = ls(all.names = TRUE))
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR/stringApp") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC/stringApp") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
options(digits=2) #change digit2桁表示指定
################################################################################
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
################################################################################
#CSV入力
dat1 <- read.csv("String_twANOVA_Pq005_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
dat2 <- read.csv("String_twANOVA_Cq005_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
dat3 <- read.csv("String_twANOVA_PxCq005_filt.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
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
################################################################################
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
################################################################################
#GO抽出
dat1_GO <- dat1 %>% #dat1データの
  dplyr::filter(category=="GO Component" | #category列がGO Componentまたは
                  category=="GO Process" | #category列がGO Processまたは
                  category=="GO Function") #category列がGO Function
dat1_CC <- dat1 %>% dplyr::filter(category=="GO Component") 
dat1_BP <- dat1 %>% dplyr::filter(category=="GO Process") 
dat1_MF <- dat1 %>% dplyr::filter(category=="GO Function")

dat2_GO <- dat2 %>% #dat2データの
  dplyr::filter(category=="GO Component" | #category列がGO Componentまたは
                  category=="GO Process" | #category列がGO Processまたは
                  category=="GO Function") #category列がGO Function
dat2_CC <- dat2 %>% dplyr::filter(category=="GO Component") 
dat2_BP <- dat2 %>% dplyr::filter(category=="GO Process") 
dat2_MF <- dat2 %>% dplyr::filter(category=="GO Function")

dat3_GO <- dat3 %>% #dat3データの
  dplyr::filter(category=="GO Component" | #category列がGO Componentまたは
                  category=="GO Process" | #category列がGO Processまたは
                  category=="GO Function") #category列がGO Function
dat3_CC <- dat3 %>% dplyr::filter(category=="GO Component") 
dat3_BP <- dat3 %>% dplyr::filter(category=="GO Process") 
dat3_MF <- dat3 %>% dplyr::filter(category=="GO Function")

#KEGG抽出
dat4_KEGG <- dat1 %>% dplyr::filter(category=="KEGG Pathways")#dat1データのcategory列がKEGG Pathways
dat5_KEGG <- dat2 %>% dplyr::filter(category=="KEGG Pathways")#dat2データのcategory列がKEGG Pathways
dat6_KEGG <- dat3 %>% dplyr::filter(category=="KEGG Pathways")#dat3データのcategory列がKEGG Pathways
################################################################################
#topデータの抽出
NUM <- 10
dat1_GO_t <- dat1_GO %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p) #log10pで降順ソートし上位のデータを保存
dat1_BP_t <- dat1_BP %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat1_CC_t <- dat1_CC %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat1_MF_t <- dat1_MF %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat2_GO_t <- dat2_GO %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p) 
dat2_BP_t <- dat2_BP %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat2_CC_t <- dat2_CC %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat2_MF_t <- dat2_MF %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat3_GO_t <- dat3_GO %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat3_BP_t <- dat3_BP %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat3_CC_t <- dat3_CC %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat3_MF_t <- dat3_MF %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat4_KEGG_t <- dat4_KEGG %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat5_KEGG_t <- dat5_KEGG %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)
dat6_KEGG_t <- dat6_KEGG %>% arrange(desc(log10p)) %>% top_n(n = NUM, wt = log10p)


################################################################################
#topデータをラベル
dat1_GO_t <- dat1_GO_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat1_BP_t <- dat1_BP_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat1_CC_t <- dat1_CC_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat1_MF_t <- dat1_MF_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat2_GO_t <- dat2_GO_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat2_BP_t <- dat2_BP_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat2_CC_t <- dat2_CC_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat2_MF_t <- dat2_MF_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat3_GO_t <- dat3_GO_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat3_BP_t <- dat3_BP_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat3_CC_t <- dat3_CC_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat3_MF_t <- dat3_MF_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat4_KEGG_t <- dat4_KEGG_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat5_KEGG_t <- dat5_KEGG_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
dat6_KEGG_t <- dat6_KEGG_t %>% mutate(top = if_else(log10p > 0, true = TRUE, false = FALSE))
################################################################################
#dat1,2,3(GO)とdat4,5,6(KEGG)を1つのデータフレーム(GO)に結合
GO_t <- rbind(dat1_GO_t, dat2_GO_t, dat3_GO_t) #data1,2,3結合
KEGG_t <- rbind(dat4_KEGG_t, dat5_KEGG_t, dat6_KEGG_t) #data4,5,6結合
################################################################################
#重複削除,topデータの抽出
GO_t <- GO_t %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = NUM, wt = log10p) #log10p上位
KEGG_t <- KEGG_t %>%
  arrange(desc(log10p)) %>% #log10pで降順ソート
  distinct(description,.keep_all=TRUE) %>% #重複削除
  top_n(n = NUM, wt = log10p) #log10p上位
################################################################################
#topの絞り込み
t(colnames(GO_t))
GO_term <- data.frame(GO_t[,c(3,5,12)])
t(colnames(KEGG_t))
KEGG_term <- data.frame(KEGG_t[,c(3,5,12)])
#Annotation作成
t(colnames(dat1_GO))
#データ処理
dat1_GO_rm <- dat1_GO[,c(-3,-5)]
dat2_GO_rm <- dat2_GO[,c(-3,-5)]
dat3_GO_rm <- dat3_GO[,c(-3,-5)]
dat4_KEGG_rm <- dat4_KEGG[,c(-3,-5)]
dat5_KEGG_rm <- dat5_KEGG[,c(-3,-5)]
dat6_KEGG_rm <- dat6_KEGG[,c(-3,-5)]
#データにAnnotationを結合
#欠損値
dat1_t_GO <- left_join(GO_term, dat1_GO_rm, by = "term.name") #term.name識別子としてmerge
dat2_t_GO <- left_join(GO_term, dat2_GO_rm, by = "term.name") #term.name識別子としてmerge
dat3_t_GO <- left_join(GO_term, dat3_GO_rm, by = "term.name") #term.name識別子としてmerge
dat4_t_KEGG <- left_join(KEGG_term, dat4_KEGG_rm, by = "term.name") #term.name識別子としてmerge
dat5_t_KEGG <- left_join(KEGG_term, dat5_KEGG_rm, by = "term.name") #term.name識別子としてmerge
dat6_t_KEGG <- left_join(KEGG_term, dat6_KEGG_rm, by = "term.name") #term.name識別子としてmerge
################################################################################
#グループ挿入
dat1_t_GO$group <- "PCP"
dat2_t_GO$group <- "CLZ"
dat3_t_GO$group <- "Interaction"
dat4_t_KEGG$group <- "PCP"
dat5_t_KEGG$group <- "CLZ"
dat6_t_KEGG$group <- "Interaction"
################################################################################
#dat1,2,3(GO)とdat4,5,6(KEGG)を1つのデータフレーム(GO)に再結合
GO_t <- rbind(dat1_t_GO, dat2_t_GO, dat3_t_GO) #dat1,2,3結合
KEGG_t <- rbind(dat4_t_KEGG, dat5_t_KEGG, dat6_t_KEGG) #dat4,5,6結合
################################################################################
#NA(欠損値)の置換
#ifelse(is.na(GO_t$log10p), 0, GO_t$log10p)
#replace(GO_t$log10p, which(is.na(GO_t$log10p)), 0)
GO_t$log10p[is.na(GO_t$log10p)] <- 0
#ifelse(is.na(KEGG_t$log10p), 0, KEGG_t$log10p)
#replace(KEGG_t$log10p, which(is.na(KEGG_t$log10p)), 0)
KEGG_t$log10p[is.na(KEGG_t$log10p)] <- 0
################################################################################
#NA(欠損値)の置換
#ifelse(is.na(GO_t$`Overlap`), 0, GO_t$`Overlap`)
#replace(GO_t$`Overlap`, which(is.na(GO_t$`Overlap`)), 0)
GO_t$`Overlap`[is.na(GO_t$`Overlap`)] <- 0
#ifelse(is.na(KEGG_t$`Overlap`), 0, KEGG_t$`Overlap`)
#replace(KEGG_t$`Overlap`, which(is.na(KEGG_t$`Overlap`)), 0)
KEGG_t$`Overlap`[is.na(KEGG_t$`Overlap`)] <- 0
################################################################################
#descriptionの順序の入れ替えのために
#データをファクタに変換する
#library(dplyr)
dat1_GO_t <- dat1_GO_t %>% mutate(description = as.factor(description)) #descriptionをファクタに変換
dat2_GO_t <- dat2_GO_t %>% mutate(description = as.factor(description))
dat3_GO_t <- dat3_GO_t %>% mutate(description = as.factor(description))
dat4_KEGG_t <- dat4_KEGG_t %>% mutate(description = as.factor(description)) 
dat5_KEGG_t <- dat5_KEGG_t %>% mutate(description = as.factor(description))
dat6_KEGG_t <- dat6_KEGG_t %>% mutate(description = as.factor(description))
GO_t <- GO_t %>% mutate(description = as.factor(description))
KEGG_t <- KEGG_t %>% mutate(description = as.factor(description))
################################################################################
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
################################################################################
#ggplot2,GO
plotGO <- ggplot()+
  theme_gray()+ #フォント指定base_family = "Arial", "HiraKakuPro-W3"など
  geom_point(data = GO_t, aes(x = log10p, 
                              y = reorder(x = term, X = log10p), 
                              #position = "jitter", #geom_jitter
                              color = group,
                              alpha = 0.9,
                              size = `Overlap`))
plotGO1 <- ggplot() + theme_gray() + geom_point(data = dat1_GO_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotGO2 <- ggplot() + theme_gray() + geom_point(data = dat2_GO_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotGO3 <- ggplot() + theme_gray() + geom_point(data = dat3_GO_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotBP1 <- ggplot() + theme_gray() + geom_point(data = dat1_BP_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotBP2 <- ggplot() + theme_gray() + geom_point(data = dat2_BP_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotBP3 <- ggplot() + theme_gray() + geom_point(data = dat3_BP_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotCC1 <- ggplot() + theme_gray() + geom_point(data = dat1_CC_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotCC2 <- ggplot() + theme_gray() + geom_point(data = dat2_CC_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotCC3 <- ggplot() + theme_gray() + geom_point(data = dat3_CC_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotMF1 <- ggplot() + theme_gray() + geom_point(data = dat1_MF_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotMF2 <- ggplot() + theme_gray() + geom_point(data = dat2_MF_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotMF3 <- ggplot() + theme_gray() + geom_point(data = dat3_MF_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
################################################################################
#ggplot2,KEGG
plotKEGG <- ggplot() + theme_gray() + geom_point(data = KEGG_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotKEGG4 <- ggplot() + theme_gray() + geom_point(data = dat4_KEGG_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotKEGG5 <- ggplot() + theme_gray() + geom_point(data = dat5_KEGG_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
plotKEGG6 <- ggplot() + theme_gray() + geom_point(data = dat6_KEGG_t, aes(x = log10p, y = reorder(x = term, X = log10p), color = group, alpha = 0.9, size = `Overlap`))
################################################################################
#svg出力
#GO
#svg(file="stringApp_GO_top_all.svg") #ファイル名指定
svg(file="stringApp_GO_top.svg") #ファイル名指定
print(plotGO) #プロット作成
dev.off() #svg出力
#svg(file="stringApp_GO_top_P.svg") #ファイル名指定
#print(plotGO1) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_GO_top_C.svg") #ファイル名指定
#print(plotGO2) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_GO_top_PxC.svg") #ファイル名指定
#print(plotGO3) #プロット作成
#dev.off() #svg出力
#BP
#svg(file="stringApp_BP_top_P.svg") #ファイル名指定
#print(plotBP1) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_BP_top_C.svg") #ファイル名指定
#print(plotBP2) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_BP_top_PxC.svg") #ファイル名指定
#print(plotBP3) #プロット作成
#dev.off() #svg出力
#CC
#svg(file="stringApp_CC_top_P.svg") #ファイル名指定
#print(plotCC1) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_CC_top_C.svg") #ファイル名指定
#print(plotCC2) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_CC_top_PxC.svg") #ファイル名指定
#print(plotCC3) #プロット作成
#dev.off() #svg出力
#MF
#svg(file="stringApp_MF_top_P.svg") #ファイル名指定
#print(plotMF1) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_MF_top_C.svg") #ファイル名指定
#print(plotMF2) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_MF_top_PxC.svg") #ファイル名指定
#print(plotMF3) #プロット作成
#dev.off() #svg出力
#KEGG
#svg(file="stringApp_KEGG_top_all.svg") #ファイル名指定
svg(file="stringApp_KEGG_top.svg") #ファイル名指定
print(plotKEGG) #プロット作成
dev.off() #svg出力
#svg(file="stringApp_KEGG_top_P.svg") #ファイル名指定
#print(plotKEGG4) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_KEGG_top_C.svg") #ファイル名指定
#print(plotKEGG5) #プロット作成
#dev.off() #svg出力
#svg(file="stringApp_KEGG_top_PxC.svg") #ファイル名指定
#print(plotKEGG6) #プロット作成
#dev.off() #svg出力
#P359###########################################################################
#レシピ14.1 PDFベクタファイルへの出力
#GO
pdf(file="stringApp_GO_top_all.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
#pdf(file="stringApp_GO_top.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotGO) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_GO_top_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotGO1) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_GO_top_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotGO2) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_GO_top_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotGO3) #プロット作成
#dev.off() #pdf出力
#BP
#pdf(file="stringApp_BP_top_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotBP1) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_BP_top_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotBP2) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_BP_top_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotBP3) #プロット作成
#dev.off() #pdf出力
#CC
#pdf(file="stringApp_CC_top_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotCC1) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_CC_top_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotCC2) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_CC_top_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotCC3) #プロット作成
#dev.off() #pdf出力
#MF
#pdf(file="stringApp_MF_top_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotMF1) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_MF_top_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotMF2) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_MF_top_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotMF3) #プロット作成
dev.off() #pdf出力
#KEGG
pdf(file="stringApp_KEGG_top_all.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
#pdf(file="stringApp_KEGG_top.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotKEGG) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_KEGG_top_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotKEGG4) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_KEGG_top_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotKEGG5) #プロット作成
#dev.off() #pdf出力
#pdf(file="stringApp_KEGG_top_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(plotKEGG6) #プロット作成
dev.off() #pdf出力
################################################################################
#ggplot2
#volcano plot
volGO1 <- ggplot(dat1_GO_t, aes(x = `Overlap`, y = log10p, #xy値設定
                            shape = top, colour = top, fill = `Overlap` #カラー設定
                            )) + geom_point(size = dat1_GO_t$`Overlap` * 0.05, #プロットサイズ設定
                                            alpha = 0.8) + #プロット透明度設定
  theme_bw()+ scale_shape_manual(values = c(21,21)) + #プロット形状設定
  scale_colour_manual(values = c("black","red")) + scale_fill_gradient2(low = "dark grey", mid = "black", high  = "dark red", midpoint = 1) + 
  scale_x_continuous(limits = c(0, 150)) + #x軸の範囲設定
  scale_y_continuous(limits = c(0, 35)) #y軸の範囲設定
volGO1 #グラフをみて範囲を調整
volGO2 <- ggplot(dat2_GO_t, aes(x = `Overlap`, y = log10p, shape = top, colour = top, fill = `Overlap`)) + geom_point(size = dat1_GO_t$`Overlap` * 0.05, alpha = 0.8) + theme_bw()+ scale_shape_manual(values = c(21,21)) + scale_colour_manual(values = c("black","red")) + scale_fill_gradient2(low = "dark grey", mid = "black", high  = "dark red", midpoint = 1) + 
  scale_x_continuous(limits = c(0, 150)) + 
  scale_y_continuous(limits = c(0, 35))
volGO2
volGO3 <- ggplot(dat3_GO_t, aes(x = `Overlap`, y = log10p, shape = top, colour = top, fill = `Overlap`)) + geom_point(size = dat1_GO_t$`Overlap` * 0.05, alpha = 0.8) + theme_bw()+ scale_shape_manual(values = c(21,21)) + scale_colour_manual(values = c("black","red")) + scale_fill_gradient2(low = "dark grey", mid = "black", high  = "dark red", midpoint = 1) + 
  scale_x_continuous(limits = c(0, 150)) + 
  scale_y_continuous(limits = c(0, 35))
volGO3
volKEGG4 <- ggplot(dat4_KEGG_t, aes(x = `Overlap`, y = log10p, shape = top, colour = top, fill = `Overlap`)) + geom_point(size = dat1_GO_t$`Overlap` * 0.05, alpha = 0.8) + theme_bw()+ scale_shape_manual(values = c(21,21)) + scale_colour_manual(values = c("black","red")) + scale_fill_gradient2(low = "dark grey", mid = "black", high  = "dark red", midpoint = 1) + 
  scale_x_continuous(limits = c(0, 150)) + 
  scale_y_continuous(limits = c(0, 35))
volKEGG4
volKEGG5 <- ggplot(dat5_KEGG_t, aes(x = `Overlap`, y = log10p, shape = top, colour = top, fill = `Overlap`)) + geom_point(size = dat1_GO_t$`Overlap` * 0.05, alpha = 0.8) + theme_bw()+ scale_shape_manual(values = c(21,21)) + scale_colour_manual(values = c("black","red")) + scale_fill_gradient2(low = "dark grey", mid = "black", high  = "dark red", midpoint = 1) + 
  scale_x_continuous(limits = c(0, 150)) + 
  scale_y_continuous(limits = c(0, 35))
volKEGG5
volKEGG6 <- ggplot(dat6_KEGG_t, aes(x = `Overlap`, y = log10p, shape = top, colour = top, fill = `Overlap`)) + geom_point(size = dat1_GO_t$`Overlap` * 0.05, alpha = 0.8) + theme_bw()+ scale_shape_manual(values = c(21,21)) + scale_colour_manual(values = c("black","red")) + scale_fill_gradient2(low = "dark grey", mid = "black", high  = "dark red", midpoint = 1) + 
  scale_x_continuous(limits = c(0, 150)) + 
  scale_y_continuous(limits = c(0, 35))
volKEGG6
################################################################################
#P359###########################################################################
#レシピ14.1 PDFベクタファイルへの出力
#GO
pdf("stringApp_volGO_t_all.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
#pdf("stringApp_volGO_t_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(volGO1) #プロット作成
#dev.off() #PDF出力
#pdf("stringApp_volGO_t_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(volGO2) #プロット作成
#dev.off() #PDF出力
#pdf("stringApp_volGO_t_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(volGO3) #プロット作成
dev.off() #PDF出力
#KEGG
pdf("stringApp_volKEGG_t_all.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
#pdf("stringApp_volKEGG_t_P.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(volKEGG4) #プロット作成
#dev.off() #PDF出力
#pdf("stringApp_volKEGG_t_C.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(volKEGG5) #プロット作成
#dev.off() #PDF出力
#pdf("stringApp_volKEGG_t_PxC.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE) #幅,高さcm(inch/2.54)
print(volKEGG6) #プロット作成
dev.off() #PDF出力
################################################################################
#xslx出力
smp <- list("GO"=GO_t, "GO_P"=dat1_GO_t, "GO_C"=dat2_GO_t, "GO_PxC"=dat3_GO_t, "KEGG"=KEGG_t, "KEGG_P"=dat4_KEGG_t, "KEGG_C"=dat5_KEGG_t, "KEGG_PxC"=dat6_KEGG_t) #GO,KEGGのtop_tリスト作成
write.xlsx(smp, "GO_KEGG_top_t.xlsx") #GOシート,KEGGシート出力
################################################################################
#エクセルでRの特殊記号"|"を","に置換し、GO_KEGGr.xlsxとして保存！！！

