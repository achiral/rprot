#ClueGO_04_GO_KEGG_top20_proteinlist.R
#GO,KEGGのTop20に含まれる遺伝子リストの作成 -> GO_KEGG_list.csv
#Rでデータフレームからデータを抽出(検索)
#https://kazutan.github.io/kazutanR/hands_on_170730/filter.html#na処理
#https://qiita.com/gigatune/items/f3aa0afef7f50ab791cd
#http://www.restorative-pt.tokyo/archives/18893392.html
#行列抽出
#https://qiita.com/hitsujisuke/items/d71ee40daa0786ae1680
#http://www.restorative-pt.tokyo/archives/r_select_slice_filter.html
#http://tips-r.blogspot.com/2018/02/r.html
#http://www.restorative-pt.tokyo/archives/r_select_slice_filter.html
#分割
#https://www.it-swarm.dev/ja/r/データフレーム文字列列を複数の列に分割/970509188/
#重複削除
#https://a-habakiri.hateblo.jp/entry/2016/11/29/215013
#積集合
#https://ja.coder.work/so/r/3195287
#データフレーム ⇔ リスト の変換
#https://qiita.com/U25CE/items/bf30ea2fcc79ba399dd0
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/SWATH")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
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
data1 <- read_excel("GO_KEGG_top2.xlsx", 1) #シート1の読み込み
data2 <- read_excel("GO_KEGG_top2.xlsx", sheet = 2) #シート2の読み込み
options(digits=2) #change digit2桁表示指定
#str(ex1) #データ構造確認
#############################################################
dat1 <- filter(data1, !is.na(`Associated Genes Found`), group2 == "PCP") #NA除くPCPのtop20リスト
dat2 <- filter(data1, !is.na(`Associated Genes Found`), group2 == "CLZ") #NA除くCLZのtop20リスト
dat3 <- filter(data1, !is.na(`Associated Genes Found`), group2 == "Interaction") #NA除くInteractionのtop20リスト
dat4 <- filter(data2, !is.na(`Associated Genes Found`), group2 == "PCP") #NA除くPCPのtop20リスト
dat5 <- filter(data2, !is.na(`Associated Genes Found`), group2 == "CLZ") #NA除くCLZのtop20リスト
dat6 <- filter(data2, !is.na(`Associated Genes Found`), group2 == "Interaction") #NA除くInteractionのtop20リスト

#行の抽出(CNS, Endcrine_Metabolism, Cardiovascular)
dat1[,2] #2列目表示
dat2[,2] #2列目表示
index1 <- c(1,3,5:9,12,17) #CNSrelated
index2 <- c(1,3,4,8,18) #CNSrelated
index3 <- c(14,16) #CVrelated
index4 <- c(2,5,9,12,15,20) #EMrelated
cns1 <- dat1 %>% 
  slice(index1) %>% #CNSrelated
  select(41) #genelist列抽出
cns2 <- dat2 %>% 
  slice(index1) %>% #CNSrelated
  select(41) #genelist列抽出
cns3 <- dat3 %>% 
  slice(index1) %>% #CNSrelated
  select(41) #genelist列抽出
cns4 <- dat4 %>% 
  slice(index2) %>% #CNSrelated
  select(41) #genelist列抽出
cns5 <- dat5 %>% 
  slice(index2) %>% #CNSrelated
  select(41) #genelist列抽出
cns6 <- dat6 %>% 
  slice(index2) %>%#CNSrelated
  select(41) #genelist列抽出
cv4 <- dat4 %>% 
  slice(index3) %>% #CVrelated
  select(41) #genelist列抽出
cv5 <- dat5 %>% 
  slice(index3) %>% #CVrelated
  select(41) #genelist列抽出
cv6 <- dat6 %>% 
  slice(index3) %>% #CVrelated
  select(41) #genelist列抽出
em4 <- dat4 %>% 
  slice(index4) %>% #CVrelated
  select(41) #genelist列抽出
em5 <- dat5 %>% 
  slice(index4) %>% #CVrelated
  select(41) #genelist列抽出
em6 <- dat6 %>% 
  slice(index4) %>% #CVrelated
  select(41) #genelist列抽出
colnames(cns1) <- "AssociatedGenesFound" #列名変更
colnames(cns2) <- "AssociatedGenesFound" #列名変更
colnames(cns3) <- "AssociatedGenesFound" #列名変更
colnames(cns4) <- "AssociatedGenesFound" #列名変更
colnames(cns5) <- "AssociatedGenesFound" #列名変更
colnames(cns6) <- "AssociatedGenesFound" #列名変更
colnames(cv4) <- "AssociatedGenesFound" #列名変更
colnames(cv5) <- "AssociatedGenesFound" #列名変更
colnames(cv6) <- "AssociatedGenesFound" #列名変更
colnames(em4) <- "AssociatedGenesFound" #列名変更
colnames(em5) <- "AssociatedGenesFound" #列名変更
colnames(em6) <- "AssociatedGenesFound" #列名変更

#列抽出
#t(names(data1_P)) #行列の名前確認
dat1 <- dat1 %>% select(41)
dat2 <- dat2 %>% select(41)
dat3 <- dat3 %>% select(41)
dat4 <- dat4 %>% select(41)
dat5 <- dat5 %>% select(41)
dat6 <- dat6 %>% select(41)
colnames(dat1) <- "AssociatedGenesFound" #列名変更
colnames(dat2) <- "AssociatedGenesFound" #列名変更
colnames(dat3) <- "AssociatedGenesFound" #列名変更
colnames(dat4) <- "AssociatedGenesFound" #列名変更
colnames(dat5) <- "AssociatedGenesFound" #列名変更
colnames(dat6) <- "AssociatedGenesFound" #列名変更

#縦に結合
GOcns <- rbind(cns1, cns2)
GOcns <- rbind(GOcns, cns3)
KEGGcns <- rbind(cns4, cns5)
KEGGcns <- rbind(KEGGcns, cns6)
KEGGcv <- rbind(cv4, cv5)
KEGGcv <- rbind(KEGGcv, cv6)
KEGGem <- rbind(em4, em5)
KEGGem <- rbind(KEGGem, em6)
GO <- rbind(dat1, dat2)
GO <- rbind(GO, dat3)
KEGG <- rbind(dat4, dat5)
KEGG <- rbind(KEGG, dat6)
GO_KEGG <- rbind(GO, KEGG)

#関数作成
split_into_multiple <- function(column, pattern = ", ", into_prefix){
  cols <- str_split_fixed(column, pattern, n = Inf)
  # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
  cols[which(cols == "")] <- NA
  cols <- as.tibble(cols)
  # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m' 
  # where m = # columns of 'cols'
  m <- dim(cols)[2]
  names(cols) <- paste(into_prefix, 1:m, sep = "_")
  return(cols)
}

#gene list抽出
dat1 <- dat1 %>% 
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "gene1")) %>% 
  # selecting those that start with 'gene_' will remove the original 'AssociatedGenesFound' column
  # select(AssociatedGenesFound, starts_with("type_"))
  gather(key, value = gene1, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  #select(gene) %>% #gene listのみ抽出, 省略可能
  distinct(gene1,.keep_all=FALSE) #gene重複削除しgeneのみ返す
dat2 <- dat2 %>% 
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "gene2")) %>% 
  gather(key, value = gene2, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(gene2,.keep_all=FALSE) #gene重複削除しgeneのみ返す
dat3 <- dat3 %>% 
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "gene3")) %>% 
  gather(key, value = gene3, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(gene3,.keep_all=FALSE) #gene重複削除しgeneのみ返す
dat4 <- dat4 %>% 
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "gene4")) %>% 
  gather(key, value = gene4, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(gene4,.keep_all=FALSE) #gene重複削除しgeneのみ返す
dat5 <- dat5 %>% 
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "gene5")) %>% 
  gather(key, value = gene5, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(gene5,.keep_all=FALSE) #gene重複削除しgeneのみ返す
dat6 <- dat6 %>% 
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "gene6")) %>% 
  gather(key, value = gene6, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(gene6,.keep_all=FALSE) #gene重複削除しgeneのみ返す

GOcns <- GOcns %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "GOcnsgene")) %>% 
  gather(key, value = GOcnsgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(GOcnsgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す
KEGGcns <- KEGGcns %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "KEGGcnsgene")) %>% 
  gather(key, value = KEGGcnsgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(KEGGcnsgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す
KEGGcv <- KEGGcv %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "KEGGcvgene")) %>% 
  gather(key, value = KEGGcvgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(KEGGcvgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す
KEGGem <- KEGGem %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "KEGGcvgene")) %>% 
  gather(key, value = KEGGemgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(KEGGemgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す
GO <- GO %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "GOgene")) %>% 
  gather(key, value = GOgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(GOgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す
KEGG <- KEGG %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "KEGGgene")) %>% 
  gather(key, value = KEGGgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(KEGGgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す
GO_KEGG <- GO_KEGG %>%
  bind_cols(split_into_multiple(.$AssociatedGenesFound, ", ", "GO_KEGGgene")) %>% 
  gather(key, value = GO_KEGGgene, -AssociatedGenesFound, na.rm = T) %>% #geneを1列にまとめる
  distinct(GO_KEGGgene,.keep_all=FALSE) #gene重複削除しgeneのみ返す

#通し番号列追加
dat1 <- dat1 %>% mutate(num = row_number())
dat2 <- dat2 %>% mutate(num = row_number())
dat3 <- dat3 %>% mutate(num = row_number())
dat4 <- dat4 %>% mutate(num = row_number())
dat5 <- dat5 %>% mutate(num = row_number())
dat6 <- dat6 %>% mutate(num = row_number())
GOcns <- GOcns %>% mutate(num = row_number())
KEGGcns <- KEGGcns %>% mutate(num = row_number())
KEGGcv <- KEGGcv %>% mutate(num = row_number())
KEGGem <- KEGGem %>% mutate(num = row_number())
GO <- GO %>% mutate(num = row_number())
KEGG <- KEGG %>% mutate(num = row_number())
GO_KEGG <- GO_KEGG %>% mutate(num = row_number())

#結合
m <- merge(dat1,dat2, by = "num", all = T)
m <- merge(m, dat3, by = "num", all = T)
m <- merge(m, dat4, by = "num", all = T)
m <- merge(m, dat5, by = "num", all = T)
m <- merge(m, dat6, by = "num", all = T)
m <- merge(m, GOcns, by = "num", all = T)
m <- merge(m, KEGGcns, by = "num", all = T)
m <- merge(m, KEGGcv, by = "num", all = T)
m <- merge(m, KEGGem, by = "num", all = T)
m <- merge(m, GO, by = "num", all = T)
m <- merge(m, KEGG, by = "num", all = T)
m <- merge(m, GO_KEGG, by = "num", all = T)
#積集合
t(names(m)) #行列の名前確認
GOint <- intersect(m$gene1, m$gene2)
GOint <- intersect(GOint, m$gene3)
KEGGint <- intersect(m$gene4, m$gene5)
KEGGint <- intersect(KEGGint, m$gene6)
#データフレーム変換
GOint <- data.frame(GOint=GOint) #通し番号列追加
GOint <- GOint %>%
  distinct(GOint,.keep_all=FALSE) #gene重複削除しgeneのみ返す
GOint <- GOint %>% mutate(num = row_number())
KEGGint <- data.frame(KEGGint=KEGGint) #通し番号列追加
KEGGint <- KEGGint %>%
  distinct(KEGGint,.keep_all=FALSE) #gene重複削除しgeneのみ返す
KEGGint <- KEGGint %>% mutate(num = row_number()) #通し番号列追加
m <- merge(m, GOint, by = "num", all = T)
m <- merge(m, KEGGint, by = "num", all = T)
#出力
write.csv(m, file="GO_KEGG_list.csv", na="", quote=F, row.names=F)

#############################################################
#置換(うまくいかない)
#data1$`Associated Genes Found`[data1$`Associated Genes Found` == "["] <- ""
#data1$`Associated Genes Found`[data1$`Associated Genes Found` == "]"] <- ""
#分割
#split <- str_split(data1$`Associated Genes Found`, pattern = ",", simplify = TRUE)
#############################################################
#エクセル出力
#smp <- list("GO_P"=data1_P, "GO_C"=data1_C, "GO_I"=data1_I,
#            "KEGG_P"=data2_P, "KEGG_C"=data2_C, "KEGG_I"=data2_I) #GOリスト,KEGGリスト作成
#write.xlsx(smp, "GO_KEGG_top3.xlsx", rowNames = F) #GOシート,KEGGシート出力
#write.xlsx(data1_P, "data1_P.xlsx", sheetName="sheet1") #出力