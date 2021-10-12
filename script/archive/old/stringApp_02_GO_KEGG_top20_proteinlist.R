#stringApp_02_GO_KEGG_top20_proteinlist.R
#5_GO_KEGG_top20に含まれるDEP list作成
##############################################################
#setwd("/Users/user/Dropbox/0_Work/R/SWATH") #作業ディレクトリ設定
setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
options(digits=2) #change digit2桁表示指定
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
##############################################################
###         #5_GO_KEGG_top20に含まれるDEP list作成         ###
##############################################################
#xlsx入力
data1 <- read_excel("GO_KEGGrep.xlsx", 1) #シート1入力
data2 <- read_excel("GO_KEGGrep.xlsx", 2) #シート2入力
#列名の置換
t(names(data1)) #行列の名前確認
names(data1)[which(names(data1)=="genes.y.y") ] <- "genes"
names(data2)[which(names(data2)=="genes.y.y") ] <- "genes"
##############################################################
#各群のデータフレーム作成
dat7 <- filter(data1, group2 == "PCP") #PCPのtop20,NA除外：!is.na(`genes`)
dat8 <- filter(data1, group2 == "CLZ") #CLZのtop20
dat9 <- filter(data1, group2 == "Interaction") #Interactionのtop20
dat10 <- filter(data2, group2 == "PCP") #PCPのtop20
dat11 <- filter(data2, group2 == "CLZ") #CLZのtop20
dat12 <- filter(data2, group2 == "Interaction") #Interactionのtop20
##############################################################
#行の抽出(CNS, Signaling,EndcrineMetabolism)
dat7[,5] #5列目(description.x.x)表示
dat8[,5] #5列目(description.x.x)表示
dat9[,5] #5列目(description.x.x)表示
dat10[,5] #5列目(description.x.x)表示
dat11[,5] #5列目(description.x.x)表示
dat12[,5] #5列目(description.x.x)表示
index1 <- c(1,2,9,14,16) #CNSrelated
index2 <- c(2,8,9,12,16) #CNSrelated
index3 <- c(6,18,19) #Signalrelated
index4 <- c(1,3,4,5,13,15,17,20) #EMrelated
##############################################################
t(names(dat7)) #行列の名前確認
t(names(dat8)) #行列の名前確認
t(names(dat9)) #行列の名前確認
t(names(dat10)) #行列の名前確認
t(names(dat11)) #行列の名前確認
t(names(dat12)) #行列の名前確認
#CNSrelated
cns7 <- dat7 %>% slice(index1) %>% select(51)
cns8 <- dat8 %>% slice(index1) %>% select(51)
cns9 <- dat9 %>% slice(index1) %>% select(51)
cns10 <- dat10 %>% slice(index2) %>% select(51)
cns11 <- dat11 %>% slice(index2) %>% select(51)
cns12 <- dat12 %>% slice(index2) %>% select(51)
#Signalrelated
sgn10 <- dat10 %>% slice(index3) %>% select(51)
sgn11 <- dat11 %>% slice(index3) %>% select(51)
sgn12 <- dat12 %>% slice(index3) %>% select(51)
#EMrelated
em10 <- dat10 %>% slice(index4) %>% select(51)
em11 <- dat11 %>% slice(index4) %>% select(51)
em12 <- dat12 %>% slice(index4) %>% select(51)
#列抽出
t(names(dat7)) #行列の名前確認
dat7 <- dat7 %>% select(51)
dat8 <- dat8 %>% select(51)
dat9 <- dat9 %>% select(51)
dat10 <- dat10 %>% select(51)
dat11 <- dat11 %>% select(51)
dat12 <- dat12 %>% select(51)
#縦に結合
GOcns <- rbind(cns7, cns8)
GOcns <- rbind(GOcns, cns9)
KEGGcns <- rbind(cns10, cns11)
KEGGcns <- rbind(KEGGcns, cns12)
KEGGsgn<- rbind(sgn10, sgn11)
KEGGsgn <- rbind(KEGGsgn, sgn12)
KEGGem <- rbind(em10, em11)
KEGGem <- rbind(KEGGem, em12)
GO <- rbind(dat7, dat8)
GO <- rbind(GO, dat9)
KEGG <- rbind(dat10, dat11)
KEGG <- rbind(KEGG, dat12)
GO_KEGG <- rbind(GO, KEGG)
##############################################################
#関数作成
split_into_multiple <- function(column, pattern = ",", into_prefix){
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
##############################################################
#gene list抽出
dat7 <- dat7 %>% 
  bind_cols(split_into_multiple(.$genes, ",", "gene1")) %>% 
  # selecting those that start with 'gene_' will remove the original 'AssociatedGenesFound' column
  # select(AssociatedGenesFound, starts_with("type_"))
  gather(key, value = gene1, -genes, na.rm = T) %>% #geneを1列にまとめる
  #select(gene) %>% #gene listのみ抽出, 省略可能
  distinct(gene1,.keep_all=FALSE) #gene重複削除しgeneのみ返す
dat8 <- dat8 %>% 
  bind_cols(split_into_multiple(.$genes, ",", "gene2")) %>%
  gather(key, value = gene2, -genes, na.rm = T) %>%
  distinct(gene2,.keep_all=FALSE)
dat9 <- dat9 %>% 
  bind_cols(split_into_multiple(.$genes, ",", "gene3")) %>%
  gather(key, value = gene3, -genes, na.rm = T) %>%
  distinct(gene3,.keep_all=FALSE)
dat10 <- dat10 %>% 
  bind_cols(split_into_multiple(.$genes, ",", "gene4")) %>%
  gather(key, value = gene4, -genes, na.rm = T) %>%
  distinct(gene4,.keep_all=FALSE)
dat11 <- dat11 %>% 
  bind_cols(split_into_multiple(.$genes, ",", "gene5")) %>%
  gather(key, value = gene5, -genes, na.rm = T) %>%
  distinct(gene5,.keep_all=FALSE)
dat12 <- dat12 %>% 
  bind_cols(split_into_multiple(.$genes, ",", "gene6")) %>%
  gather(key, value = gene6, -genes, na.rm = T) %>%
  distinct(gene6,.keep_all=FALSE)
##############################################################
GOcns <- GOcns %>%
  bind_cols(split_into_multiple(.$genes, ",", "GOcnsgene")) %>% 
  gather(key, value = GOcnsgene, -genes, na.rm = T) %>%
  distinct(GOcnsgene,.keep_all=FALSE)
KEGGcns <- KEGGcns %>%
  bind_cols(split_into_multiple(.$genes, ",", "KEGGcnsgene")) %>% 
  gather(key, value = KEGGcnsgene, -genes, na.rm = T) %>%
  distinct(KEGGcnsgene,.keep_all=FALSE)
KEGGsgn <- KEGGsgn %>%
  bind_cols(split_into_multiple(.$genes, ",", "KEGGsgngene")) %>% 
  gather(key, value = KEGGsgngene, -genes, na.rm = T) %>%
  distinct(KEGGsgngene,.keep_all=FALSE)
KEGGem <- KEGGem %>%
  bind_cols(split_into_multiple(.$genes, ",", "KEGGcvgene")) %>% 
  gather(key, value = KEGGemgene, -genes, na.rm = T) %>%
  distinct(KEGGemgene,.keep_all=FALSE)
GO <- GO %>%
  bind_cols(split_into_multiple(.$genes, ",", "GOgene")) %>% 
  gather(key, value = GOgene, -genes, na.rm = T) %>%
  distinct(GOgene,.keep_all=FALSE)
KEGG <- KEGG %>%
  bind_cols(split_into_multiple(.$genes, ",", "KEGGgene")) %>% 
  gather(key, value = KEGGgene, -genes, na.rm = T) %>%
  distinct(KEGGgene,.keep_all=FALSE)
GO_KEGG <- GO_KEGG %>%
  bind_cols(split_into_multiple(.$genes, ",", "GO_KEGGgene")) %>% 
  gather(key, value = GO_KEGGgene, -genes, na.rm = T) %>%
  distinct(GO_KEGGgene,.keep_all=FALSE)
##############################################################
#通し番号列追加
dat1 <- dat7 %>% mutate(num = row_number())
dat2 <- dat8 %>% mutate(num = row_number())
dat3 <- dat9 %>% mutate(num = row_number())
dat4 <- dat10 %>% mutate(num = row_number())
dat5 <- dat11 %>% mutate(num = row_number())
dat6 <- dat12 %>% mutate(num = row_number())
GOcns <- GOcns %>% mutate(num = row_number())
KEGGcns <- KEGGcns %>% mutate(num = row_number())
KEGGsgn <- KEGGsgn %>% mutate(num = row_number())
KEGGem <- KEGGem %>% mutate(num = row_number())
GO <- GO %>% mutate(num = row_number())
KEGG <- KEGG %>% mutate(num = row_number())
GO_KEGG <- GO_KEGG %>% mutate(num = row_number())
##############################################################
#結合
m <- merge(dat1,dat2, by = "num", all = T)
m <- merge(m, dat3, by = "num", all = T)
m <- merge(m, dat4, by = "num", all = T)
m <- merge(m, dat5, by = "num", all = T)
m <- merge(m, dat6, by = "num", all = T)
m <- merge(m, GOcns, by = "num", all = T)
m <- merge(m, KEGGcns, by = "num", all = T)
m <- merge(m, KEGGsgn, by = "num", all = T)
m <- merge(m, KEGGem, by = "num", all = T)
m <- merge(m, GO, by = "num", all = T)
m <- merge(m, KEGG, by = "num", all = T)
m <- merge(m, GO_KEGG, by = "num", all = T)
##############################################################
#積集合
t(names(m)) #行列の名前確認
GOint <- intersect(m$gene1, m$gene2)
GOint <- intersect(GOint, m$gene3)
KEGGint <- intersect(m$gene4, m$gene5)
KEGGint <- intersect(KEGGint, m$gene6)
#データフレーム変換
GOint <- data.frame(GOint=GOint)
GOint <- GOint %>% distinct(GOint,.keep_all=FALSE) #gene重複削除しgeneのみ返す
GOint <- GOint %>% mutate(num = row_number()) #通し番号列追加
KEGGint <- data.frame(KEGGint=KEGGint)
KEGGint <- KEGGint %>% distinct(KEGGint,.keep_all=FALSE) #gene重複削除しgeneのみ返す
KEGGint <- KEGGint %>% mutate(num = row_number()) #通し番号列追加
m <- merge(m, GOint, by = "num", all = T)
m <- merge(m, KEGGint, by = "num", all = T)
##############################################################
#csv出力
write.csv(m, file="GO_KEGG_top20_proteinlist.csv", na="", quote=F, row.names=F)
##############################################################
