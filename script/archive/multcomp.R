#multcomp
#https://astatsa.com/OneWay_Anova_with_TukeyHSD/_Rcode_tutorial/
#エクセル入力
#http://www.yujitakenoshita.org/post/read-excel-in-r/
#列削除
#http://byungdugjun.blogspot.com/2014/07/r-x-image-disease-1-1-1-2-2-0-3-0-1-4-1.html
#BH法
#https://stats.biopapyrus.jp/stats/fdr-bh.html
#THSD
#https://qiita.com/hfu62/items/f9f4803828fd7e1a5cec
#My anova loop prints out the same results in R
#https://stackoverflow.com/questions/50914023/my-anova-loop-prints-out-the-same-results-in-r
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/multcomp") #作業ディレクトリ設定
#setwd("~/GoogleDrive/マイドライブ/0_Work/R/SWATH") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
##############################################################
rm(list = ls(all.names = TRUE))
#パッケージインストール
#install.packages("multcomp") #多重比較検定
#install.packages("tidyverse")
#install.packages("dplyr")
#ライブラリ読み込み
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
#library(dplyr)
library(readxl) #エクセル入力(read_excel)
library(openxlsx) #エクセル入出力(write.xlsx)
library(xlsx) #エクセル入力
library(multcomp) #多重比較検定
library(BH) #FDR
##############################################################
#xlsx入力
#rm(list = ls(all.names = TRUE))
library(readxl) #エクセル入力(read_excel) #エラーが出るため再読み込み
data <- read_excel("data4-6.xlsx", 2) #シート2入力
#data <- read.xlsx("data4-6.xlsx", 2, row.names=1) #シート2入力,row.names, 動くときと動かない時がある
list <- data$row.names
data_rm  <- data
data_rm[,1:3] <- NULL #列削除
rownames(data_rm) <- list
#############################################################
PC <- gl(6,12,72, label=c("aSC0", "aSC10", "aSC30", "bPC0", "bPC10", "bPC30")) #カテゴリー6,繰返し12,要素72,ラベル6種
PCP <- gl(2,36,72, label=c("PCP", "SAL")) #PCP/SAL 2,繰返し36,要素72(2x36リピート),ラベル2種
CLZ <- gl(3,12,72, label=c("CLZ0", "CLZ10", "CLZ30")) #CLZ 3doses,繰返し12,要素72,ラベル3種
#############################################################
#1wANOVA function
#############################################################
aof <- function(x) { 
  m <- data.frame(PC, x); 
  anova(aov(x ~ PC, m))
}
# apply analysis to the data and get the pvalues.
onewayANOVA <- apply(data_rm, 1, aof)
onewayANOVAp <- data.frame(lapply(onewayANOVA, function(x) { x["Pr(>F)"][1,] }))
onewayANOVAp2 <- data.frame(t(onewayANOVAp))
colnames(onewayANOVAp2) <- "p_PC" #列名の変更
write.table(onewayANOVAp2, file="1wanova-results.txt", quote=F, sep='\t')
#onewayANOVAp3 <- read.delim2("1wanova-results.txt") #以下でsdataに統合
#onewayANOVAp3 <- cbind(onewayANOVAp3, rownames(onewayANOVAp3)) #以下でsdataに統合
#############################################################
#2wANOVA function
#############################################################
aof2 <- function(x) { 
  n <- data.frame(PCP,CLZ, x); 
  anova(aov(x ~ PCP + CLZ + PCP*CLZ, n))
}
# apply analysis to the data and get the pvalues.
twowayANOVA <- apply(data_rm, 1, aof2)
twowayANOVAp <- data.frame(lapply(twowayANOVA, function(x) { x["Pr(>F)"][1:3,] }))
twowayANOVAp2 <- data.frame(t(twowayANOVAp))
colnames(twowayANOVAp2) <- c("p_PCP","p_CLZ","p_interaction") #列名の取得
sapply(twowayANOVAp2, class)
write.table(twowayANOVAp2, file="2wanova-results.txt", quote=F, sep='\t')
#twowayANOVAp3 <- read.delim2("2wanova-results.txt") #以下でsdataに統合
#sapply(twowayANOVAp3, class)
sdata <- cbind(data_rm,c(onewayANOVAp2, twowayANOVAp2))
sapply(sdata, class) #p値がnumericであることを確認
#write.table(sdata, file="anova-results.txt", quote=F, sep='\t')
#############################################################
#2wANOVA BH-FDR
#############################################################
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み,#エラーが出るため再読み込み
#sort <- sdata %>% dplyr::arrange(p_PCP) #データフレーム上でソート(increasing order),FDR計算に使用できない
#sort_list <- list(sort$p_PCP) #p値のみリスト化,FDR計算に使用できない
#PCPp_sort <- sort(sdata$p_PCP)
#sapply(PCPp_sort, class) #p値がnumericであることを確認
#dfではBH-FDR計算できない
#PCPp <- data.frame(sort(t(twowayANOVAp[1,]))) # the p values must be sorted in increasing order! 
#CLZp <- data.frame(sort(t(twowayANOVAp[2,]))) # the p values must be sorted in increasing order! 
#PCPCLZp <- data.frame(sort(t(twowayANOVAp[3,]))) # the p values must be sorted in increasing order! 
#dfにしないとBH-FDR計算できるが、昇順にしなくても良い!!
#PCPp <- sort(t(twowayANOVAp[1,])) # the p values must be sorted in increasing order! 
#CLZp <- sort(t(twowayANOVAp[2,])) # the p values must be sorted in increasing order! 
#PCPCLZp <- sort(t(twowayANOVAp[3,])) # the p values must be sorted in increasing order! 
PCp <- sdata$p_PC
PCPp <- sdata$p_PCP 
CLZp <- sdata$p_CLZ
PCPCLZp <- sdata$p_interaction
checkP <- data.frame(cbind(PCp, PCPp, CLZp, PCPCLZp))

PCq <- data.frame(p.adjust(PCp, method = "BH"))
PCPq <- data.frame(p.adjust(PCPp, method = "BH"))
CLZq <- data.frame(p.adjust(CLZp, method = "BH"))
PCPCLZq <- data.frame(p.adjust(PCPCLZp, method = "BH"))
checkQ <- data.frame(cbind(PCq, PCPq, CLZq, PCPCLZq))

#fdr.result <- bh(PCPp, 0.05) #bh()関数が動かない
bhthresh <- cbind(PCq, PCPq, CLZq, PCPCLZq) # put the bh results in our table. 
write.table(bhthresh, "bhthresh.txt", sep='\t', quote=F) # print to a file.
sdata <- cbind(sdata, bhthresh)
sapply(sdata, class) #p,q値がnumericであることを確認
write.table(sdata, file="anova-results.txt", quote=F, sep='\t')
#############################################################
#############################################################
#############################################################
# Get the genes with good PCP effect pvalues.
#reg.hi.pdata <- filter(sdata, p_PCP < 0.0001 & p_interaction > 0.1)
#reg.hi.p <- sort(t(data.frame(pvalues3[1, pvalues3[1,] < 0.05 & pvalues3[3,] > 0.1]))) #上手くいかない→dplyrで実施
#reg.hi.pdata <- data_rm[row.names(reg.hi.p),] #上手くいかない→dplyrで実施
#############################################################
#TukeyHSD function
#diff群間の平均値の差(例)B-Aが-127.3であればデータBの平均がデータAの平均より-127.3大きい
#lwr,upr=下方信頼限界,情報信頼限界:信頼区間の下限値 (lower) と上限値 (upper)
#0を含まない場合 (例)B-A は含まず D-A は含む=2群間差は0ではないので有意差あり
#p.adj < 0.05=2群間に有意差あり(信頼区間内に0を含まない)
#############################################################
#THSD1w <- function(x) { 
#  nnn <- data.frame(PC, x); 
#  TukeyHSD(aov(x ~ PC, nnn))
#}
#THSDresults1w <- apply(data_rm, 1, THSD1w) 
THSD <- function(x) { 
  nn <- data.frame(PCP,CLZ, x); 
  TukeyHSD(aov(x ~ PCP + CLZ + PCP*CLZ, nn))
}
THSDresults <- apply(data_rm, 1, THSD) 
#plot(THSDresults[["SPTN1"]]) #SPTN1のプロット
#xx <- data.frame(THSDresults[["SPTN1"]][["PCP:CLZ"]])$`p.adj`
#xx$`p.adj`
THSD_PCP <- data.frame(lapply(THSDresults, function(x) {x["PCP"]}))
THSD_CLZ <- data.frame(lapply(THSDresults, function(x) {x["CLZ"]}))
THSD_PCPCLZ <- data.frame(lapply(THSDresults, function(x) {x["PCP:CLZ"]}))
write.table(THSD_PCP, file="THSD_PCP.txt", quote=F, sep='\t')
write.table(THSD_CLZ, file="THSD_CLZ.txt", quote=F, sep='\t')
write.table(THSD_PCPCLZ, file="THSD_PCPCLZ.txt", quote=F, sep='\t')
#THSD_PCP2 <- read.table("THSD_PCP.txt",header=T, row.names=1)
#THSD_CLZ2 <- read.table("THSD_CLZ.txt",header=T, row.names=1)
#THSD_PCPCLZ2 <- read.table("THSD_PCPCLZ.txt",header=T, row.names=1)
#エラーが出るのでtidyverse再インストール,再読み込み
install.packages("tidyverse")
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
THSDp_PCP <- select(THSD_PCP, ends_with("p.adj")) #p値抽出
THSDp_CLZ <- select(THSD_CLZ, ends_with("p.adj")) #p値抽出
THSDp_PCPCLZ <- select(THSD_PCPCLZ, ends_with("p.adj")) #p値抽出
#transpose
THSDp_PCP2 <- data.frame(t(THSDp_PCP))
#THSDp_PCP2 <- cbind(THSDp_PCP2, rownames(THSDp_PCP2)) #以下でsdataに統合
THSDp_CLZ2 <- data.frame(t(THSDp_CLZ))
#THSDp_CLZ2 <- cbind(THSDp_CLZ2, rownames(THSDp_CLZ2)) #以下でsdataに統合
THSDp_PCPCLZ2 <- data.frame(t(THSDp_PCPCLZ))
#THSDp_PCPCLZ2 <- cbind(THSDp_PCPCLZ2, rownames(THSDp_PCPCLZ2)) #以下でsdataに統合
sdata <- cbind(sdata,c(THSDp_PCP2, THSDp_CLZ2, THSDp_PCPCLZ2))
#############################################################
#xlsx出力
#############################################################
#シートに分かれない!!!
#library(openxlsx)
#smp <- list("sheet1" = sdata, "sheet2" = checkP, "sheet3" = checkQ) #リスト作成
#write.xlsx(smp, "stat.xlsx") #1-3シート出力
library(writexl)
sdata2 <- cbind(rownames(sdata),sdata)
sheets <- list("sheet1" = sdata2, "sheet2" = checkP, "sheet3" = checkQ) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets, "stat2.xlsx", format_headers = FALSE)
#############################################################
#############################################################
#############################################################
# correlation:検討中
#############################################################
#比較するグループ(列)の指定
PCPtempl <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) #PCP
#cor.test(t(data_rm[14,]),PCPtempl) #14行目の相関解析
template.match <- function(x, template) {
  k <- cor.test(x,template)
  k$p.value
}
#t(names(reg.hi.pdata)) #行列の名前確認
#reg.hi.pdata_rm <- reg.hi.pdata
#t(names(reg.hi.pdata_rm)) #行列の名前確認
#reg.hi.pdata_rm[,73:75] <- NULL #列削除
#reg.hi.pdata_rm

#Filtered
#PCPtempl.results <- apply(reg.hi.pdata_rm, 1, template.match, PCPtempl) 
#PCPtempl.results2 <- cbind(reg.hi.pdata, PCPtempl.results)

#All
PCPtempl.results3 <- apply(data_rm, 1, template.match, PCPtempl) #PCPtempleにて2群の設定が必要
PCPtempl.results4 <- cbind(sdata, PCPtempl.results3)
#PCPtempl.results4_filt <- filter(PCPtempl.results4, p_PCP < 0.0001 & p_interaction > 0.1)
names(PCPtempl.results4)[which(names(PCPtempl.results4)=="PCPtempl.results3")] <- "p_cortest"
write.table(PCPtempl.results4, file="PCPtempl-results.txt", quote=F, sep='\t')
d <- read.delim2("PCPtempl-results.txt")
#############################################################
# ttest:検討中
#############################################################
#t.test(t(reg.hi.pdata[1,PCP!="PCP"]), t(reg.hi.pdata[1,PCP=="SAL"]), var.equal=T)
template.match2 <- function(x, template) {
  l <- t.test(x,template)
  l$p.value
}
PCPtempl.results5 <- data.frame(apply(data_rm, 1, template.match2, PCPtempl))
colnames(PCPtempl.results5) <- "p_ttest_PCP" #列名の変更
#PCPtempl.results6 <- cbind(sdata, PCPtempl.results5)
#names(PCPtempl.results6)[which(names(PCPtempl.results6)=="PCPtempl.results5")] <- "p_ttest_PCP"
#############################################################
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
#以下,不要
############################################################
############################################################
############################################################
#xlsx入力
#rm(list = ls(all.names = TRUE))
#expt1 <- read.table("demo.txt",header=T)
data <- read_excel("data4-6.xlsx", 2) #シート2入力
name <- as.list(t(data[,1])) # first remember the names
data_t <- as.data.frame(t(data[,-(1:3)])) # transpose all but the first three columns
colnames(data_t) <- name
str(data_t) # Check the column types
#############################################################
#txt出力
write.table (data_t, file = "data_t.txt", sep = "\t",quote = FALSE, row.names = TRUE)
#txt入力
data_tt <- read.table("data_t.txt",header=T, row.names=NULL)
#xlsx出力
#smp <- list("data4"=data4,"data5"=data5,"data6"=data6) #リスト作成,Subtract後,Subtract前,元データ+Subtract前
#write.xlsx(smp, "data4-6.xlsx") #シート出力
#############################################################
#データ整形
split <- str_split(data_tt[,1], pattern = "\\_", simplify = TRUE) #_で分割
colnames(split) <- c("group", "no") #列名変更
y <- data.frame(split)
#置換
y$group <- sub("PC","bPC", y$group) #sortのためPCの接頭語b挿入
y$group <- sub("SC","aSC", y$group) #sortのためSCの接頭語a挿入
#as.factor(y$group)
#as.data.frame(y)
data <- cbind(data_tt, y) #dataとxを列ベクトル単位で結合
#as.factor(data$Group)
#############################################################
#1wANOVA
#amod <- aov(data[,2]~data$group,data=data)
#amod2 <- aov(data$SPTN1~data$group,data=data)
amod <- aov(SPTN1~group,data=data)
#amod2 <- aov(SPTN1~group,data=data)
summary(amod)
#############################################################
#Tukey HSD multiple comparison
#workしない
#data$group = as.factor(data$group)
#data$Group
#class(data$group)
tmod <- glht(amod,linfct = mcp(group = "Tukey"))
summary(tmod)
TukeyHSD(amod, conf.level = 0.95)
#TukeyHSD(amod, conf.level = 0.99)

T <- summary(tmod)$test$tstat
as.data.frame(T)

k <- amod$rank
v <- amod$df.residual
pValSheffe <- 1-pf(T**2/(k-1),k-1,v)
as.matrix(pValSheffe)
#############################################################
#Bonferroni simultaneous multiple comparison of q pairs
data$group #aSC0 aSC10 aSC30 bPC0 bPC10 bPC30
contrasts <- rbind(
  "aSC10 - aSC0"  = c(-1,1,0,0,0,0),
  "aSC30 - aSC0"  = c(-1,0,1,0,0,0),
  "bPC0 - aSC0"   = c(-1,0,0,1,0,0),
  "bPC10 - aSC0"  = c(-1,0,0,0,1,0),
  "bPC30 - aSC0"  = c(-1,0,0,0,0,1),
  "aSC30 - aSC10" = c(0,-1,1,0,0,0),
  "bPC0 - aSC10"  = c(0,-1,0,1,0,0),
  "bPC10 - aSC10" = c(0,-1,0,0,1,0),
  "bPC30 - aSC10" = c(0,-1,0,0,0,1),
  "bPC0 - aSC30"  = c(0,0,-1,1,0,0),
  "bPC10 - aSC30" = c(0,0,-1,0,1,0),
  "bPC30 - aSC30" = c(0,0,-1,0,0,1),
  "bPC10 - bPC0"  = c(0,0,0,-1,1,0),
  "bPC30 - bPC0"  = c(0,0,0,-1,0,1),
  "bPC30 - bPC10" = c(0,0,0,0,-1,1)
  )
contrasts
bmod <- glht(amod,linfct = mcp(group = contrasts))
summary(bmod,test = adjusted("bonferroni"))
#############################################################
#Holm simultaneous multiple comparison for all pairs
summary(bmod,test=adjusted("holm"))
#############################################################
#Bonferroni simultaneous multiple comparison for the relevant q=5 pairs of contrasts relative to aSC0 only
contrasts2 <- rbind(
  "aSC10 - aSC0"  = c(-1,1,0,0,0,0),
  "aSC30 - aSC0"  = c(-1,0,1,0,0,0),
  "bPC0 - aSC0"   = c(-1,0,0,1,0,0),
  "bPC10 - aSC0"  = c(-1,0,0,0,1,0),
  "bPC30 - aSC0"  = c(-1,0,0,0,0,1)
  )
contrasts2
bmod2 <- glht(amod,linfct=mcp(group = contrasts2))
summary(bmod2,test=adjusted("bonferroni")) 
summary(bmod2,test=adjusted("holm")) 
#############################################################
z <- as.data.frame(pValSheffe)
z <- cbind(z,rownames(z))
z <- cbind(z,contrasts)
rownames(z) <- NULL
colnames(z) <- c("pValSheffe","contrasts","aSC0", "aSC10", "aSC30", "bPC0", "bPC10", "bPC30")
z[,1] <- rownames(z)
names(z)[1] <- "No"
zz <- z[,-1]
#############################################################
#txt出力
write.table (z, file = "z.txt", sep = "\t",quote = FALSE, row.names = FALSE)
write.table (zz, file = "zz.txt", sep = "\t",quote = FALSE, row.names = FALSE)
#csv出力
write.csv(z, file = "z.csv", row.names = FALSE)
write.csv(zz, file = "zz.csv", row.names = FALSE)
#txt入力
z2 <- read.table("z.txt",header=T)#, row.names=NULL
zz2 <- read.table("zz.txt",header=T)
#csv入力
z3 <- read.csv("z.csv",row.names=1)
zz3 <- read.csv("zz.csv",row.names=1)
#xlsx出力
smp <- list("z"=z,"zz"=zz) #z通し番号あり,zz通し番号なし
write.xlsx(smp, "z-zz.xlsx") #シート出力
#xlsx入力
z4 <- read_excel("z-zz.xlsx", 1) #シート1入力
zz4 <- read_excel("z-zz.xlsx", 2) #シート2入力
#############################################################
