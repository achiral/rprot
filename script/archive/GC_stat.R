#GC_stat
#基本統計量
#https://note.com/kotoko_tyo/n/ne09e4398093d
#カテゴリカルデータ集計
#https://www1.doshisha.ac.jp/~mjin/R/Chap_45/45.html
#http://www.eeso.ges.kyoto-u.ac.jp/emm/wp-content/uploads/2010/11/questionary01.pdf
#http://www.housecat442.com/?p=346
#正規性####
#シャピロウィルク検定(Shapiro-Wilk test) 
#Kolmogorov-Smirnov test
#https://data-science.gr.jp/implementation/ist_r_kolmogorov_smirnov_test.html
#等分散性####
#Bartlett test
#https://data-science.gr.jp/implementation/ist_r_bartlett_test.html
#Levene's test
#https://data-science.gr.jp/implementation/ist_r_levene_test.html
#Hartley test
#https://data-science.gr.jp/implementation/ist_r_hartley_test.html
#推奨はバートレット検定,非正規分布はルビーン検定,高検出力が必要な時はハートレイ検定

#ノンパラ検定
#https://data-science.gr.jp/implementation/ist_r_kruskal_wallis_test.html
#http://mizumot.com/handbook/?page_id=422
#Steel-Dwass（スティールドゥワス）検定
#https://jojoshin.hatenablog.com/entry/2016/05/29/222528

#Fisher's exact test
#https://data-science.gr.jp/implementation/ist_r_fisher_exact_probability_test.html
#Pearson’s chi-squared test
#https://www.datascienceblog.net/post/statistical_test/contingency_table_tests/
#カプランマイヤー曲線
#https://note.com/maxwell/n/nc3ab385128d8
#正規分布に従うものか否かを調べる検定法
#https://data-science.gr.jp/implementation/ist_r_shapiro_wilk_test.html
#回帰分析
#############################################################
rm(list = ls(all = TRUE))
setwd("/Users/user/Dropbox/My Mac (AkiranoMacBook.local)/Desktop/GC_stat")
getwd()
#############################################################
#Excelファイル読み込み
#install.packages("readxl")
#install.packages("psych")
library(readxl)
data1 <- read_excel("Book2.xlsx", 3)
data2 <- read_excel("Book2.xlsx", 4)
data3 <- read_excel("Book2.xlsx", 5)
data4 <- read_excel("Book2.xlsx", 6)
data5 <- read_excel("Book2.xlsx", 7)
data6 <- read_excel("Book2.xlsx", 8)
data7 <- read_excel("Book2.xlsx", 9)
data8 <- read_excel("Book2.xlsx", 10)
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
#library(openxlsx)
#write.xlsx(data4, sheetname = "sheet1", file = "out_file.xlsx") 
#############################################################
#変数のクラス確認####
sapply(data1, class)
sapply(data2, class)
sapply(data3, class)
#characterをfactorに変換
#data[,1] <- as.factor(data[,1])
#class(data[,1])
#基本統計量の一覧表####
library(psych)
(out1 <- describeBy(data1, group = data1$Regimen))
(out2 <- describeBy(data2, group = data2$Regimen))
(out3 <- describeBy(data3, group = data3$Regimen))

#ダミー変数加工####
#install.packages("caret")
library(caret)
library(ggplot2)
tmp <- dummyVars(~.-ID, data = data) #全質的変数を対象
data.dummy <- as.data.frame(predict(tmp, data)) #tmpの質的変数をダミー変数に変換
str(data.dummy) #Df構造表示


#Shapiro正規性の検定 #p<0.05:非正規→ノンパラ検定####
library(tidyverse)
#eGFR,Ccr,CIN,CIV(courses)
data1_GC <- data1 %>% filter(Regimen == "GC")
data1_GCa <- data1 %>% filter(Regimen == "GCa")
data1_GCs <- data1 %>% filter(Regimen == "GCs")
shapiro.test(x=data1$eGFR)
shapiro.test(x=data1_GC$eGFR) 
shapiro.test(x=data1_GCa$eGFR) 
shapiro.test(x=data1_GCs$eGFR) 
shapiro.test(x=data1$Ccr)
shapiro.test(x=data1_GC$Ccr) 
shapiro.test(x=data1_GCa$Ccr) 
shapiro.test(x=data1_GCs$Ccr) 
shapiro.test(x=data1$CIN)
shapiro.test(x=data1_GC$CIN) 
shapiro.test(x=data1_GCa$CIN) 
shapiro.test(x=data1_GCs$CIN) 
shapiro.test(x=data1$CIV)
shapiro.test(x=data1_GC$CIV) 
shapiro.test(x=data1_GCa$CIV) 
shapiro.test(x=data1_GCs$CIV) 
#Age(patients)
data2_GC <- data2 %>% filter(Regimen == "GC")
data2_GCa <- data2 %>% filter(Regimen == "GCa")
data2_GCs <- data2 %>% filter(Regimen == "GCs")
shapiro.test(x=data2$Age)
shapiro.test(x=data2_GC$Age)
shapiro.test(x=data2_GCa$Age)
shapiro.test(x=data2_GCs$Age)
#Survival(courses)
data3_GC <- data3 %>% filter(Regimen == "GC")
data3_GCa <- data3 %>% filter(Regimen == "GCa")
data3_GCs <- data3 %>% filter(Regimen == "GCs")
shapiro.test(x=data3$Time)
shapiro.test(x=data3_GC$Time)
shapiro.test(x=data3_GCa$Time)
shapiro.test(x=data3_GCs$Time)

#Kolmogorov-Smirnov test####
#Shapiro-Wilk test(N<2000), Kolmogorov-Smirnov test(N>=2000)???
#eGFR,Ccr,CIN,CIV(courses)
ks.test(x=data1$eGFR,y="pnorm",mean=mean(data1$eGFR),sd=sd(data1$eGFR))
ks.test(x=data1_GC$eGFR,y="pnorm",mean=mean(data1_GC$eGFR),sd=sd(data1_GC$eGFR))
ks.test(x=data1_GCa$eGFR,y="pnorm",mean=mean(data1_GCa$eGFR),sd=sd(data1_GCa$eGFR))
ks.test(x=data1_GCs$eGFR,y="pnorm",mean=mean(data1_GCs$eGFR),sd=sd(data1_GCs$eGFR))
ks.test(x=data1$Ccr,y="pnorm",mean=mean(data1$Ccr),sd=sd(data1$Ccr))
ks.test(x=data1_GC$Ccr,y="pnorm",mean=mean(data1_GC$Ccr),sd=sd(data1_GC$Ccr))
ks.test(x=data1_GCa$Ccr,y="pnorm",mean=mean(data1_GCa$Ccr),sd=sd(data1_GCa$Ccr))
ks.test(x=data1_GCs$Ccr,y="pnorm",mean=mean(data1_GCs$Ccr),sd=sd(data1_GCs$Ccr))
ks.test(x=data1$CIN,y="pnorm",mean=mean(data1$CIN),sd=sd(data1$CIN))
ks.test(x=data1_GC$CIN,y="pnorm",mean=mean(data1_GC$CIN),sd=sd(data1_GC$CIN))
ks.test(x=data1_GCa$CIN,y="pnorm",mean=mean(data1_GCa$CIN),sd=sd(data1_GCa$CIN))
ks.test(x=data1_GCs$CIN,y="pnorm",mean=mean(data1_GCs$CIN),sd=sd(data1_GCs$CIN))
ks.test(x=data1$CIV,y="pnorm",mean=mean(data1$CIV),sd=sd(data1$CIV))
ks.test(x=data1_GC$CIV,y="pnorm",mean=mean(data1_GC$CIV),sd=sd(data1_GC$CIV))
ks.test(x=data1_GCa$CIV,y="pnorm",mean=mean(data1_GCa$CIV),sd=sd(data1_GCa$CIV))
ks.test(x=data1_GCs$CIV,y="pnorm",mean=mean(data1_GCs$CIV),sd=sd(data1_GCs$CIV))
#Age(patients)
ks.test(x=data2$Age,y="pnorm",mean=mean(data2$Age),sd=sd(data2$Age))
ks.test(x=data2_GC$Age,y="pnorm",mean=mean(data2_GC$Age),sd=sd(data2_GC$Age))
ks.test(x=data2_GCa$Age,y="pnorm",mean=mean(data2_GCa$Age),sd=sd(data2_GCa$Age))
ks.test(x=data2_GCs$Age,y="pnorm",mean=mean(data2_GCs$Age),sd=sd(data2_GCs$Age))
#Survival(courses)
ks.test(x=data3$Time,y="pnorm",mean=mean(data3$Time),sd=sd(data3$Time))
ks.test(x=data3_GC$Time,y="pnorm",mean=mean(data3_GC$Time),sd=sd(data3_GC$Time))
ks.test(x=data3_GCa$Time,y="pnorm",mean=mean(data3_GCa$Time),sd=sd(data3_GCa$Time))
ks.test(x=data3_GCs$Time,y="pnorm",mean=mean(data3_GCs$Time),sd=sd(data3_GCs$Time))

#2群間の等分散性の検定(F検定) #p<0.05:非等分散####
library(tidyverse)
data1_GCGCs <- data1 %>% filter(Regimen == c("GC","GCs"))
data1_GCGCa <- data1 %>% filter(Regimen == c("GC","GCa"))
data1_GCaGCs <- data1 %>% filter(Regimen == c("GCa","GCs"))
var.test(eGFR~Regimen, data = data1_GCGCs)
var.test(eGFR~Regimen, data = data1_GCGCa)
var.test(eGFR~Regimen, data = data1_GCaGCs)
var.test(Ccr~Regimen, data = data1_GCGCs)
var.test(Ccr~Regimen, data = data1_GCGCa)
var.test(Ccr~Regimen, data = data1_GCaGCs)
var.test(CIN~Regimen, data = data1_GCGCs)
var.test(CIN~Regimen, data = data1_GCGCa)
var.test(CIN~Regimen, data = data1_GCaGCs)
var.test(CIV~Regimen, data = data1_GCGCs)
var.test(CIV~Regimen, data = data1_GCGCa)
var.test(CIV~Regimen, data = data1_GCaGCs)

data2_GCGCs <- data2 %>% filter(Regimen == c("GC","GCs"))
data2_GCGCa <- data2 %>% filter(Regimen == c("GC","GCa"))
data2_GCaGCs <- data2 %>% filter(Regimen == c("GCa","GCs"))
var.test(Age~Regimen, data = data2_GCGCs)
var.test(Age~Regimen, data = data2_GCGCa)
var.test(Age~Regimen, data = data2_GCaGCs)

data3_GCGCs <- data3 %>% filter(Regimen == c("GC","GCs"))
data3_GCGCa <- data3 %>% filter(Regimen == c("GC","GCa"))
data3_GCaGCs <- data3 %>% filter(Regimen == c("GCa","GCs"))
var.test(Time~Regimen, data = data3_GCGCs)
var.test(Time~Regimen, data = data3_GCGCa)
var.test(Time~Regimen, data = data3_GCaGCs)

#Bartlett test(多群間の等分散性の検定) #p<0.05:非等分散####
bartlett.test(formula=data1$eGFR~data1$Regimen)
bartlett.test(formula=data1$Ccr~data1$Regimen)
bartlett.test(formula=data1$CIN~data1$Regimen)
bartlett.test(formula=data1$CIV~data1$Regimen)
bartlett.test(formula=data2$Age~data2$Regimen)
bartlett.test(formula=data3$Time~data3$Regimen)

#Levene's test（多群間の等分散性の検定:非正規分布データでも頑健）#p<0.05:非等分散####
#install.packages("lawstat", repos="http://cran.ism.ac.jp/")
library(lawstat)
levene.test(y=data1$eGFR,group=data1$Regimen)
levene.test(y=data1$Ccr,group=data1$Regimen)
levene.test(y=data1$CIN,group=data1$Regimen)
levene.test(y=data1$CIV,group=data1$Regimen)
#levene.test(y=data2$Age,group=data2$Regimen)
data2_rm <- data2 %>% filter(!is.na(Age))
levene.test(y=data2_rm$Age,group=data2_rm$Regimen)
levene.test(y=data3$Time,group=data3$Regimen)

#Hartley test（多群間の等分散性の検定:検出力(1-β)に優れる)#p<0.05:非等分散####
#各水準におけるサンプルサイズが揃っている必要あり→適用不可
#install.packages("SuppDists", repos="http://cran.ism.ac.jp/")
#library(SuppDists)

#Kruskal-Wallis one-way analysis of variance #p<0.05: 3群の代表値に差がある####
kruskal.test(x=list(data1_GC$eGFR,data1_GCa$eGFR,data1_GCs$eGFR))
kruskal.test(x=list(data1_GC$Ccr,data1_GCa$Ccr,data1_GCs$Ccr))
kruskal.test(x=list(data1_GC$CIN,data1_GCa$CIN,data1_GCs$CIN))
kruskal.test(x=list(data1_GC$CIV,data1_GCa$CIV,data1_GCs$CIV))
kruskal.test(x=list(data2_GC$Age,data2_GCa$Age,data2_GCs$Age))
kruskal.test(x=list(data3_GC$Time,data3_GCa$Time,data3_GCs$Time))

KW_eGFR <- kruskal.test(data1$eGFR~data1$Regimen)
KW_Ccr <- kruskal.test(data1$Ccr~data1$Regimen)
KW_CIN <- kruskal.test(data1$CIN~data1$Regimen)
KW_CIV <- kruskal.test(data1$CIV~data1$Regimen)
KW_Age <- kruskal.test(data2$Age~data2$Regimen)
KW_Time <- kruskal.test(data3$Time~data3$Regimen)

#効果量イータ2乗
eta2_eGFR <- KW_eGFR$statistic/(length(data1$eGFR)-1)
eta2_Ccr <- KW_Ccr$statistic/(length(data1$Ccr)-1)
eta2_CIN <- KW_CIN$statistic/(length(data1$CIN)-1)
eta2_CIV <- KW_CIV$statistic/(length(data1$CIV)-1)
eta2_Age <- KW_Age$statistic/(length(data2$Age)-1)
eta2_Time <- KW_Time$statistic/(length(data3$Time)-1)
names(eta2_eGFR) <- "eta2"
names(eta2_Ccr) <- "eta2"
names(eta2_CIN) <- "eta2"
names(eta2_CIV) <- "eta2"
names(eta2_Age) <- "eta2"
names(eta2_Time) <- "eta2"

#Bonferroni's multiple comparison test####
(BF_eGFR <- pairwise.wilcox.test(data1$eGFR, data1$Regimen, p.adj="bonferroni", exact=F))
(BF_Ccr <- pairwise.wilcox.test(data1$Ccr, data1$Regimen, p.adj="bonferroni", exact=F))
(BF_CIN <- pairwise.wilcox.test(data1$CIN, data1$Regimen, p.adj="bonferroni", exact=F))
(BF_CIV <- pairwise.wilcox.test(data1$CIV, data1$Regimen, p.adj="bonferroni", exact=F))
(BF_Age <- pairwise.wilcox.test(data2$Age, data2$Regimen, p.adj="bonferroni", exact=F))
(BF_Time <- pairwise.wilcox.test(data3$Time, data3$Regimen, p.adj="bonferroni", exact=F))

#Holm####
(Holm_eGFR <- pairwise.wilcox.test(data1$eGFR, data1$Regimen, p.adj="holm", exact=F))
(Holm_Ccr <- pairwise.wilcox.test(data1$Ccr, data1$Regimen, p.adj="holm", exact=F))
(Holm_CIN <- pairwise.wilcox.test(data1$CIN, data1$Regimen, p.adj="holm", exact=F))
(Holm_CIV <- pairwise.wilcox.test(data1$CIV, data1$Regimen, p.adj="holm", exact=F))
(Holm_Age <- pairwise.wilcox.test(data2$Age, data2$Regimen, p.adj="holm", exact=F))
(Holm_Time <- pairwise.wilcox.test(data3$Time, data3$Regimen, p.adj="holm", exact=F))

#マン・ホイットニーnのU検定(Wilcoxonの順位和検定)####
C1.C2 <- subset(data1, Regimen=="GC" | Regimen=="GCa")
C1.C3 <- subset(data1, Regimen=="GC" | Regimen=="GCs")
C2.C3 <- subset(data1, Regimen=="GCa" | Regimen=="GCs")
WC_GCGCa_eGFR <- wilcox.test(C1.C2$eGFR~C1.C2$Regimen, correct=FALSE) 
WC_GCGCs_eGFR <- wilcox.test(C1.C3$eGFR~C1.C3$Regimen, correct=FALSE) 
WC_GCaGCs_eGFR <- wilcox.test(C2.C3$eGFR~C2.C3$Regimen, correct=FALSE)
WC_GCGCa_Ccr <- wilcox.test(C1.C2$Ccr~C1.C2$Regimen, correct=FALSE) 
WC_GCGCs_Ccr <- wilcox.test(C1.C3$Ccr~C1.C3$Regimen, correct=FALSE) 
WC_GCaGCs_Ccr <- wilcox.test(C2.C3$Ccr~C2.C3$Regimen, correct=FALSE)
WC_GCGCa_CIN <- wilcox.test(C1.C2$CIN~C1.C2$Regimen, correct=FALSE) 
WC_GCGCs_CIN <- wilcox.test(C1.C3$CIN~C1.C3$Regimen, correct=FALSE) 
WC_GCaGCs_CIN <- wilcox.test(C2.C3$CIN~C2.C3$Regimen, correct=FALSE)
WC_GCGCa_CIV <- wilcox.test(C1.C2$CIV~C1.C2$Regimen, correct=FALSE) 
WC_GCGCs_CIV <- wilcox.test(C1.C3$CIV~C1.C3$Regimen, correct=FALSE) 
WC_GCaGCs_CIV <- wilcox.test(C2.C3$CIV~C2.C3$Regimen, correct=FALSE)

C1.C2 <- subset(data2, Regimen=="GC" | Regimen=="GCa")
C1.C3 <- subset(data2, Regimen=="GC" | Regimen=="GCs")
C2.C3 <- subset(data2, Regimen=="GCa" | Regimen=="GCs")
WC_GCGCa_Age <- wilcox.test(C1.C2$Age~C1.C2$Regimen, correct=FALSE) 
WC_GCGCs_Age <- wilcox.test(C1.C3$Age~C1.C3$Regimen, correct=FALSE) 
WC_GCaGCs_Age <- wilcox.test(C2.C3$Age~C2.C3$Regimen, correct=FALSE)

C1.C2 <- subset(data3, Regimen=="GC" | Regimen=="GCa")
C1.C3 <- subset(data3, Regimen=="GC" | Regimen=="GCs")
C2.C3 <- subset(data3, Regimen=="GCa" | Regimen=="GCs")
WC_GCGCa_Time <- wilcox.test(C1.C2$Time~C1.C2$Regimen, correct=FALSE) 
WC_GCGCs_Time <- wilcox.test(C1.C3$Time~C1.C3$Regimen, correct=FALSE) 
WC_GCaGCs_Time <- wilcox.test(C2.C3$Time~C2.C3$Regimen, correct=FALSE)

#Steel-Dwass（スティールドゥワス）検定####
Steel.Dwass <- function(data,group){
  OK <- complete.cases(data, group)
  data <- data[OK]
  group <- group[OK]
  n.i <- table(group)
  ng <- length(n.i)
  t <- combn(ng, 2, function(ij) {
    i <- ij[1]
    j <- ij[2]
    r <- rank(c(data[group == i], data[group == j]))
    R <- sum(r[1:n.i[i]])
    N <- n.i[i]+n.i[j]
    E <- n.i[i]*(N+1)/2
    V <- n.i[i]*n.i[j]/(N*(N-1))*(sum(r^2)-N*(N+1)^2/4)
    return(abs(R-E)/sqrt(V))
  })
  p <- ptukey(t*sqrt(2), ng, Inf, lower.tail=FALSE)
  result <- cbind(t, p)
  rownames(result) <- combn(ng, 2, paste, collapse=":")
  return(result)
}

data1 <- mutate(data1, RegimenNo = recode(Regimen, 
                                 `GC` = 1L, `GCa` = 2L, `GCs` = 3L, 
                                 .missing = -1L))
Steel.Dwass(data1$eGFR, as.numeric(data1$RegimenNo))
Steel.Dwass(data1$Ccr, as.numeric(data1$RegimenNo))
Steel.Dwass(data1$CIN, as.numeric(data1$RegimenNo))
Steel.Dwass(data1$CIV, as.numeric(data1$RegimenNo))

data2 <- mutate(data2, RegimenNo = recode(Regimen, 
                                          `GC` = 1L, `GCa` = 2L, `GCs` = 3L, 
                                          .missing = -1L))
Steel.Dwass(data2$Age, as.numeric(data2$RegimenNo))

data3 <- mutate(data3, RegimenNo = recode(Regimen, 
                                          `GC` = 1L, `GCa` = 2L, `GCs` = 3L, 
                                          .missing = -1L))
Steel.Dwass(data3$Time, as.numeric(data3$RegimenNo))

#Fisher's exact test####
library(tidyverse)
#Patient characteristics####
GCGCsGCa <- data4 %>% filter(Regimen == c("GC","GCs","GCa")) %>% select(c("Male","Female"))
GCGCs <- data4 %>% filter(Regimen == c("GC","GCs")) %>% select(c("Male","Female"))
GCGCa <- data4 %>% filter(Regimen == c("GC","GCa")) %>% select(c("Male","Female"))
GCsGCa <- data4 %>% filter(Regimen == c("GCs","GCa")) %>% select(c("Male","Female"))
#GCGCsGCa <- data4[,grep("^GC$|^GCs$|^GCa$",colnames(data4))]
#GCGCs <- data4[,grep("^GC$|^GCs$",colnames(data4))]
#GCGCa <- data4[,grep("^GC$|^GCa$",colnames(data4))]
#GCsGCa <- data4[,grep("^GCs$|^GCa$",colnames(data4))]
fisher.test(GCGCsGCa)
fisher.test(GCGCs)
fisher.test(GCGCa)
fisher.test(GCsGCa)
#Serotonin 5-HT3 receptor antagonist####
#ALL
GCGCsGCa5 <- data5 %>% filter(Regimen == c("GC","GCs","GCa")) %>% select(-Regimen)
GCGCs5 <- data5 %>% filter(Regimen != "GCa") %>% select(-Regimen)
GCGCa5 <- data5 %>% filter(Regimen != "GCs") %>% select(-Regimen)
GCsGCa5 <- data5 %>% filter(Regimen != "GC") %>% select(-Regimen)
fisher.test(GCGCsGCa5)
fisher.test(GCGCs5)
fisher.test(GCGCa5)
fisher.test(GCsGCa5)
#1st_vs_2nd
GCGCsGCa5_2 <- GCGCsGCa5 %>% mutate(FGAE = Ramosetron + Granisetron) %>% mutate(SGAE = Palonosetron) %>% select(FGAE, SGAE)
GCGCs5_2 <- GCGCs5 %>% mutate(FGAE = Ramosetron + Granisetron) %>% mutate(SGAE = Palonosetron) %>% select(FGAE, SGAE)
GCGCa5_2 <- GCGCa5 %>% mutate(FGAE = Ramosetron + Granisetron) %>% mutate(SGAE = Palonosetron) %>% select(FGAE, SGAE)
GCsGCa5_2 <- GCsGCa5 %>% mutate(FGAE = Ramosetron + Granisetron) %>% mutate(SGAE = Palonosetron) %>% select(FGAE, SGAE)
fisher.test(GCGCsGCa5_2)
fisher.test(GCGCs5_2)
fisher.test(GCGCa5_2)
fisher.test(GCsGCa5_2)
#Additional antiemetics####
GCGCsGCa6 <- data6 %>% filter(Regimen == c("GC","GCs","GCa")) %>% select(-Regimen)
GCGCs6 <- data6 %>% filter(Regimen != "GCa") %>% select(-Regimen)
GCGCa6 <- data6 %>% filter(Regimen != "GCs") %>% select(-Regimen)
GCsGCa6 <- data6 %>% filter(Regimen != "GC") %>% select(-Regimen)
fisher.test(GCGCsGCa6)
fisher.test(GCGCs6)
fisher.test(GCGCa6)
fisher.test(GCsGCa6)
#Treatment response####
GCGCsGCa7 <- data7 %>% filter(Regimen == c("GC","GCs","GCa")) %>% select(-Regimen)
GCGCs7 <- data7 %>% filter(Regimen != "GCa") %>% select(-Regimen)
GCGCa7 <- data7 %>% filter(Regimen != "GCs") %>% select(-Regimen)
GCsGCa7 <- data7 %>% filter(Regimen != "GC") %>% select(-Regimen)
fisher.test(GCGCsGCa7)
fisher.test(GCGCs7)
fisher.test(GCGCa7)
fisher.test(GCsGCa7)
#Prevalence of CINV on each treatment day####
df <- data8
for(i in 9:10){
  #df2 <- df %>% filter(Day == i) %>% select(CIN, noCIN)
  #print(fisher.test(df2))
  df2 <- df %>% filter(Regimen != "GCa") %>% filter(Day == i) %>% select(CIN, noCIN)
  print(fisher.test(df2))
  df2 <- df %>% filter(Regimen != "GCs") %>% filter(Day == i) %>% select(CIN, noCIN)
  print(fisher.test(df2))
  df2 <- df %>% filter(Regimen != "GC") %>% filter(Day == i) %>% select(CIN, noCIN)
  print(fisher.test(df2))
  #df2 <- df %>% filter(Day == i) %>% select(CIV, noCIV)
  #print(fisher.test(df2))
  df2 <- df %>% filter(Regimen != "GCa") %>% filter(Day == i) %>% select(CIV, noCIV)
  print(fisher.test(df2))
  df2 <- df %>% filter(Regimen != "GCs") %>% filter(Day == i) %>% select(CIV, noCIV)
  print(fisher.test(df2))
  df2 <- df %>% filter(Regimen != "GC") %>% filter(Day == i) %>% select(CIV, noCIV)
  print(fisher.test(df2))
}
#############################################################
#カプランマイヤー曲線####
#install.packages("survminer")
library(MASS)
library(survival)
library(survminer)
#1 gehan データを load
#data(gehan)
data3 <- read_excel("Book2.xlsx", 5)
#2 生存時間分析のための Surv オブジェクトを survival::Surv で作成
#surv.obj <- Surv(time = gehan$time, event = gehan$cens)
surv.obj <- Surv(time = data3$Time, event = data3$Cens)
#3 作成した Surv オブジェクトを引数とし，survival::survfit で投薬（treat）の有無毎に KM 曲線を作成
#ge.sf <- survfit(surv.obj ~ treat, data = gehan)
ge.sf <- survfit(surv.obj ~ Regimen, data = data3)
#4 作成した KM 曲線の結果の概要を survminer::surv_summary でデータフレームとして作成
#ge.sf.df <- surv_summary(ge.sf, data = gehan)
ge.sf.df <- surv_summary(ge.sf, data = data3)
View(ge.sf.df)
#5
plot <- ggsurvplot(
  fit = ge.sf,
  data = data3,
  conf.int = F,
  pval = T,
  risk.table = F,
  cumevents = F,
  cumcensor = F,
  ggtheme = ggplot2::theme_light(),
  tables.height = 0.15
)
svg(file="plot.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
#############################################################
#violin plot
library(tidyverse)
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)
#install.packages("devtools")
library(devtools)
devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")
library(reshape2)
dat <- melt(data1) #Using Regimen as id variables
dat_CIN <- dat %>% filter(variable == "CIN")
dat_CIV <- dat %>% filter(variable == "CIV")

vio_dot_plot <- function(x){
  dat <- x
  ggplot(data = dat, mapping = aes(x = Regimen, y = value, fill = Regimen)) +
    geom_violin(scale = "count", trim = F) +
    stat_summary(fun.data = mean_sdl, 
                 fun.args = list(mult = 1), 
                 geom = "pointrange",
                 position = position_nudge(0.05),
                 color = c("dark red", "dark green", "dark blue")) +
    geom_dotplot(binaxis = "y", dotsize = 0.3, stackdir = "down", binwidth = 0.05, 
                 position = position_nudge(-0.025),
                 #position = position_jitterdodge(jitter.width=0, dodge.width=0),
    ) + 
    facet_wrap(~variable, ncol = 5, scales = "free") +
    theme(legend.position = "none") + #theme_bw()
    labs(x = "Regimen", y = "CINV (days)")
}
plot <- vio_dot_plot(dat_CIN)
svg(file="plot_CIN.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
plot <- vio_dot_plot(dat_CIV)
svg(file="plot_CIV.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
#############################################################
#violin plot
vio_median_plot <- function(x){
  dat <- x
  ggplot(data = dat, mapping = aes(x = Regimen, y = value, fill = Regimen)) +
    geom_violin(scale = "count", trim = F) +
    stat_summary(fun.data = median_hilow,
                 geom = "pointrange",
                 size = 0.7,
                 linetype = 1,
                 position = position_nudge(0)) + #,
                 #color = c("dark red", "dark green", "dark blue")) + 
    #geom_poiontrange(data = dat, width=0.2, size=1) +
    facet_wrap(~variable, ncol = 5, scales = "free") +
    theme(legend.position = "none") + #theme_bw()
    labs(x = "Regimen", y = "CINV (days)")
}
vio_median_plot(dat_CIN)
vio_median_plot(dat_CIV)
plot <- vio_median_plot(dat_CIN)
svg(file="plot_median_CIN.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
plot <- vio_median_plot(dat_CIV)
svg(file="plot_median_CIV.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
#############################################################
#violin_box plot
vio_box_plot <- function(x){
  dat <- x
  ggplot(data = dat, mapping = aes(x = Regimen, y = value, fill = Regimen)) +
    geom_violin(scale = "count", trim = F, width = 2.5, fill = "dark gray") +
    geom_boxplot(width = 0.05, fill = "white", alpha = 1) +
    stat_boxplot(geom = "errorbar", width = 0.1) +
    #stat_summary(fun.data = median_hilow, geom = "pointrange", size = 0.7, linetype = 2, position = position_nudge(0)) +
    #geom_count() +
    #geom_poiontrange(data = dat, width=0.2, size=1) +
    ylim(c(0, NA)) +
    facet_wrap(~variable, ncol = 5, scales = "free") +
    theme(legend.position = "none") + #theme_bw()
    labs(x = "Regimen", y = "CINV (days)")
}
vio_box_plot(dat_CIN)
vio_box_plot(dat_CIV)
plot <- vio_box_plot(dat_CIN)
svg(file="vio_box_plot_CIN.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
plot <- vio_box_plot(dat_CIV)
svg(file="vio_box_plot_CIV.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
#############################################################
#box plot
box_plot <- function(x){
  dat <- x
  ggplot(data = dat, mapping = aes(x = Regimen, y = value, fill = Regimen)) +
    geom_boxplot(width = 0.05, fill = "gray", alpha = 1) +
    #facet_wrap(~variable, ncol = 5, scales = "free") +
    theme(legend.position = "none") + #theme_bw()
    labs(x = "Regimen", y = "CINV (days)")
}
box_plot(dat_CIN)
box_plot(dat_CIV)
plot <- box_plot(dat_CIN)
svg(file="boxplot_CIN.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
plot <- box_plot(dat_CIV)
svg(file="boxplot_CIV.svg") #ファイル名指定
print(plot) #プロット作成
dev.off() #svg出力
#############################################################
#############################################################
#他のプロジェクトの解析スクリプト
#############################################################
#############################################################
#二標本t検定(等分散を仮定)
t.test(scz_f$Age, con_f$Age, var.equal=T) #女性
t.test(Age ~ Diagnosis1.Control, data = female, var.equal=T) #女性
t.test(scz_m$Age, con_m$Age, var.equal=T) #男性
t.test(Age ~ Diagnosis1.Control, data = male, var.equal=T) #男性
t.test(Education ~ Diagnosis1, data = data, var.equal=T) #教育p<0.001

#二標本t検定(ウェルチの検定)
t.test(Age~Diagnosis1, data = data) #Welch's t test(等分散を問わない)
t.test(scz_f$Age, con_f$Age, var.equal=F) #女性
t.test(Age ~ Diagnosis1.Control, data = female, var.equal=F) #女性
t.test(scz_m$Age, con_m$Age, var.equal=F) #男性
t.test(Age ~ Diagnosis1.Control, data = male, var.equal=F) #男性
t.test(scz$Education, con$Education, var.equal=F) #教育p<0.001
t.test(Education ~ Diagnosis1, data = data, var.equal=F) #教育p<0.001

#ウィルコクソン順位和検定・マン・ホイットニーU検定
#install.packages("coin")
#install.packages("tidyverse")
#パッケージの読み込み
library(coin)
library(tidyverse)
wilcox.test(scz_f$Age, con_f$Age,) #女性
wilcox.test(scz_m$Age, con_m$Age,) #男性
wilcox_test(Age ~ Diagnosis1.Control, data = female) #女性
wilcox_test(Age ~ Diagnosis1.Control, distribution="exact", data = female) #女性
#変数のクラス確認
sapply(female, class)
#numericをfactorに変換
female[,2] <- as.factor(female[,2])
class(data[,2])
wilcox_test(Education ~ Diagnosis1, distribution="exact", data = data) #教育p<0.001


#############################################################
#回帰分析
#Excelファイル読み込み
#install.packages("readxl")
#y
#library(readxl)
#data1 <- read_excel("GxInfo2.xlsx", 1) #シート1の読み込み
#str(data1)
#ダミー変数加工
#install.packages("caret")
#library(caret)
#library(ggplot2)
#tmp <- dummyVars(~., data = data1) #全質的変数を対象
#data1.dummy <- as.data.frame(predict(tmp, data1)) #tmpの質的変数をダミー変数に変換
#str(data1.dummy) #Df構造表示
#head(data1.dummy)
#############################################################
#CSVファイル読み込み
data2 <- read.csv("GxInfo2.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(data2)
#ダミー変数加工
#install.packages("caret")
library(caret)
library(ggplot2)
tmp2 <- dummyVars(~.-ID, data = data2) #全質的変数を対象
data2.dummy <- as.data.frame(predict(tmp2, data2)) #tmpの質的変数をダミー変数に変換
str(data2.dummy) #Df構造表示
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(data2.dummy, sheetname = "sheet1", file = "out_file.xlsx") 
#SCZ,CON#####################################################
#CSVファイル読み込み
data <- read.csv("GxInfo3.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(data)
#ダミー変数加工
#install.packages("caret")
library(caret)
library(ggplot2)
tmp <- dummyVars(~.-ID, data = data) #全質的変数を対象(ID除外)
data.dummy <- as.data.frame(predict(tmp, data)) #tmpの質的変数をダミー変数に変換
str(data.dummy) #Df構造表示
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(data.dummy, sheetname = "sheet1", file = "out_file.xlsx") 
#############################################################
#散布図行列
#install.packages("car")
library("car")
scatterplotMatrix(data)
#############################################################
#単回帰分析
result1 <- lm(ddCt~data.dummy[,1],data = data.dummy)
result2 <- lm(ddCt~data.dummy[,2],data = data.dummy)
result3 <- lm(ddCt~data.dummy[,3],data = data.dummy)
result4 <- lm(ddCt~data.dummy[,4],data = data.dummy)
result5 <- lm(ddCt~data.dummy[,5],data = data.dummy)
result6 <- lm(ddCt~data.dummy[,6],data = data.dummy)
result7 <- lm(ddCt~data.dummy[,7],data = data.dummy)
result8 <- lm(ddCt~data.dummy[,8],data = data.dummy)
result9 <- lm(ddCt~data.dummy[,9],data = data.dummy)
result10 <- lm(ddCt~data.dummy[,10],data = data.dummy)
result11 <- lm(ddCt~data.dummy[,11],data = data.dummy)
result12 <- lm(ddCt~data.dummy[,12],data = data.dummy)
result13 <- lm(ddCt~data.dummy[,13],data = data.dummy)
result14 <- lm(ddCt~data.dummy[,14],data = data.dummy)
result15 <- lm(ddCt~data.dummy[,15],data = data.dummy)
result16 <- lm(ddCt~data.dummy[,16],data = data.dummy)
result17 <- lm(ddCt~data.dummy[,17],data = data.dummy)
result18 <- lm(ddCt~data.dummy[,18],data = data.dummy)
result19 <- lm(ddCt~data.dummy[,19],data = data.dummy)
result20 <- lm(ddCt~data.dummy[,20],data = data.dummy)
result21 <- lm(ddCt~data.dummy[,21],data = data.dummy)
result22 <- lm(ddCt~data.dummy[,22],data = data.dummy)
result23 <- lm(ddCt~data.dummy[,23],data = data.dummy)
result24 <- lm(ddCt~data.dummy[,24],data = data.dummy)
result25 <- lm(ddCt~data.dummy[,25],data = data.dummy)
summary(result1)
summary(result2)
summary(result3)
summary(result4)
summary(result5)
summary(result6)
summary(result7)
summary(result8)
summary(result9)
summary(result10)
summary(result11)
summary(result12)
summary(result13)
summary(result14)
summary(result15)
summary(result16)
summary(result17)
summary(result18)
summary(result19)
summary(result20)
summary(result21)
summary(result22)
summary(result23)
summary(result24)
summary(result25)
#############################################################
#重回帰分析
result <- lm(ddCt~.-No-ddCt,data = data.dummy)
summary(result)
#result <- lm(ddCt~Label.CNVdup+Label.CON+Label.SCZ,data = data.dummy)
#summary(result)
#result <- lm(ddCt~Diagnosis1.Schizophrenia+Label.CNVdup+Label.SCZ
#             +Gender.M+Age+Education+DominantHand.Ambidextrous
#             +DominantHand.Righthanded,data = data.dummy)
#summary(result)
#欠損値除外(R step na.omit)
data.dummy[complete.cases(data.dummy), ]
#ステップワイズ回帰(後退)
full.model <- lm(ddCt ~.,data = data.dummy)
reduced.model <- step(full.model, direction = "backward")
#ステップワイズ回帰（前進）
min.model <- lm(ddCt ~ 1,data = data.dummy)
fwd.model <- step(min.model, direction = "forward", 
                  scope = (~Diagnosis1.Schizophrenia+Label.CNVdup+Label.SCZ
                           +Gender.M+Age+Education+DominantHand.Ambidextrous
                           +DominantHand.Righthanded))
summary(result)
#SCZ#########################################################
#CSVファイル読み込み
data <- read.csv("GxInfo4.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(data)
#ダミー変数加工
#install.packages("caret")
library(caret)
library(ggplot2)
tmp <- dummyVars(~.-ID, data = data) #全質的変数を対象(ID除外)
data.dummy <- as.data.frame(predict(tmp, data)) #tmpの質的変数をダミー変数に変換
str(data.dummy) #Df構造表示
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(data.dummy, sheetname = "sheet1", file = "out_file.xlsx") 
#############################################################
#散布図行列
#install.packages("car")
library("car")
scatterplotMatrix(data)
#############################################################
#単回帰分析
result1 <- lm(ddCt~data.dummy[,1],data = data.dummy)
result2 <- lm(ddCt~data.dummy[,2],data = data.dummy)
result3 <- lm(ddCt~data.dummy[,3],data = data.dummy)
result4 <- lm(ddCt~data.dummy[,4],data = data.dummy)
result5 <- lm(ddCt~data.dummy[,5],data = data.dummy)
result6 <- lm(ddCt~data.dummy[,6],data = data.dummy)
result7 <- lm(ddCt~data.dummy[,7],data = data.dummy)
result8 <- lm(ddCt~data.dummy[,8],data = data.dummy)
result9 <- lm(ddCt~data.dummy[,9],data = data.dummy)
result10 <- lm(ddCt~data.dummy[,10],data = data.dummy)
result11 <- lm(ddCt~data.dummy[,11],data = data.dummy)
result12 <- lm(ddCt~data.dummy[,12],data = data.dummy)
result13 <- lm(ddCt~data.dummy[,13],data = data.dummy)
result14 <- lm(ddCt~data.dummy[,14],data = data.dummy)
result15 <- lm(ddCt~data.dummy[,15],data = data.dummy)
result16 <- lm(ddCt~data.dummy[,16],data = data.dummy)
result17 <- lm(ddCt~data.dummy[,17],data = data.dummy)
result18 <- lm(ddCt~data.dummy[,18],data = data.dummy)
result19 <- lm(ddCt~data.dummy[,19],data = data.dummy)
result20 <- lm(ddCt~data.dummy[,20],data = data.dummy)
result21 <- lm(ddCt~data.dummy[,21],data = data.dummy)
result22 <- lm(ddCt~data.dummy[,22],data = data.dummy)
result23 <- lm(ddCt~data.dummy[,23],data = data.dummy)
result24 <- lm(ddCt~data.dummy[,24],data = data.dummy)
result25 <- lm(ddCt~data.dummy[,25],data = data.dummy)
result26 <- lm(ddCt~data.dummy[,26],data = data.dummy)
result27 <- lm(ddCt~data.dummy[,27],data = data.dummy)
result28 <- lm(ddCt~data.dummy[,28],data = data.dummy)
result29 <- lm(ddCt~data.dummy[,29],data = data.dummy)
result30 <- lm(ddCt~data.dummy[,30],data = data.dummy)
result31 <- lm(ddCt~data.dummy[,31],data = data.dummy)
result32 <- lm(ddCt~data.dummy[,32],data = data.dummy)
result33 <- lm(ddCt~data.dummy[,33],data = data.dummy)
result34 <- lm(ddCt~data.dummy[,34],data = data.dummy)
result35 <- lm(ddCt~data.dummy[,35],data = data.dummy)
result36 <- lm(ddCt~data.dummy[,36],data = data.dummy)
result37 <- lm(ddCt~data.dummy[,37],data = data.dummy)
result38 <- lm(ddCt~data.dummy[,38],data = data.dummy)
result39 <- lm(ddCt~data.dummy[,39],data = data.dummy)
result40 <- lm(ddCt~data.dummy[,40],data = data.dummy)
result41 <- lm(ddCt~data.dummy[,41],data = data.dummy)
result42 <- lm(ddCt~data.dummy[,42],data = data.dummy)
result43 <- lm(ddCt~data.dummy[,43],data = data.dummy)
result44 <- lm(ddCt~data.dummy[,44],data = data.dummy)
result45 <- lm(ddCt~data.dummy[,45],data = data.dummy)
result46 <- lm(ddCt~data.dummy[,46],data = data.dummy)
result47 <- lm(ddCt~data.dummy[,47],data = data.dummy)
result48 <- lm(ddCt~data.dummy[,48],data = data.dummy)
result49 <- lm(ddCt~data.dummy[,49],data = data.dummy)
result50 <- lm(ddCt~data.dummy[,50],data = data.dummy)
result51 <- lm(ddCt~data.dummy[,51],data = data.dummy)
result52 <- lm(ddCt~data.dummy[,52],data = data.dummy)
result53 <- lm(ddCt~data.dummy[,53],data = data.dummy)
result54 <- lm(ddCt~data.dummy[,54],data = data.dummy)
result55 <- lm(ddCt~data.dummy[,55],data = data.dummy)
result56 <- lm(ddCt~data.dummy[,56],data = data.dummy)
result57 <- lm(ddCt~data.dummy[,57],data = data.dummy)
result58 <- lm(ddCt~data.dummy[,58],data = data.dummy)
summary(result1)
summary(result2)
summary(result3)
summary(result4)
summary(result5)
summary(result6)
summary(result7)
summary(result8)
summary(result9)
summary(result10)
summary(result11)
summary(result12)
summary(result13)
summary(result14)
summary(result15)
summary(result16)
summary(result17)
summary(result18)
summary(result19)
summary(result20)
summary(result21)
summary(result22)
summary(result23)
summary(result24)
summary(result25)
summary(result26)
summary(result27)
summary(result28)
summary(result29)
summary(result30)
summary(result31)
summary(result32)
summary(result33)
summary(result34)
summary(result35)
summary(result36)
summary(result37)
summary(result38)
summary(result39)
summary(result40)
summary(result41)
summary(result42)
summary(result43)
summary(result44)
summary(result45)
summary(result46)
summary(result47)
summary(result48)
summary(result49)
summary(result50)
summary(result51)
summary(result52)
summary(result53)
summary(result54)
summary(result55)
summary(result56)
summary(result57)
summary(result58)
#############################################################
#重回帰分析
result <- lm(ddCt~.-No-ddCt,data = data.dummy)
summary(result)
#result <- lm(ddCt~Label.CNVdup+Label.SCZ+Gender.F+Gender.M+Age+OnSet+DUrationOfIllness+Education
#              +Episode.Multiple_episodes+DominantHand.Ambidextrous+DominantHand.Lefthanded+DominantHand.Righthanded
#              +Type.Hebephrenic+Type.Paranoid+Type.Residual+Type.Undifferentiated
#              +CPZ+CPZ_FGA+CPZ_SGA
#              +AP1.FGA+AP1.FGA_SGA+AP1.SGA
#              +AP2.APZ+AP2.BPD+AP2.BPD_LPZ_OLZ+AP2.BPD_RIS+AP2.CPZ_OLZ_RIS+AP2.CPZ_QTP_RIS+AP2.HPD_OLZ+AP2.LPZ_RIS+AP2.OLZ+AP2.PAL+AP2.PER_RIS+AP2.QTP+AP2.QTP_RIS+AP2.RIS
#              +FGA1.BPD+FGA1.BPD_LPZ+FGA1.CPZ+FGA1.HPD+FGA1.LPZ
#              +SGA1.APZ+SGA1.OLZ+SGA1.OLZ_RIS+SGA1.PAL+SGA1.PER_RIS+SGA1.QTP+SGA1.QTP_RIS+SGA1.RIS,data = data4.dummy)
#summary(result)
result4 <- lm(ddCt~Label.CNVdup+Gender.M+Age+OnSet+DUrationOfIllness+Education
              +Episode.Multiple_episodes+DominantHand.Ambidextrous+DominantHand.Lefthanded+DominantHand.Righthanded
              +Type.Hebephrenic+Type.Paranoid+Type.Residual
              +CPZ
              +FGA1.
              +SGA1.,data = data4.dummy)
summary(result4)
#欠損値除外(R step na.omit)
data4.dummy[complete.cases(data4.dummy), ]
#ステップワイズ回帰(後退)
full.model <- lm(ddCt ~.,data = data4.dummy)
reduced.model <- step(full.model, direction = "backward")
#ステップワイズ回帰（前進）
min.model <- lm(ddCt ~ 1,data = data4.dummy)
fwd.model <- step(min.model, direction = "forward", scope = (~Label.CNVdup+Label.SCZ+Gender.F+Gender.M+Age+OnSet+DUrationOfIllness+Education
                                                             +Episode.Multiple_episodes+DominantHand.Ambidextrous+DominantHand.Lefthanded+DominantHand.Righthanded
                                                             +Type.Hebephrenic+Type.Paranoid+Type.Residual+Type.Undifferentiated
                                                             +CPZ+CPZ_FGA+CPZ_SGA
                                                             +AP1.FGA+AP1.FGA_SGA+AP1.SGA
                                                             +AP2.APZ+AP2.BPD+AP2.BPD_LPZ_OLZ+AP2.BPD_RIS+AP2.CPZ_OLZ_RIS+AP2.CPZ_QTP_RIS+AP2.HPD_OLZ+AP2.LPZ_RIS+AP2.OLZ+AP2.PAL+AP2.PER_RIS+AP2.QTP+AP2.QTP_RIS+AP2.RIS
                                                             +FGA1.BPD+FGA1.BPD_LPZ+FGA1.CPZ+FGA1.HPD+FGA1.LPZ
                                                             +SGA1.APZ+SGA1.OLZ+SGA1.OLZ_RIS+SGA1.PAL+SGA1.PER_RIS+SGA1.QTP+SGA1.QTP_RIS+SGA1.RIS
))
summary(result4)
min.model <- lm(ddCt ~ 1,data = data4.dummy)
fwd.model <- step(min.model, direction = "forward", scope = (~Label.CNVdup+Gender.M+Age+OnSet+DUrationOfIllness+Education
                                                             +Episode.Multiple_episodes+DominantHand.Ambidextrous+DominantHand.Lefthanded+DominantHand.Righthanded
                                                             +Type.Hebephrenic+Type.Paranoid+Type.Residual
                                                             +CPZ
                                                             +FGA1.
                                                             +SGA1.
))
summary(result4)
#############################################################
#重回帰分析
result2 <- lm(ddCt ~. ,data = data2)
result2 <- lm(ddCt ~. ,data = data2.dummy)
summary(result2)
#欠損値除外(R step na.omit)
data2.dummy[complete.cases(data2.dummy), ]
#ステップワイズ回帰
full.model <- lm(ddCt ~.,data = data2.dummy)
reduced.model <- step(full.model, direction = "backward")
#############################################################
#plot(data.dummy)
#abline(result15, lwd = 1, col ="red")

