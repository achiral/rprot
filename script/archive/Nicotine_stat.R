#Nicotine_stat
#基本統計量
#https://note.com/kotoko_tyo/n/ne09e4398093d
#カテゴリカルデータ集計
#https://www1.doshisha.ac.jp/~mjin/R/Chap_45/45.html
#http://www.eeso.ges.kyoto-u.ac.jp/emm/wp-content/uploads/2010/11/questionary01.pdf
#http://www.housecat442.com/?p=346
#シャピロウィルク検定(Shapiro-Wilk test) 
#正規分布に従うものか否かを調べる検定法
#https://data-science.gr.jp/implementation/ist_r_shapiro_wilk_test.html
#回帰分析
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/Nicotine")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#############################################################
#Excelファイル読み込み
#install.packages("readxl")
#y
library(readxl)
data1 <- read_excel("Info.xlsx", 1) #シート1の読み込み
data2 <- read_excel("Gx.xlsx", sheet = 1) #シート1の読み込み
#dataシートのマージ
#data3 <- merge(data1, data2, by="ID",incomparables=NA) #識別子IDとしてdata1とdata2に共通する情報を統合,NA対象外(full_joinと同じ)
#dplyrパッケージによる共通列マージ(RCookbook2p171)
library(dplyr)
data4 <- left_join(data1, data2, by="ID")#識別子IDとしてdata1にdata2の情報を統合
view(data4)
#Excelに保存(RCookbook2p103)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(data4, sheetname = "sheet1", file = "out_file.xlsx") 
#SCZ,CON#####################################################
#Excelファイル読み込み
#install.packages("readxl")
#y
#library(readxl)
#data <- read_excel("GxInfo2.xlsx", 1) #シート1の読み込み
#str(data1)
#CSVファイル読み込み
data <- read.csv("GxInfo2.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(data)
#############################################################
#変数のクラス確認
sapply(data, class)
#characterをfactorに変換
#data[,1] <- as.factor(data[,1])
#class(data[,1])
#基本統計量の一覧表
library(psych)
out <- describeBy(data,group = data$Diagnosis1)
print(out)

#カテゴリカルデータの集計
tab <- xtabs(~Diagnosis1+Gender+Episode+DominantHand+Type+AP1+AP2+FGA1+SGA1, data = data)

tab1 <- xtabs(~Diagnosis1+Gender, data = data)
tab2 <- xtabs(~Diagnosis1+Episode, data = data)
tab3 <- xtabs(~Diagnosis1+DominantHand, data = data)
tab4 <- xtabs(~Diagnosis1+Type, data = data)
tab5 <- xtabs(~Diagnosis1+AP1, data = data)
tab6 <- xtabs(~Diagnosis1+AP2, data = data)
tab7 <- xtabs(~Diagnosis1+FGA1, data = data)
tab8 <- xtabs(~Diagnosis1+SGA1, data = data)
print(tab1)
print(tab2)
print(tab3)
print(tab4)
print(tab5)
print(tab6)
print(tab7)
print(tab8)

#filter機能でデータを絞る
library(dplyr)
cpz300 <- data %>% filter(CPZ < 300)
cpz300_600 <- data %>% filter(CPZ >= 300, CPZ < 600)
cpz600_1000 <- data %>% filter(CPZ >= 600, CPZ < 1000)
cpz1000 <- data %>% filter(CPZ >= 1000)



#ダミー変数加工
#install.packages("caret")
library(caret)
library(ggplot2)
tmp <- dummyVars(~.-ID, data = data) #全質的変数を対象
data.dummy <- as.data.frame(predict(tmp, data)) #tmpの質的変数をダミー変数に変換
str(data.dummy) #Df構造表示

#filter機能で統合失調症と健常者、男性と女性を絞る
scz <- data.dummy %>% filter(Diagnosis1.Schizophrenia == 1)
scz_f <- scz %>% filter(Gender.F == 1)
scz_m <- scz %>% filter(Gender.M == 1)
con <- data.dummy %>% filter(Diagnosis1.Control == 1)
con_f <- con %>% filter(Gender.F == 1)
con_m <- con %>% filter(Gender.M == 1)
female <- data.dummy %>% filter(Gender.F == 1)
male <- data.dummy %>% filter(Gender.M == 1)

#正規性の検定
shapiro.test(x=data$Age) 
shapiro.test(x=scz$Age) 
shapiro.test(x=scz_f$Age) #非正規性→ノンパラ
shapiro.test(x=scz_m$Age) 
shapiro.test(x=con$Age) 
shapiro.test(x=con_f$Age) 
shapiro.test(x=con_m$Age) 
shapiro.test(x=scz[,12]) #非正規性→ノンパラ
shapiro.test(x=con[,12]) #非正規性→ノンパラ
#変数のクラス確認
#sapply(con, class)
#numericをfactorに変換
#con[,12] <- as.numeric(con[,12])
#class(con[,12])

#等分散性の検定
var.test(Age~Diagnosis1, data = data)
var.test(Age~Diagnosis1.Control, data = female)
var.test(Age~Diagnosis1.Control, data = male)
var.test(Education~Diagnosis1, data = data) #p<0.05のため非等分散

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

