#RCookbook2
#Rグラフィックスクックブック第2版
#############################################################
#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/RGraphicsCookbook")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#############################################################
#P1
#パッケージインストール
#install.packages(c("tidyverse","gcookbook"))  #install.packages("tidyverse")と#install.packages("gcookbook")を一度に実行
#install.packages("ggplot2")
#install.packages("tidyverse")で読み込まれる

#ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(tidyverse)
#library(ggplot2) #library(tidyverse)で読み込まれる
#library(dplyr)   #library(tidyverse)で読み込まれる
library(gcookbook)

#パッケージアップデート
#update.packages()
#update.packages(ask = FALSE)

#テキスト読み込み
data <- read.csv("datafile.csv")
#data <- read_csv("datafile.csv") #readrパッケージ,高速
#ヘッダなし
#data <- read.csv("datafile.csv", header = FALSE) 
#手動で見出しを割り当てる
#names(data) <- c("Column1","Column2","Column3","Column4")
#デリミタ文字の指定
#data <- read.csv("datafile.csv", sep = "\t")  #タブ区切りテキスト
#data <- read.csv("datafile.csv", sep = " ")  #スペース区切りテキスト
#characterで読み込む
data <- read.csv("datafile.csv", stringsAsFactors = FALSE) 
#ファクタに変換
data$Sex <- factor(data$Sex)
str(data) #データフレームの要約表示
#read.table() #テキストデータの読み込み

#P6
#Excelファイル読み込み
#install.packages("readxl")
library(readxl)
data <- read_excel("datafile.xlsx", 1)
data2 <- read_excel("datafile.xlsx", sheet = 2) #シート2の読み込み
data3 <- read_excel("datafile.xlsx", sheet = "datafile3") #シート名で読み込み
data4 <- read_excel("datafile.xlsx", sheet = "datafile4", col_names = TRUE, col_types = c("blank", "text", "date", "text", "numeric")) 

#7
#SPSS/SAS/Stataファイルの読み込み
#install.packages("haven")
library(foreign)
#data <- read_sav("datafile.sav")  #SPSSファイルの読み込み
#data <- read_sas("datafile.sas")  #SASファイルの読み込み
#data <- read_dta("datafile.stata")  #Stataファイルの読み込み
#data <- read.octave() #OctaveファイルMATLABファイルの読み込み
#data <- read.systat() #SYSTATファイルの読み込み
#data <- read.xport() #SAS XPORTファイルの読み込み
#data <- read.spss() #SPSSファイルの読み込み

#P8
#パイプ演算子 library(magrittr)由来
library(dplyr)  #パイプ演算子magrittrをインポート
morley  #morleyデータセットの確認
morley %>%
  filter(Expt == 1) %>%  #Expt列=1のデータを抽出
  summary()  #データの要約表示
#summary(filter(morley, Expt == 1)) #パイプ演算子なし
f(x) #x %>% f()と同じ
h(g(f(x))) #以下と同じ
#x %>%
#  f() %>%
#  g() %>%
#  h()
#x <- x %>%
#  f() %>%
#  g() %>%
#  h()     #結果をxに保存
filter(morley, Expt == 1)
x <- filter(morley, Expt == 1) #Exptが1のデータを抽出
summary(x) #xの要約統計量
morley %>% filter(Expt == 1) #上記と同じ

#############################################################
#P83
#散布図
library(gcookbook) #hightweightデータセット読み込み
library(dplyr)
#2つの列の概要表示
heightweight %>%
  select(ageYear, heightIn)
ggplot(heightweight, aes(x = ageYear, y = heightIn)) +
  geom_point() #デフォルト
ggplot(heightweight, aes(x = ageYear, y = heightIn)) +
  geom_point(shape = 21) #白抜きプロット
ggplot(heightweight, aes(x = ageYear, y = heightIn)) +
  geom_point(size = 1.5) #プロットのsizeを2から1.5に小さくする
#使用する3つの列を指定する
heightweight %>%
  select(sex, ageYear, heightIn)
ggplot(heightweight, aes(x = ageYear, y = heightIn, color = sex)) +
  geom_point
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) +
  geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, color = sex)) +
  geom_point() #shapeとcolorをsexで割り当てる
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, color = sex)) +
  geom_point() +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_brewer(palette = "Set1") #shapeとcolorをsexで割り当て、任意の色と形を指定する


#P87#######################################################
#レシピ5.3 点の形を指定する
library(gcookbook) #hightweightデータセット読み込み
ggplot(heightweight, aes(x = ageYear, y = heightIn)) +
  geom_point(shape = 3) #shapeを一括で変更する
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) +
  geom_point(size = 3) #shapeをsexで割り当て、sizeを2から3に大きくする
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1, 4)) #shapeをsexで割り当て、sizeを2から3に大きくし、任意の形を指定する