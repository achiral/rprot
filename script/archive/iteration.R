#iteration
#イテレーション(反復)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ASSIGN")
library(assign)
#aに1を代入する
assign("a",1)
a
#文字をつなぐ
c <- "ゴレライ"
paste("ラッスン",c,sep="") 
#bに1〜10を繰り返し代入する
b <- "x"
assign(b,1:10)
x
#
for(i in 1:10){
  nam <- paste("y",i,sep="")
  assign(nam,i*i)
}
y1;y2;y3;y8;y9;y10


#SCZ#########################################################
#CSVファイル読み込み
data <- read.csv("GxInfo4.csv", stringsAsFactors = T, fileEncoding = "UTF-8-BOM") #CSVの読み込み
str(data)
#############################################################
#変数のクラス確認
sapply(data, class)
#characterをfactorに変換
data[,1] <- as.factor(data[,1])
class(data[,1])
#for関数による変換
colc <- grep("character", sapply(data,class))
for(i in 1:length(colc)){
  data[,colc[i]] <- as.factor(data[,colc[i]])
}
#データの分割
colf <- grep("factor",sapply(data, class)) #列番号で返す
colf
datan <- data4[,-colf] #factor列を除く
dataf <- data4[,colf] #factor列
#plotとhist
plot(dataf[,列番号],
     main = colname(dataf)[列番号],
     cex.main =3 #グラフタイトル文字サイズ
     col = "black", #バーの色
     space = 0.3) #バーの間隔
hist(datan[,列番号],
     cex.main =3 #グラフタイトル文字サイズ
     xlab = NULL, #x軸ラベル
     ylab = NULL, #y軸ラベル
     breaks = 20, #階級の分割数
     col = "gray") #バーの色
#for関数
for(i in 1:ncol(dataf)){
  plot(dataf[,i])
}
リストの中身:(1,2,3,4)
for(i in 1:ncol(dataf)){
  plot(dataf[,i],
       main = colnames(dataf)[i],
       col ="black",
       space = 0.3,
       cex.main = 3)
}
for(i in 1:ncol(datan)){
  hist(datan[,i],
       main = colnames(datan)[i],
       cex.main = 3,
       xlab = NULL,
       ylab = NULL,
       breaks = 20,
       col = "gray")
}


for(i in 1:ncol(data)){
  result[[i]] <- lm(ddCt~data[,i],data = data)
}
