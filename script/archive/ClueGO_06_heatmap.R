#ClueGO_06_heatmap.R
#https://note.com/shigeruhanano/n/naca751f9a547
#https://datator.exblog.jp/23736047/
##https://www.datanovia.com/en/blog/how-to-normalize-and-standardize-data-in-r-for-great-heatmap-visualization/#heatmap-of-the-raw-data
#https://www.r-graph-gallery.com/215-the-heatmap-function.html
#############################################################
setwd("/Users/user/Dropbox/0_Work/R/SWATH") #作業ディレクトリ設定
getwd() #作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
#############################################################
#パッケージの読み込み
library(EnhancedVolcano)
#library(airway)
library(magrittr)
#library("DESeq2")
library(tidyverse) #ggplot2,dplyr使用
#tidyverse(ggplot2,dplyr),gcookbook読み込み
#library(dplyr)
#library(stringr)
library(scales) #muted()関数使用のため
#library(ggcorrplot) #グラフ作成のため
#library(cowplot)
library(rJava)
library(genefilter) #ヒートマップ
library(gplots) #ヒートマップ
library(ComplexHeatmap) #ヒートマップ
library(RColorBrewer) #色
library(openxlsx) #エクセル入出力(write.xlsx)
#library(XLConnect) #エクセル入出力,動作しない
#library(xlsx) #エクセル出力
#library(xlsx2) #エクセル出力,動作しない
library(readxl) #エクセル入力(read_excel)
#library(tablaxlsx) #エクセル表出力
library(gridExtra) #svg出力
#library(rvg) #svg出力
#library(rsvg) #svg出力
library(sets) #集合演算
#############################################################
#前処理：データが測定値の対数→対数を戻して平均し、log2変換→済
#txt入力
dat1 <- read.table("dat1.txt", header=TRUE, sep="\t", row.names=1) #行名入力
dat2 <- read.delim("dat2.txt", header=TRUE, row.names=1) #行名入力
dat3 <- read.delim("dat3.txt", header=TRUE, row.names=1) #行名入力
dat4 <- read.delim("dat4.txt", header=TRUE, row.names=1) #行名入力
dat5 <- read.delim("dat5.txt", header=TRUE, row.names=1) #行名入力
dat6 <- read.delim("dat6.txt", header=TRUE, row.names=1) #行名入力
#xlsx入力
#dat1 <- read_excel("heatCL.xlsx", 1) #シート1の読み込み
#dat2 <- read_excel("heatCL.xlsx", 2) #シート2の読み込み
#dat3 <- read_excel("heatCL.xlsx", 3) #シート3の読み込み
#dat4 <- read_excel("heatCL.xlsx", 4) #シート4の読み込み
#dat5 <- read_excel("heatCL.xlsx", 5) #シート5の読み込み
#dat6 <- read_excel("heatCL.xlsx", 6) #シート6の読み込み
#############################################################
#Raw data
png("hmr1.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
#heatmap(as.matrix(dat1), margin=c(4,8), main="C1 (Raw Data)")
heatmap.2(as.matrix(dat1), margin=c(4,8), 
          main="C1 (Raw Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmr2.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat2), margin=c(4,8), 
          main="C2 (Raw Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmr3.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat3), margin=c(4,8), 
          main="C3 (Raw Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmr4.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat4), margin=c(4,8), 
          main="C4 (Raw Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmr5.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat5), margin=c(4,8), 
          main="C5 (Raw Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()
#############################################################
#Z-score化
dat1z <- genescale(dat1, axis=1, method="Z")
dat2z <- genescale(dat2, axis=1, method="Z")
dat3z <- genescale(dat3, axis=1, method="Z")
dat4z <- genescale(dat4, axis=1, method="Z")
dat5z <- genescale(dat5, axis=1, method="Z")

png("hmz1.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat1z), margin=c(4,8), 
          main="C1 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmz2.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat2z), margin=c(4,8), 
          main="C2 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmz3.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat3z), margin=c(4,8), 
          main="C3 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmz4.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat4z), margin=c(4,8), 
          main="C4 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()

png("hmz5.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat5z), margin=c(4,8), 
          main="C5 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none")
dev.off()
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#txt入力
dat1 <- read.table("dat1.txt", header=TRUE, sep="\t") #行名入力しない
dat2 <- read.delim("dat2.txt", header=TRUE) #行名入力しない
dat3 <- read.delim("dat3.txt", header=TRUE) #行名入力しない
dat4 <- read.delim("dat4.txt", header=TRUE) #行名入力しない
dat5 <- read.delim("dat5.txt", header=TRUE) #行名入力しない
dat6 <- read.delim("dat6.txt", header=TRUE) #行名入力しない
#############################################################
#Average(SC0)
dat1$SC0_ave <- rowMeans(dat1[,2:5]) #2-5列目の平均列を作成
dat2$SC0_ave <- rowMeans(dat2[,2:5]) #2-5列目の平均列を作成
dat3$SC0_ave <- rowMeans(dat3[,2:5]) #2-5列目の平均列を作成
dat4$SC0_ave <- rowMeans(dat4[,2:5]) #2-5列目の平均列を作成
dat5$SC0_ave <- rowMeans(dat5[,2:5]) #2-5列目の平均列を作成
dat6$SC0_ave <- rowMeans(dat6[,2:5]) #2-5列目の平均列を作成
#############################################################
#Fold change(vs SC0_ave)
m <- mutate(dat1, SC10_1fc = SC10_1 / SC0_ave)
m <- mutate(m, SC10_2fc = SC10_2 / SC0_ave)
m <- mutate(m, SC10_3fc = SC10_3 / SC0_ave)
m <- mutate(m, SC10_4fc = SC10_4 / SC0_ave)
m <- mutate(m, PC0_1fc = PC0_1 / SC0_ave)
m <- mutate(m, PC0_2fc = PC0_2 / SC0_ave)
m <- mutate(m, PC0_3fc = PC0_3 / SC0_ave)
m <- mutate(m, PC0_4fc = PC0_4 / SC0_ave)
m <- mutate(m, PC10_1fc = PC10_1 / SC0_ave)
m <- mutate(m, PC10_2fc = PC10_2 / SC0_ave)
m <- mutate(m, PC10_3fc = PC10_3 / SC0_ave)
dat1 <- mutate(m, PC10_4fc = PC10_4 / SC0_ave)

m <- mutate(dat2, SC10_1fc = SC10_1 / SC0_ave)
m <- mutate(m, SC10_2fc = SC10_2 / SC0_ave)
m <- mutate(m, SC10_3fc = SC10_3 / SC0_ave)
m <- mutate(m, SC10_4fc = SC10_4 / SC0_ave)
m <- mutate(m, PC0_1fc = PC0_1 / SC0_ave)
m <- mutate(m, PC0_2fc = PC0_2 / SC0_ave)
m <- mutate(m, PC0_3fc = PC0_3 / SC0_ave)
m <- mutate(m, PC0_4fc = PC0_4 / SC0_ave)
m <- mutate(m, PC10_1fc = PC10_1 / SC0_ave)
m <- mutate(m, PC10_2fc = PC10_2 / SC0_ave)
m <- mutate(m, PC10_3fc = PC10_3 / SC0_ave)
dat2 <- mutate(m, PC10_4fc = PC10_4 / SC0_ave)

m <- mutate(dat3, SC10_1fc = SC10_1 / SC0_ave)
m <- mutate(m, SC10_2fc = SC10_2 / SC0_ave)
m <- mutate(m, SC10_3fc = SC10_3 / SC0_ave)
m <- mutate(m, SC10_4fc = SC10_4 / SC0_ave)
m <- mutate(m, PC0_1fc = PC0_1 / SC0_ave)
m <- mutate(m, PC0_2fc = PC0_2 / SC0_ave)
m <- mutate(m, PC0_3fc = PC0_3 / SC0_ave)
m <- mutate(m, PC0_4fc = PC0_4 / SC0_ave)
m <- mutate(m, PC10_1fc = PC10_1 / SC0_ave)
m <- mutate(m, PC10_2fc = PC10_2 / SC0_ave)
m <- mutate(m, PC10_3fc = PC10_3 / SC0_ave)
dat3 <- mutate(m, PC10_4fc = PC10_4 / SC0_ave)

m <- mutate(dat4, SC10_1fc = SC10_1 / SC0_ave)
m <- mutate(m, SC10_2fc = SC10_2 / SC0_ave)
m <- mutate(m, SC10_3fc = SC10_3 / SC0_ave)
m <- mutate(m, SC10_4fc = SC10_4 / SC0_ave)
m <- mutate(m, PC0_1fc = PC0_1 / SC0_ave)
m <- mutate(m, PC0_2fc = PC0_2 / SC0_ave)
m <- mutate(m, PC0_3fc = PC0_3 / SC0_ave)
m <- mutate(m, PC0_4fc = PC0_4 / SC0_ave)
m <- mutate(m, PC10_1fc = PC10_1 / SC0_ave)
m <- mutate(m, PC10_2fc = PC10_2 / SC0_ave)
m <- mutate(m, PC10_3fc = PC10_3 / SC0_ave)
dat4 <- mutate(m, PC10_4fc = PC10_4 / SC0_ave)

m <- mutate(dat5, SC10_1fc = SC10_1 / SC0_ave)
m <- mutate(m, SC10_2fc = SC10_2 / SC0_ave)
m <- mutate(m, SC10_3fc = SC10_3 / SC0_ave)
m <- mutate(m, SC10_4fc = SC10_4 / SC0_ave)
m <- mutate(m, PC0_1fc = PC0_1 / SC0_ave)
m <- mutate(m, PC0_2fc = PC0_2 / SC0_ave)
m <- mutate(m, PC0_3fc = PC0_3 / SC0_ave)
m <- mutate(m, PC0_4fc = PC0_4 / SC0_ave)
m <- mutate(m, PC10_1fc = PC10_1 / SC0_ave)
m <- mutate(m, PC10_2fc = PC10_2 / SC0_ave)
m <- mutate(m, PC10_3fc = PC10_3 / SC0_ave)
dat5 <- mutate(m, PC10_4fc = PC10_4 / SC0_ave)

m <- mutate(dat6, SC10_1fc = SC10_1 / SC0_ave)
m <- mutate(m, SC10_2fc = SC10_2 / SC0_ave)
m <- mutate(m, SC10_3fc = SC10_3 / SC0_ave)
m <- mutate(m, SC10_4fc = SC10_4 / SC0_ave)
m <- mutate(m, PC0_1fc = PC0_1 / SC0_ave)
m <- mutate(m, PC0_2fc = PC0_2 / SC0_ave)
m <- mutate(m, PC0_3fc = PC0_3 / SC0_ave)
m <- mutate(m, PC0_4fc = PC0_4 / SC0_ave)
m <- mutate(m, PC10_1fc = PC10_1 / SC0_ave)
m <- mutate(m, PC10_2fc = PC10_2 / SC0_ave)
m <- mutate(m, PC10_3fc = PC10_3 / SC0_ave)
dat6 <- mutate(m, PC10_4fc = PC10_4 / SC0_ave)
#############################################################
#列抽出
t(names(dat1)) #行列の名前確認
dat1 <- dat1[,c(1,19:30)] #列抽出
dat2 <- dat2[,c(1,19:30)]  #列抽出
dat3 <- dat3[,c(1,19:30)]  #列抽出
dat4 <- dat4[,c(1,19:30)]  #列抽出
dat5 <- dat5[,c(1,19:30)]  #列抽出
dat6 <- dat6[,c(1,19:30)]  #列抽出
#列名用txt入力
dat <- read.table("dat1.txt", header=TRUE, sep="\t") #行名入力しない
t(names(dat)) #行列の名前確認
dat <- dat[,c(1,6:17)] #列抽出
# Copy the column names from dat to dat1-6:
colnames(dat1) <- colnames(dat)
colnames(dat2) <- colnames(dat)
colnames(dat3) <- colnames(dat)
colnames(dat4) <- colnames(dat)
colnames(dat5) <- colnames(dat)
colnames(dat6) <- colnames(dat)
#############################################################
#txt出力
write.table(dat1, file="dat1fc.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat2, file="dat2fc.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat3, file="dat3fc.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat4, file="dat4fc.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat5, file="dat5fc.txt", sep="\t", row.names=F) #行名出力しない
write.table(dat6, file="dat6fc.txt", sep="\t", row.names=F) #行名出力しない
#############################################################
#txt入力
dat1 <- read.table("dat1fc.txt", header=TRUE, sep="\t", row.names=1) #行名入力
dat2 <- read.delim("dat2fc.txt", header=TRUE, row.names=1) #行名入力
dat3 <- read.delim("dat3fc.txt", header=TRUE, row.names=1) #行名入力
dat4 <- read.delim("dat4fc.txt", header=TRUE, row.names=1) #行名入力
dat5 <- read.delim("dat5fc.txt", header=TRUE, row.names=1) #行名入力
dat6 <- read.delim("dat6fc.txt", header=TRUE, row.names=1) #行名入力
#############################################################
#Raw data(Log2 FC)
png("hmr1fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat1), margin=c(4,8), #col = greenred(75), colオプション上手く行かない
          scale="row", #Normarize
          key=T, keysize=1.5, #不明
          main="C1 (Log2 FC)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=0.9, cexRow=0.5
)
dev.off()

png("hmr2fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat2), margin=c(4,8), #col = greenred(75), colオプション上手く行かない
          scale="row", #Normarize
          key=T, keysize=1.5, #不明
          main="C2 (Log2 FC)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=0.9, cexRow=0.5
)
dev.off()

png("hmr3fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat3), margin=c(4,8), #col = greenred(75), colオプション上手く行かない
          scale="row", #Normarize
          key=T, keysize=1.5, #不明
          main="C3 (Log2 FC)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=0.9, cexRow=0.5
)
dev.off()

png("hmr4fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat4), margin=c(4,8), #col = greenred(75), colオプション上手く行かない
          scale="row", #Normarize
          key=T, keysize=1.5, #不明
          main="C4 (Log2 FC)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=0.9, cexRow=0.5
)
dev.off()

png("hmr5fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat5), margin=c(4,8), #col = greenred(75), colオプション上手く行かない
          scale="row", #Normarize
          key=T, keysize=1.5, #不明
          main="C5 (Log2 FC)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=0.9, cexRow=0.5
)
dev.off()

png("hmr6fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat6), margin=c(4,8), #col = greenred(75), colオプション上手く行かない
          scale="row", #Normarize
          key=T, keysize=1.5, #不明
          main="All (Log2 FC)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=0.9, cexRow=0.5
)
dev.off()
#############################################################
#Z-score化
dat1z <- genescale(dat1, axis=1, method="Z")
dat2z <- genescale(dat2, axis=1, method="Z")
dat3z <- genescale(dat3, axis=1, method="Z")
dat4z <- genescale(dat4, axis=1, method="Z")
dat5z <- genescale(dat5, axis=1, method="Z")
dat6z <- genescale(dat6, axis=1, method="Z")

Colv= TRUE
png("hmz1fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat1z), 
          main="C1 (Z score Data)", #グラフタイトル
          Rowv = TRUE, #Fだと樹形図書けない
          Colv= NA, #Fだと樹形図書けないColv=if(symm)"Rowv" else TRUE, #Colだとサンプル順番ずれる
          distfun = dist, hclustfun = hclust, #樹形図設定
          dendrogram = c("row"), #"both","row","column","none" #Colだとサンプル順番ずれる
          #reorderfun = function(d, w) reorder(d, w), #不明
          #symm = FALSE, #不明
          scale = c("row"), #scale="row", #Normarize "none","row", "column"
          #na.rm=TRUE,
          #revC = identical(Colv, "Rowv"), # image plot
          #add.expr, #不明 # image plot
          #breaks, # mapping data to colors
          #symbreaks=any(x < 0, na.rm=TRUE) || scale!="none", # mapping data to colors
          col = greenred(75), #col = brewer.pal(11,"RdBu"), #col="heat.colors",
          #colsep, # block sepration
          #rowsep, # block sepration
          #sepcolor="white", # block sepration
          #sepwidth=c(0.05,0.05), # block sepration
          #cellnote, # cell labeling
          #notecex=1.0, # cell labeling
          #notecol="cyan", # cell labeling
          #na.color=par("bg"), # cell labeling
          trace=c("none"), #"column","row","both","none", トレース無し
          #tracecol="cyan", #トレース色
          #hline=median(breaks), # level trace
          #vline=median(breaks), # level trace
          #linecol=tracecol, # level trace
          margins = c(5,5), # Row/Column Labelingラベル記載スペースの設定
          #ColSideColors, # Row/Column Labeling
          #RowSideColors, # Row/Column Labeling
          cexRow=1, cexCol=1, # Row/Column Labeling
          #labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, # Row/Column Labeling
          #adjRow = c(0,NA), adjCol = c(NA,0), offsetRow = 0.5, offsetCol = 0.5, # Row/Column Labeling
          #colRow = NULL, colCol = NULL, # Row/Column Labeling
          key = TRUE, keysize = 1.5, # color key + density info
          #density.info=c("none"), #"histogram","density","none" # color key + density info
          #denscol=tracecol, # color key + density info
          #symkey = any(x < 0, na.rm=TRUE) || symbreaks, # color key + density info
          #densadj = 0.25, # color key + density info
          #key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, # color key + density info
          #key.ytickfun = NULL, # color key + density info
          #key.par=list(), # color key + density info
          #main = NULL, xlab = NULL, ylab = NULL, # plot labels
          #lmat = NULL, lhei = NULL, lwid = NULL, # plot layout
          #extrafun=NULL # extras
)
dev.off()


png("hmz1fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat1z), margin=c(5,5), 
          scale="row", #Normarize 
          key=T, keysize=1.5, #不明
          main="C1 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col = greenred(75), #col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=1, cexRow=1
)
dev.off()

png("hmz2fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat2z), margin=c(5,5), 
          scale="row", #Normarize 
          key=T, keysize=1.5, #不明
          main="C2 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col = greenred(75), #col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=1, cexRow=1
)
dev.off()

png("hmz3fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat3z), margin=c(5,5),
          scale="row", #Normarize 
          key=T, keysize=1.5, #不明
          main="C3 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col = greenred(75), #col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=1, cexRow=1
)
dev.off()

png("hmz4fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat4z), margin=c(5,5),
          scale="row", #Normarize 
          key=T, keysize=1.5, #不明
          main="C4 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col = greenred(75), #col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=1, cexRow=1
)
dev.off()

png("hmz5fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat5z), margin=c(5,5), 
          scale="row", #Normarize 
          key=T, keysize=1.5, #不明
          main="C5 (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col = greenred(75), #col=brewer.pal(11,"RdBu"),  
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=1, cexRow=1
)
dev.off()

png("hmz6fc.png", height = 2400, width = 2400, res = 300) ##ファイル名、画像サイズ、解像度を指定
heatmap.2(as.matrix(dat6z), margin=c(5,5), 
          scale="row", #Normarize 
          key=T, keysize=1.5, #不明
          main="All (Z score Data)", #グラフタイトル
          Rowv = TRUE, Colv= NA, 
          distfun = dist, hclustfun = hclust, 
          dendrogram = c("row"), col = greenred(75), #col=brewer.pal(11,"RdBu"), 
          trace="none", #density.info="none", #densityは注釈の分布
          cexCol=1, cexRow=1
)
dev.off()
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#limmaパッケージによるt検定
# Draw a boxplot.
boxplot(dat1, cex.axis=0.8, las=2, main="Original distribution of data",
        ylab="Log2(Intensity)") 
# Normalization
library(preprocessCore)
dat1m <- normalize.quantiles(as.matrix(dat1)) # A quantile normalization. 
# Copy the row and column names from data to data2:
rownames(dat1m) <- rownames(dat1)
colnames(dat1m) <- colnames(dat1)
# Draw a boxplot.
boxplot(dat1m, cex.axis=0.8, las=2, main="Distribution after normalization",
        ylab="Log2(Intensity)")
# t-test using the limma package:
library(limma)
design = cbind(SC10 = c(1,1,1,1,0,0,0,0,0,0,0,0), # First 4 columns->SC10
               PC0 = c(0,0,0,0,1,1,1,1,0,0,0,0),
               PC10 = c(0,0,0,0,0,0,0,0,1,1,1,1)) # Last 4 columns->PC10
fit <- lmFit(dat1m, design=design) # Fit the original matrix to the above design.
# We want to compare SC10 vs. PC0, SC10 vs. PC10 and PC0 vs. PC10
contrastsMatrix <- makeContrasts("SC10-PC0","SC10-PC10","PC0-PC10",
                                 levels = design) 
fit2 <- contrasts.fit(fit, contrasts = contrastsMatrix) # Making the comparisons.
fit2 <- eBayes(fit2) # Moderating the t-tetst by eBayes method.
df <- data.frame(fit2)
#############################################################