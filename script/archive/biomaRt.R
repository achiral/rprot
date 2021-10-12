#biomaRt
#biomaRtインストール
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("biomaRt")

#biomaRt開く
browseVignettes("biomaRt")
library(biomaRt)

#遺伝子にGeneOntologyのアノテーションをつける
ensid <- c("ADH1B", "ALDH2")
db <- useMart("ensembl")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("hgnc_symbol", "go_id", "name_1006"),filters = "hgnc_symbol",values = ensid,mart = hg)
res

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RamiGO")

#作業ディレクトリ設定
setwd("/Users/user/Dropbox/0_Work/R/Directory_DEF")
#作業ディレクトリ確認
getwd()
#作業ディレクトリ内のファイル表示
dir() 
#データの読み込み
res <- read.table("Stat_DEF2h.txt", header=TRUE, stringsAsFactors=F)
head(res)
#データの読み込み
res3 <- read.table("Stat_LCL.txt", header=TRUE, stringsAsFactors=F)
head(res3)




#実数に変換
p.num <- as.numeric(res$p) 
res[,"p.num"] <- p.num
FC.num <- as.numeric(res$FC) 
res[,"FC.num"] <- FC.num
res2 <- read.table("Stat_DEF24h.txt", header=TRUE, stringsAsFactors=F)
head(res2)
# Make a basic volcano plot
with(res3, plot(log2FC, -log10(p), pch=20, main="Volcano plot", xlim=c(-2,2)))
# Add colored points: red if p<0.05, orange of log2FC>0.58, green if both)
with(subset(res3, p<.05 ), points(log2FC, -log10(p), pch=20, col="red"))
with(subset(res3, abs(log2FC)>0.58), points(log2FC, -log10(p), pch=20, col="orange"))
with(subset(res3, p<.05 & abs(log2FC)>0.58), points(log2FC, -log10(p), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res3, p<.05 & abs(log2FC)>0.58), textxy(log2FC, -log10(p), labs=Gene, cex=.8))
with(subset(res3, p<.05 & abs(log2FC)>1), textxy(log2FC, -log10(p), labs=Gene, cex=.8))



#オリジナル
res <- read.table("results.txt", header=TRUE)
head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))




#クリップボードからのデータ読み込み
dat <- read.table(pipe("pbpaste"), header=TRUE)  #Macの場合
# header=TRUEで変数名を含む，FALSEにすると含まない
attach(dat)    # dat$COLUMN などと書かなくていいようにattach
cor(Log2FC_DEF2h, Log2FC_DEF24h)  # 何も指定しなければピアソンの積率相関係数
#相関行列表（correlation matrix）を出す
cor(cbind(dat))
#相関係数以外（p値など）も出す
cor.test(Log2FC_DEF2h, Log2FC_DEF24h)
#スピアマンの順位相関係数
cor(Log2FC_DEF2h, Log2FC_DEF24h, method="spearman")
#ケンドールの順位相関係数（データが少なくて，同率順位が多いとき）
cor(Log2FC_DEF2h, Log2FC_DEF24h, method="kendall")

#散布図
par(family="HiraKakuPro-W3") # Macで日本語表記する
plot(Log2FC_DEF2h, Log2FC_DEF24h, xlim=c(-5,5), ylim=c(-5,5)) # xlimとylimで範囲を指定
# 回帰直線を入れる場合は以下を追加
abline(lm(Log2FC_DEF2h~Log2FC_DEF24h), col="red")

#psychパッケージを使って散布図作成（相関係数，ヒストグラム，回帰直線）
#install.packages("psych")　でインストールしておく
library(psych)                # psychパッケージ
par(family="HiraKakuPro-W3")  # Macで日本語表示する
pairs.panels(dat)

#コピペから散布図
pairs(dat)
pairs(read.table(pipe("pbpaste"), header=TRUE))  #Macの場合
pairs.panels(dat,hist.col="white",rug=F,ellipses=F,lm=T)
cor(dat)

#無向グラフ
install.packages("qgraph")
library(qgraph)
qgraph(cor(dat),edge.labels=T)
qgraph(cor(dat),edge.labels=T,minimum=.2,edge.color="black")

#ヒートマップ
cor.plot(cor(dat2))
cor.plot(cor(dat2),numbers=T)

#主座標分析（古典的多次元尺度構成法）的な可視化
plot(cmdscale(dist(cor(dat3))),type="n",xlab="", ylab=""); text(cmdscale(dist(cor(dat3))),rownames(cor(dat3)))
#コリログラム
install.packages("corrgram")
library(corrgram)
#corrgram(dat3, upper.panles = XXXX, lower.panel = YYYY)
##panel.pie：相関係数の大きさを示す色付き円グラフ
##panel.shade：相関係数を色で示す
##panel.ellipse：相関楕円と当てはめ
##panel.pts ：散布図
corrgram(dat, upper.panles = panel.ellipse, lower.panel = panel.pie)
corrgram(dat,lower.panel=panel.ellipse,upper.panel=panel.pie)
corrgram(dat,lower.panel=panel.pts,upper.panel=panel.shade)

