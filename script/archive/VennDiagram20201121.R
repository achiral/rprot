#Proteome #VennDiagram
#https://www.kimoton.com/entry/2017/09/01/112812
##########################################################
rm(list = ls(all.names = TRUE))
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC") #作業ディレクトリ設定
#setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR") #作業ディレクトリ設定
getwd()
dir() 
##########################################################
library(readxl) #エクセル入力(read_excel)
#xlsx入力
P <- read_excel("DEPtwANOVA.xlsx", 1)
C <- read_excel("DEPtwANOVA.xlsx", 2)
PxC <- read_excel("DEPtwANOVA.xlsx", 3)
##########################################################
library(VennDiagram)
# グループを準備
PCP <- as.vector(P$Protein.IDs)
CLZ <- as.vector(C$Protein.IDs)
PCPxCLZ <- as.vector(PxC$Protein.IDs)
# データをリスト型に変換
dat1 <- list("PCP" = PCP, "CLZ" = CLZ, "PCPxCLZ" = PCPxCLZ)
# Chart
venn.diagram(
  x = dat1,
  #main = "AMY",
  #category.names = ,
  filename = "venn.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2),              # border line
  lty = c(1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1),              # border line color
  fill = c(4, 3, 2),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5 , 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1.5,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  #cat.dist = c(0.1, 0.1, 0.07), # カテゴリー名の枠線からの距離
  cat.dist = c(0.1, 0.1, 0.06),   # カテゴリー名の枠線からの距離NAc
  #cat.pos = c(315, 45, 180),    # カテゴリー名の位置の指定 (0-360度)
  cat.pos = c(0, 0, 0),          # カテゴリー名の位置の指定 (0-360度)NAc
  cat.default.pos = "outer",     # 
  rotation = 1,                  #
)
##########################################################
#7種類まで可能
#library(venn)
#dat2 <- list("PCP" = PCP, "CLZ" = CLZ, "PCPxCLZ" = PCPxCLZ)
#png("venn2.png", 300, 300, bg = "transparent")
#venn(dat2, ilab=TRUE, zcolor = "style", bg = "transparent")
#dev.off()
#pdf("venn2.pdf", width = 20/2.54, height = 20/2.54, useDingbats = FALSE)
#venn(dat2, ilab=TRUE, zcolor = "style", bg = "transparent")
#dev.off()
##########################################################
library(venneuler)
##########################################################
#重複部分の抽出
dup <- data.frame(intersect(intersect(PCP, CLZ), PCPxCLZ))
names(dup)[1] <- "Protein.IDs"
library(tidyverse)
dup2 <- left_join(dup, P, by = "Protein.IDs")
#xslx出力
library(xlsx)
smp <- list("intersect"=dup2)
write.xlsx(smp, "intersect.xlsx")





##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
#各部位の共通
library(readxl) #エクセル入力(read_excel)
##########################################################
rm(list = ls(all.names = TRUE))
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/AMY")
aP <- read_excel("DEPtwANOVA.xlsx", 1)
aC<- read_excel("DEPtwANOVA.xlsx", 2)
aPxC <- read_excel("DEPtwANOVA.xlsx", 3)
list_a <- list("PCP" = aP$Protein.IDs, "CLZ" = aC$Protein.IDs, "PCPxCLZ" = aPxC$Protein.IDs)
# Chart
venn.diagram(
  x = list_a,
  #main = "AMY",
  #category.names = ,
  filename = "venn.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2),              # border line
  lty = c(1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1),              # border line color
  fill = c(4, 3, 2),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5 , 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1.5,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.1, 0.1, 0.07), # カテゴリー名の枠線からの距離
  cat.pos = c(315, 45, 180),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  rotation = 1,                  #
)
##########################################################
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/HIP") #作業ディレクトリ設定
hP <- read_excel("DEPtwANOVA.xlsx", 1)
hC<- read_excel("DEPtwANOVA.xlsx", 2)
hPxC <- read_excel("DEPtwANOVA.xlsx", 3)
list_h <- list("PCP" = hP$Protein.IDs, "CLZ" = hC$Protein.IDs, "PCPxCLZ" = hPxC$Protein.IDs)
# Chart
venn.diagram(
  x = list_a,
  #main = "HIP",
  #category.names = ,
  filename = "venn.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2),              # border line
  lty = c(1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1),              # border line color
  fill = c(4, 3, 2),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5 , 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1.5,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.1, 0.1, 0.07), # カテゴリー名の枠線からの距離
  cat.pos = c(315, 45, 180),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  rotation = 1,                  #
)
##########################################################
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/NAc") #作業ディレクトリ設定
nP <- read_excel("DEPtwANOVA.xlsx", 1)
nC<- read_excel("DEPtwANOVA.xlsx", 2)
nPxC <- read_excel("DEPtwANOVA.xlsx", 3)
list_n <- list("PCP" = nP$Protein.IDs, "CLZ" = nC$Protein.IDs, "PCPxCLZ" = nPxC$Protein.IDs)
# Chart
venn.diagram(
  x = list_n,
  #main = "NAc",
  #category.names = ,
  filename = "venn.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2),              # border line
  lty = c(1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1),              # border line color
  fill = c(4, 3, 2),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5 , 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1.5,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.1, 0.1, 0.06),   # カテゴリー名の枠線からの距離NAc
  cat.pos = c(0, 0, 0),          # カテゴリー名の位置の指定 (0-360度)NAc
  cat.default.pos = "outer",     # 
  rotation = 1,                  #
)
##########################################################
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/PFC") #作業ディレクトリ設定
pP <- read_excel("DEPtwANOVA.xlsx", 1)
pC<- read_excel("DEPtwANOVA.xlsx", 2)
pPxC <- read_excel("DEPtwANOVA.xlsx", 3)
list_p <- list("PCP" = pP$Protein.IDs, "CLZ" = pC$Protein.IDs, "PCPxCLZ" = pPxC$Protein.IDs)
# Chart
venn.diagram(
  x = list_p,
  #main = "PFC",
  #category.names = ,
  filename = "venn.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2),              # border line
  lty = c(1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1),              # border line color
  fill = c(4, 3, 2),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5 , 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1.5,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.1, 0.1, 0.07), # カテゴリー名の枠線からの距離
  cat.pos = c(315, 45, 180),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  rotation = 1,                  #
)
##########################################################
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/STR") #作業ディレクトリ設定
sP <- read_excel("DEPtwANOVA.xlsx", 1)
sC<- read_excel("DEPtwANOVA.xlsx", 2)
sPxC <- read_excel("DEPtwANOVA.xlsx", 3)
list_s <- list("PCP" = sP$Protein.IDs, "CLZ" = sC$Protein.IDs, "PCPxCLZ" = sPxC$Protein.IDs)
# Chart
venn.diagram(
  x = list_s,
  #main = "STR",
  #category.names = ,
  filename = "venn.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2),              # border line
  lty = c(1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1),              # border line color
  fill = c(4, 3, 2),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5 , 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1.5,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.1, 0.1, 0.07), # カテゴリー名の枠線からの距離
  cat.pos = c(315, 45, 180),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  rotation = 1,                  #
)
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
setwd("/Users/user/Dropbox/0_Work/R/Perseus_Like_Analysis/ALL")
list_P <- list("AMY" = aP$Protein.IDs, "HIP" = hP$Protein.IDs, "NAc" = nP$Protein.IDs, "PFC" = pP$Protein.IDs, "STR" = sP$Protein.IDs)
list_C <- list("AMY" = aC$Protein.IDs, "HIP" = hC$Protein.IDs, "NAc" = nC$Protein.IDs, "PFC" = pC$Protein.IDs, "STR" = sC$Protein.IDs)
list_PxC <- list("AMY" = aPxC$Protein.IDs, "HIP" = hPxC$Protein.IDs, "NAc" = nPxC$Protein.IDs, "PFC" = pPxC$Protein.IDs, "STR" = sPxC$Protein.IDs)

list_P <- list("AMY" = aP$Protein.IDs, "HIP" = hP$Protein.IDs, "NAc" = nP$Protein.IDs, "PFC" = pP$Protein.IDs, "STR" = sP$Protein.IDs)
list_C <- list("AMY" = aC$Protein.IDs, "HIP" = hC$Protein.IDs, "NAc" = nC$Protein.IDs, "PFC" = pC$Protein.IDs, "STR" = sC$Protein.IDs)
list_PxC <- list("AMY" = aPxC$Protein.IDs, "HIP" = hPxC$Protein.IDs, "NAc" = nPxC$Protein.IDs, "PFC" = pPxC$Protein.IDs, "STR" = sPxC$Protein.IDs)

# Chart_PCP
venn.diagram(
  x = list_P,
  #main = "PCP",
  #category.names = ,
  filename = "venn_P.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2, 2, 2),              # border line
  lty = c(1, 1, 1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1, 1, 1),              # border line color
  fill = c(4, 3, 2, 5, 6),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.2, 0.25, 0.2, 0.23, 0.25), # カテゴリー名の枠線からの距離
  #cat.pos = c(288, 0, 72, 144, 216),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  #rotation = 1,                  #
)
# Chart_CLZ
venn.diagram(
  x = list_C,
  #main = "CLZ",
  #category.names = ,
  filename = "venn_C.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2, 2, 2),              # border line
  lty = c(1, 1, 1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1, 1, 1),              # border line color
  fill = c(4, 3, 2, 5, 6),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.2, 0.25, 0.2, 0.23, 0.25), # カテゴリー名の枠線からの距離
  #cat.pos = c(288, 0, 72, 144, 216),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  #rotation = 1,                  #
)
# Chart_PCPxCLZ
venn.diagram(
  x = list_PxC,
  #main = "PCPxCLZ",
  #category.names = ,
  filename = "venn_PxC.tiff",
  output = T,
  # Output features
  imagetype = "tiff",            # png, jpg, tiff
  height = 3000,                 # 高さ
  width = 3000,                  # 幅
  resolution = 600,              # 解像度
  compression = "lzw",
  margin = 0.1,                  # 余白
  # Circles
  lwd = c(2, 2, 2, 2, 2),              # border line
  lty = c(1, 1, 1, 1, 1),              # border line type, "blank"
  col = c(1, 1, 1, 1, 1),              # border line color
  fill = c(4, 3, 2, 5, 6),             # background color, myCol red(2)green(3)blue(4)sky(5)pink(6)yellow(7)
  alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),     # transparency
  #scaled = T,                   # 領域の面積に反映 (T) 
  #inverted = F,                 # 各領域の位置を逆転 (T),
  # Numbers
  label.col = 1,                 # label font color
  cex =1,                      # label font size
  fontface = 1,                  # 各領域のラベルのフォントデザイン"bold"
  fontfamily = "sans",           # 各領域のラベルのフォント
  # Set names
  cat.col = 1,                   # カテゴリー名の表示色
  cat.cex = 1.5,                 # カテゴリー名のフォントサイズ
  cat.fontface = 2,              # カテゴリー名のフォントデザイン"bold"
  cat.fontfamily = "sans",       # カテゴリー名のフォント
  cat.dist = c(0.2, 0.25, 0.2, 0.23, 0.25), # カテゴリー名の枠線からの距離
  #cat.pos = c(288, 0, 72, 144, 216),    # カテゴリー名の位置の指定 (0-360度)
  cat.default.pos = "outer",     # 
  #rotation = 1,                  #
)
##########################################################
#重複部分の抽出
dup_P <- data.frame(intersect
                    (intersect
                      (intersect
                        (intersect
                          (aP$Protein.IDs, hP$Protein.IDs),
                          nP$Protein.IDs),
                        pP$Protein.IDs),
                      sP$Protein.IDs))
names(dup_P)[1] <- "Protein.IDs"

dup_C <- data.frame(intersect
                    (intersect
                      (intersect
                        (intersect
                          (aC$Protein.IDs, hC$Protein.IDs),
                          nC$Protein.IDs),
                        pC$Protein.IDs),
                      sC$Protein.IDs))
names(dup_C)[1] <- "Protein.IDs"

dup_PxC <- data.frame(intersect
                      (intersect
                        (intersect
                          (intersect
                            (aPxC$Protein.IDs, hPxC$Protein.IDs),
                            nPxC$Protein.IDs),
                          pPxC$Protein.IDs),
                        sPxC$Protein.IDs))
names(dup_PxC)[1] <- "Protein.IDs"

library(tidyverse)
dup_P <- left_join(dup_P, aP, by = "Protein.IDs")
dup_C <- left_join(dup_C, aC, by = "Protein.IDs")
dup_PxC <- left_join(dup_PxC, aPxC, by = "Protein.IDs")

#xslx出力
library(xlsx)
#smp_PC <- list("intersect_P"=dup_P, "intersect_C"=dup_C, "intersect_PxC"=dup_PxC)
smp_P <- list("intersect_C"=dup_P)
smp_C <- list("intersect_C"=dup_C)
smp_PxC <- list("intersect_PxC"=dup_PxC)
#write.xlsx(smp_PC, "intersect_PC.xlsx")
write.xlsx(smp_P, "intersect_P.xlsx")
write.xlsx(smp_C, "intersect_C.xlsx")
write.xlsx(smp_PxC, "intersect_PxC.xlsx")
##########################################################