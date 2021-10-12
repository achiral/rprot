#Stepchart
#レベル順序(https://www.jaysong.net/RBook/factor.html)
#最大値,最小値の位置（https://bioinfo-dojo.net/2017/08/31/r-which_max/）
#ベクトルからのオブジェクト抽出（https://siguniang.wordpress.com/2010/10/04/last_element_of_vecor_in_r/）
#列名変更(https://indenkun.hatenablog.com/entry/2020/06/21/003000)
####################################################################################################################
setwd("~/Dropbox/0_Work/R/Prescription/dose")
df <- read.csv("dose.csv")
df2 <- read.csv("dose2.csv")
####################################################################################################################
library(tidyverse)
####################################################################################################################
df2 <- read.csv("dose2.csv")
df3 <- df2 %>% gather(key = AP, value = "dose", "CPZeq", "APZ", "OLZ", "RIS", "QTP", "BNS", "CPZ", "HPD", "PER", "CLZ", "TMP", "PAR", "ASE", "BRE", "ZTP")
uniID <- unique(df3$ID)

df3 <- df3 %>% mutate(AP = fct_inorder(AP)) #登場順でfactorのlevelsを再整理
levels(df3$AP)
unique(df3$AP)

my_theme <- theme_bw() +
  theme(
    text = element_text(size = 16, color = "black", family = "Helvetica", face = "bold"),
    axis.title = element_text(size = rel(2))
  )


pdf(file="step.pdf", width=6.4,height=6.4, family="Japan1")
for(i in 1:length(uniID)){
  df4 <- df3 %>% filter(ID == uniID[i]) %>% filter(!is.na(dose))
  p<-ggplot(df4,aes(day, dose)) + 
    geom_step(aes(color= AP), alpha = 0.7, size = 2) +
    geom_point(aes(shape = AP, color = AP), size = 1.5, alpha = 0.9) +
    my_theme
  print(p)
}
dev.off()
####################################################################################################################
#減量スクリーニング
setwd("~/Dropbox/0_Work/R/Prescription/dose")
df <- read.csv("dose.csv")
library(tidyverse)
MAX <- apply(df[,-1], MARGIN = 2, function(x) max(x, na.rm=TRUE))
MIN <- apply(df[,-1], MARGIN = 2, function(x) min(x, na.rm=TRUE))
MAX_d <- apply(df[,-1], MARGIN = 2, function(x) {max(which(x==max(x, na.rm=TRUE)))}) #最大用量だった最後の日
MIN_d <- apply(df[,-1], MARGIN = 2, function(x) {which(x==min(x, na.rm=TRUE))[[1]]}) #最小用量だった最初の日
RATE <- (MAX-MIN)/(MIN_d-MAX_d)

#増量スクリーニング
MAX_d2 <- apply(df[,-1], MARGIN = 2, function(x) {which(x==max(x, na.rm=TRUE))[[1]]}) #最大用量だった最初の日
MIN_d2 <- apply(df[,-1], MARGIN = 2, function(x) {max(which(x==min(x, na.rm=TRUE)))}) #最小用量だった最後の日
RATE2 <- (MAX-MIN)/(MAX_d2-MIN_d2)

#データの整形
df4 <- as.data.frame(rbind(MAX,MIN,MAX_d,MIN_d,RATE,MAX_d2,MIN_d2,RATE2))
df5 <- as.data.frame(cbind(MAX,MIN,MAX_d,MIN_d,RATE,MAX_d2,MIN_d2,RATE2))
df6 <- cbind(rownames(df5),df5)
colnames(df6)[1] <- "ID"
df7 <- df6 %>% separate(ID, c("ID", "AP"), sep="_")

#output xlsx
library(openxlsx) #入出力(write.xlsx)
smp <- list("dose"=df,"dose2"=df2,"gather"=df3,"rate1"=df4,"rate2"=df7) #リスト作成
write.xlsx(smp, "data.xlsx")
####################################################################################################################