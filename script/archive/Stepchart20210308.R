#Stepchart
#レベル順序(https://www.jaysong.net/RBook/factor.html)
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

df3 <- df3 %>% mutate(AP = fct_inorder(AP))
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









plot(c(df[,1],df[,1],df[,1],df[,1],df[,1]),c(df[,2],df[,3],df[,4],df[,5],df[,6]),type="s",col=c("red","blue","green","magenta","black"))

postscript("dose.eps", horizontal = F, onefile = F, paper = "special", height = 4, width = 4)
plot(df[,2],xlab="",ylab="",type="s",lwd="2",col="blue")
points(df[,3],type="s",lwd="2",col="magenta")
dev.off()

postscript("dose.eps", horizontal = F, onefile = F, paper = "special", height = 10, width = 10)
plot(df$S001_CP,xlab="",ylab="",type="s",col="blue")
points(df$S001_APZ,type="s",col="magenta")
points(df$S001_OLZ,type="s",col="magenta")
dev.off()


postscript("dose.eps", horizontal = F, onefile = F, paper = "special", height = 10, width = 10)
plot(df$S001_RIS,xlab="",ylab="",type="s",col="blue")
points(df$S001_APZ,type="s",col="magenta")
points(df$S001_OLZ,type="s",col="magenta")
dev.off()






library(ggplot2)
set.seed(10)
df=data.frame(id=1:100,y=rnorm(100))


xx <- grep("S001",colnames(df))

fun(x) <- function(x){x = for(i in xx)}
xxx <- ggplot(df)+geom_step(aes(day,df[,i]))


par(mfrow=c(2,2))
x <- list(1:6, sin, runif(6), dnorm)
for (i in x) {
  plot(i, xlim=c(0,2*pi))
}

par(mfrow=c(2,2))

x <- list(1:6, sin, runif(6), dnorm)
for (i in xx) {
  plot(i)
}

xxx
xx <- grep("_CP",colnames(df))




