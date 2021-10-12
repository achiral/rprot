#Tukey HSD
#How to perform Tukey HSD.test() on list of dataframes?
#https://stackoverflow.com/questions/49168959/how-to-perform-tukey-hsd-test-on-list-of-dataframes
##############################################################
setwd("/Users/user/Dropbox/0_Work/R/TukeyHSD") #作業ディレクトリ設定
getwd()#作業ディレクトリ確認
dir() #作業ディレクトリ内のファイル表示
##############################################################
#パッケージインストール
#install.packages("multcomp") #多重比較検定
#ライブラリ読み込み
library(tidyverse) #ライブラリtidyverse(ggplot2,dplyr),gcookbook読み込み
library(readxl) #エクセル入力(read_excel)
library(openxlsx) #エクセル入出力(write.xlsx)
library(multcomp) #多重比較検定
require(agricolae) #
##############################################################
#xlsx入力
rm(list = ls(all.names = TRUE))
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


category <- c(rep("young", 3), rep("Middle", 4), rep("old", 5))
fat <- c(1857.87, 1953.90, 1440.70, 1553.81, 1785.91, 1893.82, 1483.75, 1784.99, 2011.01, 2023.04, 2011.05, 1788.81)
BMI <- c(21.1, 23.2, 24.5, 25.6, 21.8, 18.0, 19.2, 20.1, 22.1, 25.0, 26.1, 25.1)
age <- c(25, 23, 27, 55, 58, 62, 45, 75, 80, 75, 83, 89)
df <- data.frame(fat, BMI, age, category)
df2 <- data.frame(fat, BMI, age, category)

lm.fat <- (lm(fat ~ as.factor(category), data = df2))
anova(lm.fat)
require(agricolae)
HSD.test(lm.fat, "as.factor(category)", group = TRUE, console = TRUE)

an <- lapply(df, function(x) aov(x~category, data = df))
sapply(an, anova, simplify=FALSE)
lapply(an, function(m) HSD.test((m), "as.factor(category)", group = TRUE, console = TRUE))
lapply(an, function(m) TukeyHSD(aov(m)))




List <- names(df2)[1:3] # select just the variables
model1 <- lapply(List, function(x) {
  lm(substitute(i~category, list(i = as.name(x))), data = df2)})
lapply(model1, summary)
letters = lapply(model1, function(m) HSD.test((m), "category", group = TRUE, console = TRUE))

tx <- with(df2, interaction(age2, BMI))  # determining the factors
model2 <- lapply(List, function(x) {
  glm(substitute(i~tx, list(i = as.name(x))), data = df2)}) # using the factors already in "tx"

lapply(model2, summary)

letters = lapply(model2, function(m) HSD.test((m), "tx", alpha = 0.05, group = TRUE, console = TRUE))




##############################################################
N = 100
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
tic = proc.time()
df = createEmptyDf( N, 2, colnames = c( "id", "value" ) )
for( i in 1:N ){
  # df[ i, ] = c( i, i * 2 )
  df[ i, 1 ] = i
  df[ i, 2 ] = i * 2
}
print( proc.time() - tic )


createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
tic = proc.time()
df = createEmptyDf( N, 3, colnames = c( "id", "value" ,"value2") )
for( i in 1:N ){
  df[ i, ] = c( i, i*2, i**2)
}
print( proc.time() - tic )
#以下同じ（時間かかる）
#tic = proc.time()
#df = data.frame()
#for( i in 1:N ){
#  df = rbind( df, c( i, i * 2 ) )}
#print( proc.time() - tic )

##############################################################
#https://www.researchgate.net/post/How_can_I_create_a_file_based_on_the_output_of_a_certain_object_in_R_program
# create a list of TukeyHSD test object from dummy data:
data.test1 <- data.frame(a=rnorm(10), b=factor(round(runif(10))), c=factor(round(runif(10))))
data.test2 <- data.frame(d=rnorm(10), e=factor(round(runif(10))), f=factor(round(runif(10))))
tukey_full <- list(TukeyHSD(aov(lm(a~b*c, data.test1))), TukeyHSD(aov(lm(d~e*f, data.test2))))
#object to store the overall table:
table_full <- NULL
# loop all elements in the list
for(i in 1:length(tukey_full)){
  tukey <- tukey_full[[i]]
  factor_table <- unlist(lapply(tukey, function(x) nrow(x)))
  factor_table <- rep(names(factor_table), factor_table)
  tukey_bound <- NULL
  for (j in 1:length(tukey))
  {
    tukey_bound <- rbind(tukey_bound, tukey[[j]])
  }
  pairs <- rownames(tukey_bound)
  rownames(tukey_bound) <- NULL
  tukey_bound <- as.data.frame(tukey_bound)
  # add a column indicating the parameter that is compaired
  tukey_bound$parameter <- factor_table
  # add a column with the pairs of levels that are compaired
  tukey_bound$pairs <- pairs
  # add column with the model number, which corresponds to the index of the list of TukeyHSD test objects
  tukey_bound$model <- i
  table_full <- rbind(table_full, tukey_bound)
  }
write.csv(table_full, "filename.csv")
##############################################################
