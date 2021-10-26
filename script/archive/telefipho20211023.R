## TeleFipho
## input multiple files(https://yyhhyy.hatenablog.com/entry/2017/04/08/120000)
################################################################################
## setup #####
################################################################################
rm(list = ls(all = TRUE))
setwd("~/Dropbox/GitHub/local/Docker/R/rprot")
# source("/home/rstudio/rproject/script/archive/functions.R")        # load functions
# source("/home/rstudio/rproject/script/archive/functions_DEPpkg.R") # load DEPpkg functions
source("script/startup.R")                  # load packages
source("script/archive/functions.R")
detach_all()
source("script/startup.R")                  # load packages
source("script/archive/functions.R")
################################################################################
## set working directory #####
################################################################################
getwd()
# setwd("/home/rstudio/rproject")
# setwd("/home/rstudio/rproject/data/Perseus_Like_Analysis/PFC")
# setwd("~/Dropbox/GitHub/local/Docker/R/rprot")

setwd("/Users/user/Dropbox/GitHub/local/Docker/R/rprot/data/telefipho")
DIR_HOST <- getwd()  # host working directory
setwd(DIR_HOST)
DATA_FOLDER <- 20210925
path <- paste(DIR_HOST, "/" ,DATA_FOLDER, sep="")
## check file name
dir(paste("./", DATA_FOLDER, sep =""))
## input data set
file_list <- list.files(path, full.names = TRUE)
df <- read_files(file_list)
## extract shock data
df1 <- shock_data(df = df, file_no = 1, shock_time = 227.29)
df2 <- shock_data(df = df, file_no = 2, shock_time = 160.26)
df3 <- shock_data(df = df, file_no = 3, shock_time = 346.46)
df4 <- shock_data(df = df, file_no = 4, shock_time = 511.30)
df5 <- shock_data(df = df, file_no = 5, shock_time = 511.30)
df6 <- shock_data(df = df, file_no = 6, shock_time = 227.29)
df7 <- shock_data(df = df, file_no = 7, shock_time = 160.26)
df8 <- shock_data(df = df, file_no = 8, shock_time = 346.46)
df9 <- shock_data(df = df, file_no = 9, shock_time = 511.30)
df10 <- shock_data(df = df, file_no = 10, shock_time = 511.30)
################################################################################
## summarize data frames
df_out <- data.frame(seq(-3, 10, by = 0.01))  # -3 から 10 まで、0.01ずつ増えるベクトル
out <- cbind(df_out, df1[2], df2[2], df3[2], df4[2])
# df_out2 <- cbind(df_out, df1[2], df2[2], df3[2], df4[2], df5[2], df6[2], df7[2], df8[2], df9[2], df10[2])

## rename columns
colnames(out) <- c("Time", "aCLZ01", "aCLZ02", "aPCP03", "aCLZ04")

## output xlsx
# library(openxlsx) #入出力(write.xlsx)
# smp <- list("out" = out)
# write.xlsx(smp, file = paste(DATA_FOLDER, ".xlsx", sep=""))  # overwrite = T
write.xlsx(out, file = paste(DATA_FOLDER, ".xlsx", sep=""))  # overwrite = T
################################################################################