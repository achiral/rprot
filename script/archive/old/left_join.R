library(readxl) #エクセル入力(read_excel)
library(tidyverse)
library(writexl) #xlsx出力
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/R")
dat1 <- read_excel("repo.xlsx", 3)
dat2 <- read_excel("repo.xlsx", 4)
dat3 <- left_join(dat1, dat2, by = "name")
write_xlsx(dat3, "repo2.xlsx", format_headers = FALSE)