#20210206ウェブ講義
#caret
#####################################################
install.packages("tidyverse")
install.packages("caret")
install.packages("doParallel")
install.packages("tictoc")
install.packages("mlbench")
install.packages("e1071")
install.packages("modelr")
install.packages("tiptoc")
#####################################################
library(tidyverse)
library(caret)
library(doParallel)
library(tictoc)
library(mlbench)
library(e1071)
library(modelr)
library(tiptoc)
#####################################################
Random forest

pred <- predict(model_rt, df_test)


moderl_rf %>%
  varImp() %>%
  ggplot(top = 4)
#####################################################
Karnel mrthod



#####################################################
Titanic_EDA
full <- as.data.frame(Titanic)
#bind_rows(train.test)
str(full) # Pandas Info

ggplot(data=full[1:LT],aes(x=Embarked, fill=Survived))+geom_bar(position=)
#####################################################
x <- col(df); x
install.packages("orrplot")
library(corrplot)
corrplot(x)

factanal(x = df, factors = 3, )



set.seed()

fitControl <- trainControl(method = "repeatedcv",
                           )

#Karnel method
install.packages("Karnelab")
library(Karnelab)

#boost method
install.packages("xgboost")
library(xgboost)

xgb.plot.multi.trees
#####################################################
library(data.table)
library()

options()



