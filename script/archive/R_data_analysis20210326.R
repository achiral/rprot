#【オンライン教室】Rデータ分析超入門202103261930-2100
#RStudio=統合開発環境
#Filesでは読み込みができる様になっている
#dplyr=データフレーム集計
#ggplot2=グラフ描写
#mlr=機械学習
#for文,if文
#Connectionでgithubと連結できる
################################################################################
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/untitled folder")
getwd()
################################################################################
print("hello world")
kion = c(10, 20, 16, 15, 18, 20, 14)
mean(kion)
max(kion)
plot(kion)
################################################################################
#install.packages("rtweet")
#install.packages("gapminder")
#install.packages("ggvis")
################################################################################
library(rtweet)
x = search_tweets("鬼滅の刃", n=100) #@kaikeitameshi
head(x$text)
################################################################################
library(ggplot2)
library(dplyr)
library(gapminder)
library(corrr)

kion = c(10, 20, 30)

gapminder_2007 <- gapminder %>% filter(year == 2007)
ggplot(gapminder_2007, aes(x = lifeExp, y = gdpPercap)) + geom_point()
gapminder_2007 %>%
  ggplot(aes(x = lifeExp, y = gdpPercap, 
             color = continent, #大陸
             size = pop #人口
             )) + 
  geom_point() #散布図(バブルチャート)
#Exportから保存できる
plot(gapminder_2007)

gapminder_2007 %>%
  group_by(continent) %>%
  summarize(mean(pop))
################################################################################
library(ggvis)
mtcars %>%
  ggvis(~wt) %>%
  layer_histograms(width = input_slider(0,2,step=0.10,label="width"),
                  center=input_slider(0,2,step=0.05,label="center"))
################################################################################








