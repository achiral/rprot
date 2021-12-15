




################################################################################
## 10 RStudio #####
################################################################################
1 + 1            # Ctr + Enter -> run
test = 1 + 1     # define object
################################################################################
## 12. package #####
################################################################################
.libPaths()
install.packages("GGally")
install.packages("patchwork")

require(GGally)
require(patchwork)
################################################################################
## 13 .Rprofile #####
################################################################################
require(stats)
require(tidyverse)
require(GGally)
require(patchwork)
require(lubridate)
################################################################################
## 14 Global option #####
################################################################################
# uncheck
# Tools > General > Basic Tab (RStudio > Preferences > General > Basic Tab)
# Workspace
# ( ) Restore .RData into worckspace at startup
# Save workspace to .RData on exit: Never
################################################################################
## 20 numerical data: integer, double #####
################################################################################
typeof(1L)      # integr
typeof(1)       # double
typeof(3.14)    # double
################################################################################
## 21 character data: character #####
################################################################################
typeof("apple") # character
typeof("10")    # character
typeof(10)      # double 
################################################################################
## 22 factor data: factor #####
################################################################################
vec_season_chr = 
 c("spring", "summer", "autumn", "winter")
vec_season_fct = 
 factor(vec_season_chr)

# change levels' order
df =
 tibble(
  season = vec_season_chr,
  avg_temp =c(22, 30, 18, 9)
 )
df %>%  # Ctrl + Shift + M
 ggplot(aes(x = season, y = avg_temp)) +
 geom_bar(stat = "identity")

vec_season_fct_rev =
 vec_season_fct %>% 
 fct_relevel("spring", "summer", "autumn", "winter")

df_rev =
 tibble(
  season = vec_season_fct_rev,
  avg_temp =c(22, 30, 18, 9)
 )

df_rev %>%
 ggplot(aes(x = season, y = avg_temp)) +
 geom_bar(stat = "identity")
################################################################################
## 23 order data: order #####
################################################################################
vec_season_chr =
 c("spring", "summer", "autumn", "winter")

vec_season_ord =
 factor(vec_season_chr, ordered = TRUE)     # ordered option

vec_season_fct =
 factor(vec_season_chr)

df_order =
 tibble(
  season =
   vec_season_ord %>% 
   fct_relevel("spring", "summer", "autumn", "winter"),
  avg_tmp =
   c(22, 20, 18, 9)
 )

df_factor =
 tibble(
  season =
   vec_season_fct %>% 
   fct_relevel("spring", "summer", "autumn", "winter"),
  avg_tmp =
   c(22, 20, 18, 9)
 )

df_order %>%
 filter(
  season >= "autumn"
 )

df_factor %>%
 filter(
  season >= "autumn"
 )
################################################################################
## 24 logical data: logical #####
################################################################################
typeof(TRUE)
typeof(FALSE)
################################################################################
## 26 vector #####
################################################################################
c(1, 5, 8, 13, 4)
c("a", "b", "CD", "EFg", "hiJK")
c(TRUE, TRUE, FALSE, TRUE, FALSE)
################################################################################
## 27 vector calculation #####
################################################################################
int_a = c(1, 5, 8, 13, 4)
int_b = c(2, 3, 5, 10, 98)
int_a + int_b

chr_a = c("a", "b", "c", "d", "e")
chr_b = c("1", "2", "3", "4", "5")
chr_a + chr_b

lgl_a = c(TRUE, TRUE, FALSE, FALSE, TRUE)
lgl_b = c(FALSE, FALSE, TRUE, FALSE, TRUE)
lgl_a + lgl_b

lgl_a_rev = c(1, 1, 0, 0, 1)
lgl_b_rev = c(0, 0, 1, 0, 1)
lgl_a_rev + lgl_b_rev
lgl_a_rev + lgl_b
################################################################################
## 28 recycle #####
################################################################################
c(1, 5, 8, 13) + c(2, 3)     # recycle
c(1, 5, 8, 13) + c(2, 3, 9)  # error
################################################################################
## 29 extract elements from vector #####
################################################################################
vec = c(1, 5, 8, 13, 4)
vec[3]                   # 3rd
vec[c(2, 4)]             # 2nd and 4th
vec[c(-1, -3)]           # not 1st, not 3rd

vec =
 c("one" =1,"two"= 5, "three" = 8, "four" = 13, "five" = 4)  # name
vec[c(1,3)]                                                  # 1st and 3rd
vec[c("one", "three")]                                       # name

vec = c(1, 5, 8, 13, 4)
vec[c(TRUE, FALSE, TRUE, TRUE, FALSE)]  # vector : TRUE
vec >= 5                                # logical: TRUE or FALSE
vec[vec >= 5]                           # vector : TRUE
################################################################################
## 30 Rbase function #####
################################################################################
vec = c(1, 5, 8, 13, 4)
sum(vec)
mean(vec)
median(vec)
max(vec)
min(vec)
quantile(vec)
var(vec)
sd(vec)

lgl_vec = c(TRUE, TRUE, FALSE, FALSE, TRUE)
sum(lgl_vec)   # count TRUE
mean(lgl_vec)  # rate of TRUE
################################################################################
## 31 make vector #####
################################################################################
# 1 to 10 same interval
1:10
# 1.5 to 10.1 range 8 same interval
seq(1.5, 10.1, length = 8)
# 1.5 to 10.1 interval 2.5
seq(1.5, 10.1, by = 2.5)
# c("red", "blue", "black") repeat three times
rep(c("red", "blue", "black"), times = 3)
# c("red", "blue", "black") repeat length 10
rep(c("red", "blue", "black"), length = 10)
# c("red", "blue", "black") repeat five times each element
rep(c("red", "blue", "black"), each = 5)
# 1 to 5 repeat five times each element
rep(1:5, each = 5)
################################################################################
## 34 make logical vector #####
################################################################################
# x: 1-10 define vector
x = 1:10
# y: c(1, 3, 4, 7, 2) define vector
y = c(1, 3, 4, 7, 2)
# Is x equal to y?
x == y
# Is x not equal to y?
x != y
# Is x same or more than y?
x >= y
# Is x less than y?
x < y
# Is x included in y?
x %in% y
################################################################################
## 35 logical vector calculation #####
################################################################################
# cabbage size vector
size_vec = 
 c(7, 13, 5, 25, 38, 57, 18, 9, 32, 28)
# same or more than 10?
lgl_vec_1 =
 size_vec >= 10
# same or less than 30?
lgl_vec_2 =
 size_vec <= 30
# AND
lgl_vec_1 & lgl_vec_2
# OR
lgl_vec_1 | lgl_vec_2
# XOR
xor(lgl_vec_1, lgl_vec_2)
# NOT
!lgl_vec_1
!lgl_vec_2
# ANY: is there TRUE?
any(lgl_vec_1)
any(lgl_vec_2)
# ALL: is it all TRUE?
all(lgl_vec_1)
all(lgl_vec_2)
################################################################################
## 37 list #####
################################################################################
list_sample =
 list(
  c(1, 3, 4, 7, 2),        # int
  c("d", "a", "t", "a"),   # chr
  8,                       # int
  list(                    # list
   c(2, 0, 2, 0),          # int
   c("c", "a", "t")        # chr
  ),
  c(TRUE, FALSE, TRUE)     # logi
 )
################################################################################
## 38 extract elements from list #####
################################################################################
list_sample[[2]]       # vector (without list)
list_sample[2]         # list
list_sample[[4]][[1]]  # vector (without list)
list_sample[4]         # list
################################################################################
## 39 named list #####
################################################################################
list_sample_name =
 list(
  "one" = c(1, 3, 4, 7, 2),
  "two" = c("d", "a", "t", "a"),
  "three" = 8,
  "four" = list(
   c(2, 0, 2, 0),
   c("c", "a", "t")
  ),
  "five" = c(TRUE, FALSE, TRUE)
 )

list_sample_name[[1]]      # 1st
list_sample_name[["one"]]  # name == "one"
list_sample_name$one       # name == "one"
################################################################################
## 41 dataframe #####
################################################################################
df =
 data.frame(
  season = 
   c("spring", "summer", "autumn", "winter"),
  avg_temp =
   c(22, 30, 18, 9)
 )
################################################################################
## 42 extract elements from dataframe #####
################################################################################
df[[1]]         # 1st : vector
df[["season"]]  # name: vector
df$season       # name: vector
################################################################################
## 43 tibble #####
################################################################################
tibble(
 season = 
  c("spring", "summer", "autumn", "winter"),
 avg_temp =
  c(22, 30, 18, 9)
)
################################################################################
## 45 define function #####
################################################################################
add_1 = function(x){
 y = x +1
 return(y)
}
add_1(10)    # 10 + 1
add_1(1:10)  # 1:10 +1
################################################################################
## 47 modify function #####
################################################################################
head(iris)
iris_setosa =
 iris %>% 
 tibble() %>% 
 filter(Species == "setosa")

iris_versicolor =
 iris %>% 
 tibble() %>% 
 filter(Species == "versicolor")

iris_virginica =
 iris %>% 
 tibble() %>% 
 filter(Species == "virginica")

iris_setosa %>% 
 ggplot(
  mapping = aes(x = Sepal.Length, y = Petal.Length)
 ) +
 geom_point() +
 geom_smooth(method = "lm", se = FALSE) +
 theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25)
 )

iris_versicolor %>% 
 ggplot(
  mapping = aes(x = Sepal.Length, y = Petal.Length)
 ) +
 geom_point() +
 geom_smooth(method = "lm", se = FALSE) +
 theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25)
 )

iris_virginica %>% 
 ggplot(
  mapping = aes(x = Sepal.Length, y = Petal.Length)
 ) +
 geom_point() +
 geom_smooth(method = "lm", se = FALSE) +
 theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25)
 )

## define funcrion
make_scatter = function(df, type){
 df_Species =
  df %>% 
  tibble() %>% 
  filter(Species == type)
 
 scatter =
  df_Species %>% 
  ggplot(
   mapping = aes(x = Sepal.Length, y = Petal.Length)
  ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme(
   axis.text = element_text(size = 25),
   axis.title = element_text(size = 25)
  )
 
 return(scatter)
}

make_scatter(df = iris, type = "setosa")
make_scatter(df = iris, type = "versicolor")
make_scatter(df = iris, type = "virginica")
################################################################################
## 48 make external file of function #####
################################################################################
## define function of "make_scatter()" into my_functions.R
## add "source(my_functions.R, encoding = "utf8")" in ".Rprofile"
################################################################################
## 52 pipe %>%  Ctrl + Shift + M #####
################################################################################
head(txhousing)
## not use %>% 
dplyr::filter(dplyr::select(txhousing, city, year, month, sales),
       sales >= 100)
## use %>% 
txhousing %>% 
 dplyr::select(city, year, month, sales) %>%
 dplyr::filter(sales >= 100)
################################################################################
## 53 pipe  %>%  Ctrl + Shift + M #####
################################################################################
## not use %>% 
filter(
 select(txhousing, city, year, month, sales),
 sales >= 100
)
## use %>% 
filter(
 txhousing %>% select(city, year, month, sales),
 sales >= 100
)
## use %>% 
txhousing %>% select(city, year, month, sales) %>% 
 filter(
  sales >= 100
  )
################################################################################
## 56 make tibble #####
################################################################################
head(iris)
iris %>% as_tibble()
as_tibble(iris)
tibble(                                # column = row(vector)
 No. = 1:3,                            # int
 laguage = c("R", "Python", "Julia")   # chr
)
## markdown
tribble(
 ~No., ~language,
 1, "R",
 2, "Python",
 3, "Julia"
)
################################################################################
## 57 usefulness of tibble #####
################################################################################
df = tibble(int = c(5,10,76),
            dbl = c(3.124, 8.425, 5.32),
            chr = c("a", "c", "E"))

iris %>% as_tibble() %>% distinct(Species)                          # unique
data_1 = iris %>% as_tibble() %>% filter(Species == "setosa")
data_2 = iris %>% as_tibble() %>% filter(Species == "versicolor")
data_3 = iris %>% as_tibble() %>% filter(Species == "virginica")

df_rev_1 = df %>% mutate(data = list(data_1, data_2, data_3))
df_rev_1 %>% pluck("data", 1)                                       # access data_1

fig_1 = make_scatter(iris, "setosa")
fig_2 = make_scatter(iris, "versicolor")
fig_3 = make_scatter(iris, "virginica")

df_rev_2 = df_rev_1 %>% mutate(fig = list(fig_1, fig_2, fig_3))
df_rev_2 %>% pluck("fig", 2)                                        # access fig_2
df_rev_2 %>% pluck("fig", 3)                                        # access fig_3
################################################################################
## 59 dplyr::filter #####
################################################################################
diamonds
diamonds %>% dplyr::filter(carat <= 0.25)
################################################################################
## 60 dplyr::filter #####
################################################################################
df = tibble(No = 1:7,
            gender = c("M", "F", "F", "M", "F", "M", "M"),
            height = c(165, 150, 170, 175, 165, 195, 180))
df %>% dplyr::filter(gender == "F")
df %>% dplyr::filter(df$gender == "F")
df %>% dplyr::filter(c(F, T, T, F, T, F, F))
################################################################################
## 60 dplyr::filter #####
################################################################################
diamonds
# price == 377
diamonds %>% dplyr::filter(price == 377)
# depth >= 62
diamonds %>% dplyr::filter(depth >= 62)
# carat from 0.23 to 0.27
diamonds %>% dplyr::filter(between(carat, 0.23, 0.27))
diamonds %>% dplyr::filter(carat >= 0.23 & carat <= 0.27)
# cut >= Good (levels)
diamonds$cut %>% levels()
diamonds %>% dplyr::filter(cut >= "Good")
# color not "E"
diamonds %>% dplyr::filter(color != "E")
# color equals "I" or "J"
diamonds %>% dplyr::filter(color %in% c("I","J"))
diamonds %>% dplyr::filter(color == "I" | color == "J")
# cut equals "NA"
diamonds %>% dplyr::filter(is.na(cut))
# make dataframe including "NA"
df = tibble(col = c(1, 3, NA, 5, 4, 8, NA, 9, NA))
# col equals "NA"
df %>% filter(is.na(col))
# col not equal "NA"
df %>% filter(!is.na(col))
# depth >= 62 and color == "H"
diamonds %>% dplyr::filter(depth >= 62 & color == "H")
diamonds %>% dplyr::filter(depth >= 62, color == "H")
diamonds %>% dplyr::filter(depth >= 62) %>% filter(color == "H")
# depth >= 62 or color == "H"
diamonds %>% dplyr::filter(depth >= 62 | color == "H")
# depth >= 62 xor color == "H"
diamonds %>% dplyr::filter(xor(depth >= 62, color == "H"))
# depth !>= 62 and color != "H"
diamonds %>% dplyr::filter(!(depth >= 62 & color == "H"))
################################################################################
## 62 dplyr::select #####
################################################################################
iris = iris %>% as_tibble()
iris %>% dplyr::select(Sepal.Length, Species)
################################################################################
## 63 dplyr::select #####
################################################################################
iris = iris %>% as_tibble()
iris %>% dplyr::select(Sepal.Length, Species)
iris %>% dplyr::select(1, 5)
################################################################################
## 64 dplyr::select #####
################################################################################
# directly select
# Sepal.Length, Species
iris = iris %>% as_tibble()
iris %>% dplyr::select(Sepal.Length, Species)
# range
# from Sepal.Length to Petal.Width
iris %>% select(Sepal.Length:Petal.Width)
# indirectly select
# without Sepal.Length, Species
iris %>% dplyr::select(-Sepal.Length, -Petal.Width)
# starts_with
# Sepal
iris %>% dplyr::select(starts_with("Sepal"))
# ends_with
# Width
iris %>% dplyr::select(ends_with("Width"))
# contains
# l.L
iris %>% dplyr::select(contains("l.L"))
# matches
# start with "Se", contain "en"
iris %>% dplyr::select(matches("^Se.*en.*"))
# all of
# Sepal.Length to Petal.Length
var = c("Sepal.Length", "Petal.Length")
iris %>% dplyr::select(all_of(var))
iris %>% dplyr::select(var)
# as character
iris %>% dplyr::select(c("Sepal.Length", "Petal.Length"))
# everything
iris %>% select(everything())
# order of columns: Species, everything
iris %>% select(Species, everything())
# where: type of data
# dbl
iris %>% select(where(is.double))
iris %>% select(where(is.factor))
iris %>% select(where(is.character))
################################################################################
## 65 dplyr::mutate #####
################################################################################
df = tibble(No. = 1:5,
            chr = c("one", "two", "three", "four", "five"))
df %>% mutate(
 a = 10:14,
 b = 15:19,
 c = 20:24,
 d = 11111:11115
)
################################################################################
## 66 dplyr::mutate #####
################################################################################




