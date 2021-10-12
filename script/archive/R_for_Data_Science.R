#『Rで始めるデータサイエンス』オライリー・ジャパン、ISBN978-4-87311-814-7．
#https://github.com/hadley/r4ds
################################################################################
library(tidyverse)
library(maps)
library(mapproj) #coord_map
library(nycflights13)
################################################################################
#まえがき####
devtools::install_github("hadley/r4ds")
install.packages("tidyverse")
install.packages(c("nycflights13", "gapminder", "Lahman"))
library(tidyverse)
tidyverse_update()
################################################################################
#1章ggplot2によるデータ可視化####
library(tidyverse)
ggplot2::mpg
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy))

#練習問題####
#1
ggplot(data = mpg)
#2
str(mtcars)
#3
?mpg #the type of drive train, where f = front-wheel drive, r = rear wheel drive, 4 = 4wd
#4
ggplot(data = mpg) +
  geom_point(mapping = aes(x = hwy, y = cyl))
#5
ggplot(data = mpg) +
  geom_point(mapping = aes(x = class, y = drv)) #カテゴリ変数のため
#1.3エステティックマッピング####
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, color = class))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, size = class))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, alpha = class))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, shape = class))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy), color = "blue")
#練習問題####
#1
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy, color = "blue"))
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy), color = "blue")
#2
str(mpg)
?mpg
#3
ggplot(data = mpg) +
  geom_point(aes(x = displ, y = hwy, color = manufacturer, size = model, shape = trans))
ggplot(data = mpg) +
  geom_point(aes(x = displ, y = hwy, color = year, size = cyl, shape = cty))
ggplot(data = mpg) +
  geom_point(aes(x = displ, y = hwy, color = year, size = cyl, shape = fl))
#4
ggplot(data = mpg) +
  geom_point(aes(x = displ, y = hwy, color = drv, size = drv, shape = drv))
#5
ggplot(data = mpg) +
  geom_point(aes(x = displ, y = hwy, stroke = 0.01))
?geom_point
#6
ggplot(data = mpg) +
  geom_point(aes(x = displ, y = hwy, color = displ < 5))
#1.5ファセット
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_wrap(~ class, nrow = 2)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_wrap(~ class, nrow = 2)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(drv ~ cyl)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(. ~ cyl)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(drv ~ .)
#練習問題####
#1
str(mpg)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_wrap(~ displ)
#2
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(drv ~ cyl)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = drv, y = cyl))
#3
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(drv ~ .)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(. ~ cyl)
#4
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_wrap( ~ class, nrow = 2)
#5
?facet_wrap
?facet_grid
#6
#視認性の問題
#1.6幾何オブジェクト####
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy))
ggplot(data = mpg) +
  geom_smooth(mapping = aes(x = displ, y = hwy))
ggplot(data = mpg) +
  geom_smooth(mapping = aes(x = displ, y = hwy, linetype = drv))
ggplot(data = mpg) +
  geom_smooth(mapping = aes(x = displ, y = hwy, group = drv))
ggplot(data = mpg) +
  geom_smooth(mapping = aes(x = displ, y = hwy, color = drv), show.legend = FALSE)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy)) +
  geom_smooth(mapping = aes(x = displ, y = hwy))
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point() +
  geom_smooth()
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point(mapping = aes(color = class)) +
  geom_smooth()
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point(mapping = aes(color = class)) +
  geom_smooth(data = filter(mpg, class == "subcompact"), se = FALSE) +
  geom_step(mapping = aes(color = class))
#練習問題####
#1
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_line(mapping = aes(color = class))
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) +
  geom_boxplot(mapping = aes(color = class))
ggplot(data = mpg, mapping = aes(displ)) +
  geom_bar()
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_area(mapping = aes(color = class))
#2
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point() +
  geom_smooth(se = FALSE)
#3
ggplot(data = mpg) +
  geom_smooth(mapping = aes(x = displ, y = hwy, color = drv), show.legend = FALSE)
#4
#標準誤差
#5
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point() +
  geom_smooth()
ggplot() +
  geom_point(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_smooth(data = mpg, mapping = aes(x = displ, y = hwy))
#6
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point() +
  geom_smooth(se = FALSE)
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point() +
  geom_smooth(mapping = aes(group = drv), se = FALSE)
ggplot(data = mpg, mapping = aes(x = displ, y = hwy, color = drv)) +
  geom_point() +
  geom_smooth(mapping = aes(group = drv), se = FALSE)
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point(mapping = aes(color = drv)) +
  geom_smooth(se = FALSE)
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point(mapping = aes(color = drv)) +
  geom_smooth(mapping = aes(linetype = drv), se = FALSE)
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_point(shape = 21, color = "white", aes(fill = drv), size = 3, stroke = 2)
#1.7統計変換####
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut))
ggplot(data = diamonds) +
  stat_count(mapping = aes(x = cut))
demo <- tribble(~a, ~b,
                "bar_1", 20,
                "bar_2", 30,
                "bar_3", 40)
ggplot(data = demo) +
  geom_bar(mapping = aes(x = a, y = b), stat = "identity")
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, y = ..prop.., group = 1))
ggplot(data = diamonds) +
  stat_summary(mapping = aes(x = cut, y = depth),
                             fun.ymin = min,
                             fun.ymax = max,
                             fun.y = median)
#練習問題####
#1
ggplot(data = diamonds) +
  stat_summary(mapping = aes(x = cut, y = depth))
#2
ggplot(data = diamonds) +
  geom_col(mapping = aes(x = cut, y = depth))

#5
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, y = ..prop..))
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, y = ..prop.., group = 1))
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = color, y = ..prop.., group = 1))
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, group = color, fill = color, y = ..prop.., group = 1))

#1.8位置調整####
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, color = cut))
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = cut))
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = clarity))

#identity 積み上げ
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "identity")
ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) +
  geom_bar(alpha = 1/5, position = "identity")
ggplot(data = diamonds, mapping = aes(x = cut, color = clarity)) +
  geom_bar(fill = NA, position = "identity")

#fill 比率
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "fill")

#dodge 隣合わせ
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "dodge")

#jitter 乱れ
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy), position = "jitter")
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) +
  geom_jitter()

#練習問題####
#1
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) +
  geom_point()
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) +
  geom_jitter() +
  geom_count()
#2
?geom_jitter
#3
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) +
  geom_jitter()
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) +
  geom_count()
#4
?geom_boxplot #dodge2

#1.9座標系
#coord_flip()x軸,y軸交換
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) +
  geom_boxplot()
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) +
  geom_boxplot() +
  coord_flip()

#coord_quickmap()縦横比
install.packages("maps")
library(maps)

nz <- map_data("nz")

ggplot(nz, aes(long, lat, group = group)) +
  geom_polygon(fill = "white", color = "black")

ggplot(nz, aes(long, lat, group = group)) +
  geom_polygon(fill = "white", color = "black") +
  coord_quickmap()

#coord_polar()極座標
bar <- ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = cut), show.legend = FALSE, width =1) +
  theme(aspect.ratio = 1) +
  labs(x = NULL, y = NULL)
bar
bar + coord_flip()
bar + coord_polar()

#練習問題####
#1
bar_ex1 <- ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, fill = clarity), show.legend = FALSE, width =1, position = "identity") +
  theme(aspect.ratio = 1) +
  labs(x = NULL, y = NULL)
bar_ex1 + coord_polar()
#2
?labs()
#3
install.packages("mapproj")
library(mapproj) #coord_map

nz <- map_data("nz")
ggplot(nz, aes(long, lat, group = group)) +
  geom_polygon(fill = "white", color = "black") +
  coord_quickmap()
ggplot(nz, aes(long, lat, group = group)) +
  geom_polygon(fill = "white", color = "black") +
  coord_map()
#4
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) +
  geom_point() +
  geom_abline() +
  coord_fixed()
?geom_abline()
?coord_fixed()

#1.10階層グラフィックス文法
#ggplot(data = <DATA>) +
#  <GEOM_FUNCTION>(mapping = aes(<MAPPINGS>),
#                  stat = <STAT>,
#                  position = <POSITION>) +
#  <COORDINATE_FUNCTION> +
#  <FACET_FUNCTION>
################################################################################
#2章ワークフロー：基礎
#2.1コーディングの基本
#"<-" : "option" + "-"
#2.3関数呼び出し
seq(1,10)
x <- "hello_world"
y <- seq(1, 10, length.out =5)
y
(y <- seq(1, 10, length.out =5))
?seq()
#練習問題####
#1
my_variable <- 10
my_varlable
my_variable
#2
library(tidyverse)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy))
filter(mpg, cyl == 8)
filter(diamonds, carat > 3)
#3
#Alt-Shift-k #shortcut一覧表示
################################################################################
#3章dplyrによるデータ変換
#3.1.1準備するもの####
library(nycflights13)
library(tidyverse)
#3.1.2nycflights13
?flights
flights
View(flights)

#3.3.2filter()で行にフィルタをかける
jan1 <- filter(flights, month == 1, day == 1)
(dec25 <- filter(flights, month ==12, day == 25))

#3.2.1比較
sqrt(2) ^ 2 == 2
near(sqrt(2) ^ 2, 2)
1/49 * 49 == 1
near(1/49 * 49, 1)

#3.2.2論理演算子
filter(flights, month == 11 | month ==12)
nov_dec <- filter(flights, month %in% c(11,12))
filter(flights, !(arr_delay > 120 | dep_delay > 120))
filter(flights, arr_delay <= 120, dep_delay <= 120)

#3.2.3欠損値
df <- tibble(x = c(1, NA, 3))
filter(df, x > 1)
filter(df, is.na(x) | x > 1)

#3.3arrange()で行を配置する
arrange(flights, year, month, day)
arrange(flights, desc(arr_delay))

#3.4select()で列を選ぶ
select(flights, year, month, day)
select(flights, year:day)
select(flights, -(year:day))

#ヘルパー関数
start_with()
ends_with()
contains()
matches()num_range()

#rename
colnames(flights)
rename(flights, tail_num = tailnum)
select(flights, time_hour, air_time, everything())

#3.5mutate()で新しい変数を追加する
flights_sml <- select(flights, year:day, ends_with("delay"), distance, air_time)
mutate(flights_sml, 
       gain = arr_delay - dep_delay,
       speed = distance / air_time * 60,
       hours = air_time / 60,
       gain_per_hour = gain / hours)
transmute(flights_sml, 
          gain = arr_delay - dep_delay,
          speed = distance / air_time * 60,
          hours = air_time / 60,
          gain_per_hour = gain / hours)

#3.5.1有用な作成関数
transmute(flights, dep_time,
          hour = dep_time %/% 100,
          minute = dep_time %% 100)

#オフセット
x <- 1:10
lead(x)
lag(x)

#累積および回転和
x
cumsum(x)
cummean(x)
cumprod(x)
cummin(x)
cummax(x)

#ランク付け
y <- c(1, 2, 2, NA, 3, 4)
min_rank(y)
min_rank(desc(y))
row_number(y)
dense_rank(y)
percent_rank(y)
cume_dist(y)
ntile(y, 5)

#3.6summarize()によるグループごとの要約
summarize(flights, delay = mean(dep_delay, na.rm = TRUE))
by_day <- group_by(flights, year, month, day)
summarize(by_day, delay = mean(dep_delay, na.rm = TRUE))

#3.6.1パイプで複数演算を結合する
by_dest <- group_by(flights, dest)
delay <- summarize(by_dest,
                   count = n(),
                   dist = mean(distance, na.rm = TRUE),
                   delay = mean(arr_delay, na.rm = TRUE))
delay <- filter(delay, count > 20, dest != "HNL")
ggplot(data = delay, mapping = aes(x = dist, y = delay)) +
  geom_point(aes(size = count), alpha = 1/3) +
  geom_smooth(se = FALSE)

delay <- flights %>% group_by(dest) %>%
  summarize(count = n(),
            dist = mean(distance, na.rm = TRUE),
            delay = mean(arr_delay, na.rm = TRUE)) %>%
  filter(count > 20, dest != "HNL")

#3.6.2欠損値
flights %>% group_by(year, month, day) %>%
  summarize(mean = mean(dep_delay))
flights %>% group_by(year, month, day) %>%
  summarize(mean = mean(dep_delay, na.rm = TRUE))

not_cancelled <- flights %>% filter(!is.na(dep_delay), !is.na(arr_delay))
not_cancelled %>% group_by(year, month, day) %>%
  summarize(mean = mean(dep_delay))

#3.6.3カウント
delays <- not_cancelled %>% group_by(tailnum) %>%
  summarize(delay = mean(arr_delay))
ggplot(data = delays, mapping = aes(x = delay)) +
  geom_freqpoly(binwidth = 10)
  
delays <- not_cancelled %>% group_by(tailnum) %>%
  summarize(delay = mean(arr_delay, na.rm = TRUE),
            n = n())
ggplot(data = delays, mapping = aes(x = n, y = delay)) +
  geom_point(alpha = 1/10)

delays %>% filter(n > 25) %>%
  ggplot(mapping = aes(x = n, y = delay)) +
  geom_point(alpha = 1/10)

batting <- as_tibble(Lahman::Batting)
batters <- batting %>% group_by(playerID) %>%
  summarize(ba = sum(H, na.rm = TRUE) / sum(AB, na.rm = TRUE),
            ab = sum(AB, na.rm = TRUE))
batters %>%
  filter(ab > 100) %>%
  ggplot(mapping = aes(x = ab, y = ba)) +
   geom_point() +
   geom_smooth(se = FALSE)
batters %>% arrange(desc(ba))

#3.6.4便利な要約関数
#中心傾向の代表値：mean(x), median(x)
not_cancelled %>%
  group_by(year, month, day) %>%
  summarize(avg_delay1 = mean(arr_delay),
            avg_delay2 = mean(arr_delay[arr_delay > 0]))
not_cancelled$arr_delay[not_cancelled$arr_delay > 0]

#散らばり（広がり）の代表値：平均二乗偏差・標準偏差sd(x),四分位範囲IQR(x),中央値絶対偏差mad(x)
not_cancelled %>%
  group_by(dest) %>%
  summarize(distance_sd = sd(distance)) %>%
  arrange(desc(distance_sd))

#ランクの代表値：min(x),quantile(x, 0.25),max(x)
not_cancelled %>%
  group_by(year, month, day) %>%
  summarize(first = min(dep_time),
            last = max(dep_time))

#位置の代表値：first(x) = x[1],nth(x, 2) = x[2],last(x) = x[length(x)]
not_cancelled %>%
  group_by(year, month, day) %>%
  summarize(first_dep = first(dep_time),
            last_dep = last(dep_time))

not_cancelled %>%
  group_by(year, month, day) %>%
  mutate(r = min_rank(desc(dep_time))) %>%
  filter(r %in% range(r))
  
#カウント：n(),n_distinct(x)





################################################################################


################################################################################
