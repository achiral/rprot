#post hoc test
#https://biolab.sakura.ne.jp/anova-post-hoc-multiple-comparison.html
#### スクリプト開始
n1 <- n2 <- n3 <- 20
np <- 0
np.anov <- 0
anova.p <- as.numeric(NULL)
k <- 100000
group <- factor(rep(1:3, c(n1, n2, n3)))

for (i in 1:k) {  
  x1 <- rnorm(n1, mean = 0, sd = 1)
  x2 <- rnorm(n2, mean = 0, sd = 1)
  x3 <- rnorm(n3, mean = 0, sd = 1)
  res <- aov(c(x1, x2, x3) ~ group)
  anova.p[i] <- anova(res)$Pr[1]
  pt <- TukeyHSD(res)$group[, 4]
  if (pt[1] < 0.05 || pt[2] < 0.05 || pt[3] < 0.05) {
    np <- np + 1
  }
  if (anova.p[i] < 0.05 && (pt[1] < 0.05 || pt[2] < 0.05 || pt[3] < 0.05)) {
    np.anov <- np.anov + 1
  }
}

hist(anova.p)
# 5% 水準での棄却率
# Tukey 検定のみ
rnp <- np/k
# ANOVA 有意後，Tukey 検定
rnp.anov <- np.anov/k

barplot(
  c(rnp, rnp.anov), ylim = c(0, 0.06), col = "#0080ff",
  names.arg =c("Tukey test only", "Tukey test after ANOVA"),
  ylab = "Type I error rates with alpha level 0.05"
)
abline(h = 0.05, col = "red", lwd = 2)
### スクリプト終了