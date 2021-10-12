#https://www.r-bloggers.com/2019/07/cpat-and-the-renyi-type-statistic-end-of-sample-change-point-detection-in-r/

#https://github.com/ntguardian/CPAT
install.packages("CPAT")
#devtools::install()
#devtools::install_github("ntguardian/CPAT")
library(CPAT)
set.seed(20180924)

(vec1 <- c(rnorm(100, 0), rnorm(100, 1)))
CUSUM.test(vec1)
DE.test(vec1)
HS.test(vec1)
HR.test(vec1)

vec2 <- as.numeric(arima.sim(200, 
                             model = list(
                               order = c(1, 0, 0),
                               ar = 0.4)))
CUSUM.test(vec2)


CUSUM.test(vec2, use_kernel_var = TRUE, kernel = "qs", bandwidth = "nw")


install.packages("strucchange")
install.packages("dynlm")
install.packages("cointReg")
library(strucchange)
library(dynlm)
library(cointReg)

data(USIncExp)
incexpres <- residuals(
  dynlm(d(log(expenditure)) ~ d(log(income)),
        data = USIncExp))
CUSUM.test(incexpres, use_kernel_var = TRUE)

HR.test(vec1)
HR.test(vec2, use_kernel_var = TRUE, kernel = "qs", bandwidth = "and")
HR.test(vec2, use_kernel_var = TRUE, kn = sqrt)


vec3 <- c(rnorm(5, mean = 0), rnorm(195, mean = 1))
CUSUM.test(vec3)
HR.test(vec3)

HR.test(incexpres, use_kernel_var = TRUE)



DE.test(vec1)
DE.test(vec2, use_kernel_var = TRUE, kernel = "qs", bandwidth = "nw")
DE.test(vec3)
DE.test(incexpres, use_kernel_var = TRUE)

HS.test(vec1)
HS.test(vec2, corr = TRUE)
HS.test(vec3)

Andrews.test(vec1, 100)
Andrews.test(vec2, 100)
Andrews.test(rev(vec3), 190)

mod <- dynlm(d(log(expenditure)) ~ d(log(income)), data = USIncExp)
X <- as.data.frame(model.frame(mod))
names(X) <- c("exp", "inc")
Andrews.test(exp ~ inc, x = X, M = 300)

