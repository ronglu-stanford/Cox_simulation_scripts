rm(list=ls())
set.seed(3859)
## -----------------------------------------------------------------------------------------------
## instantaneous risk: h(t) = h_0(t) * exp(b*x)
## cumulative risk function: H(t) = H_0(t) * exp(b*x) = integral of h(u) over 0 to t
## survival function: S(t) = exp{-H(t)} = S_0(t) ^ {exp(b*x)} = 1 - F(t) = Pr[T > t]
## probability density function: f(t) = S(t)*h(t) = S_0(t) ^ {exp(b*x)} * h_0(t) * exp(b*x)
## cummulative distribution function (CDF): F(t)

## ---------------- Cox Model:
## S(t) = exp{ - H_0(t) * exp(b*x) } where H_0(t) = lambda*(t^alpha) --> H_0^{-1}(w) = (w/lambda)^{1/alpha}  
## S(T) ~ U[0,1]
## T = H_0^{-1}( - log(U)/exp(b*x) ) =  {(1/lambda)*( - log(U)/exp(b*x) )}^{1/alpha}  

## -----------------------------------------------------------------------------------------------
library(MASS)
library(survival)
library(survminer)
library(glmnet)
## -----------------------------------------------------------------------------------------------
## ---------------- Simulation of Cox Model starts here:
## Simulate U from Uniform
## Simulate b*x
## Set parameters: lambda & alpha >0
## calculate T using T = {(1/lambda)*( - log(U)/exp(b*x) )}^{1/alpha} 
## repeat above steps n times to get a sample with size n

sim.total.n <- 10000
sim.summary <- data.frame(case.id=1:sim.total.n, lambda.value=NA,  alpha.value=NA,
                          "ST_5perc"=NA,  "ST_10perc"=NA,   "ST_15perc"=NA,   "ST_20perc"=NA,  
                          "ST_25perc"=NA,   "ST_30perc"=NA,   "ST_35perc"=NA,   "ST_40perc"=NA,   
                          "ST_45perc"=NA,   "ST_50perc"=NA,   "ST_55perc"=NA,   "ST_60perc"=NA,  
                          "ST_65perc"=NA,  "ST_70perc"=NA,   "ST_75perc"=NA,   "ST_80perc"=NA,   
                          "ST_85perc"=NA,   "ST_90perc"=NA,   "ST_95perc"=NA,
                          multivariate.Cox.FS.success.rate= NA, 
                          multivariate.Cox.effect.size.ranking.success.rate=NA)

for(case.id in 1:sim.total.n){
  ## --------- baseline hazard parameters:
  # lambda <- sample(c(10:200)/5,1) ## this is to make sure all lambda>1
  # alpha <- sample(c(10:200)/5,1)
  lambda <- sample(c(1:200)/50,1) ## ## this is to make about half of simulation data sets generated using lambda <2 
  alpha <- sample(c(1:200)/50,1)
  
  n <- 1000 ##  sample size
  p <- 10 ## number of covariates
  mu.value <- 0 ## the same mean for all features
  cov.value <- 0 ## the same covariance for all features
  
  U <- runif(n)
  b <- c(1:p) ## fixed coefficients
  
  mu.vector <- rep(mu.value, p)
  cov.matrix <- matrix(cov.value, ncol = p, nrow = p) + diag((1-cov.value), nrow = p, ncol = p)
  
  # -----------------------create multivariate normal distribution
  x <- mvrnorm(n,mu = mu.vector, Sigma = cov.matrix)
  xb <- x%*%b
  xb <- as.vector(xb)
  ST = ((1/lambda)*( - log(U)/exp(xb) ))^{1/alpha} 
  data <- data.frame(ST, x)
  
  data$OS.status <- 1   ## no censored data ## ifelse(runif(n)>0.9,  0,1)   ## 10% right-censored 
  data$OS.time <- data$ST
  
  data$OS.SurvObj <- with(data, Surv(OS.time, OS.status))
  
  sim.summary$lambda.value[case.id] <- lambda
  sim.summary$alpha.value[case.id] <- alpha
  sim.summary[case.id, 4:22] <- quantile(ST, c(1:19)/20)
  
  tryCatch({# paste(names(data), collapse = " + ")
    tmp.fit <- coxph(OS.SurvObj ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    
    feature.selection.success.rate <- length(which(fit.results[,5]<=0.05))/length(fit.results[,5])
    sim.summary$multivariate.Cox.FS.success.rate[case.id] <- feature.selection.success.rate
    sim.summary$multivariate.Cox.effect.size.ranking.success.rate[case.id]  <- length(which((rank(fit.results[,5]) == rank(c(1:p)*(-1)))))/p
  }, error = function(e){})
}

save(sim.summary, file="Cox_regression_simulation_n1000_sim_6_10082022summary_v2.RData")

sim.summary[which(sim.summary$multivariate.Cox.FS.success.rate<0.3)[1:5],]
#    case.id lambda.value alpha.value     ST_5perc    ST_10perc    ST_15perc    ST_20perc
# 9        9         2.14        0.84 9.350393e-18 3.568088e-14 1.617780e-11 7.442731e-10
# 10      10         1.66        0.64 4.026518e-23 2.777376e-17 1.003691e-14 5.436825e-12
# 11      11         1.92        0.18 3.022613e-81 2.056608e-63 7.527896e-53 1.955314e-43
# 12      12         1.18        0.24 1.775175e-64 2.796719e-50 1.383311e-40 3.493920e-33
# 28      28         2.94        0.24 2.763666e-60 1.476698e-47 1.420553e-39 2.056484e-32
#       ST_25perc    ST_30perc    ST_35perc    ST_40perc    ST_45perc    ST_50perc  ST_55perc
# 9  3.865715e-08 1.487240e-06 2.192648e-05 2.685619e-04 3.742111e-03 8.776112e-02   1.725300
# 10 1.351186e-09 4.477593e-07 2.942435e-05 8.398161e-04 1.066685e-01 3.048403e+00  88.414878
# 11 1.466334e-36 3.895275e-30 2.842997e-25 1.419553e-20 4.267093e-14 1.320888e-06  28.670434
# 12 4.890674e-28 5.685617e-23 6.840895e-18 2.404174e-12 2.684553e-07 1.569667e-03 110.522858
# 28 4.577938e-27 4.389652e-23 3.683524e-19 2.339028e-14 4.712968e-09 4.354585e-05   1.869009
#       ST_60perc    ST_65perc    ST_70perc    ST_75perc    ST_80perc    ST_85perc    ST_90perc
# 9  2.824537e+01 1.092661e+03 2.179292e+04 2.591864e+06 2.670544e+08 1.394949e+10 1.440242e+12
# 10 2.346497e+03 3.852983e+04 1.761503e+06 5.198526e+08 6.707559e+10 1.639497e+13 1.940117e+16
# 11 1.874775e+07 4.656959e+12 2.551883e+18 3.948259e+25 5.053743e+33 1.417695e+44 2.401564e+56
# 12 3.400801e+06 8.617177e+10 7.051396e+15 5.178826e+21 9.266643e+28 7.419179e+35 7.476413e+42
# 28 2.941203e+04 7.912463e+09 1.216338e+14 3.067514e+19 2.525893e+24 6.952676e+32 1.152569e+42
#       ST_95perc multivariate.Cox.FS.success.rate
# 9  1.942668e+15                                0
# 10 1.793770e+21                                0
# 11 2.841658e+72                                0
# 12 6.214227e+55                                0
# 28 4.463133e+55                                0
# multivariate.Cox.effect.size.ranking.success.rate
# 9                                                0.2
# 10                                               0.2
# 11                                               0.0
# 12                                               0.4
# 28                                               0.4
## ------------------
library(plotly)


tmp.data <- sim.summary[order(sim.summary$alpha.value),]
names(tmp.data)[23] <- "Sensitivity"

plot_ly(data=tmp.data,
        x=~alpha.value, 
        y=~lambda.value,
        z=~Sensitivity,
        type="scatter3d", mode="markers", color=sensitivity) 


