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
                          multivariate.Cox.effect.size.ranking.success.rate=NA,
                          multivariate.Cox.ph.check.p.of.max.effect.size=NA,
                          multivariate.Cox.ph.check.p.of.min.effect.size=NA,
                          multivariate.Cox.max.effect.size=NA,
                          multivariate.Cox.min.effect.size=NA,
                          ph.check.global.p=NA, 
                          multivariate.Cox.ph.false.violation.rate=NA, #  ph.check.all.p=NA,
                          ##-------:
                          univariate.Cox.FS.success.rate= NA,
                          univariate.Cox.effect.size.ranking.success.rate=NA,
                          univariate.Cox.ph.false.violation.rate=NA,
                          univariate.Cox.ph.check.p.of.max.effect.size=NA,
                          univariate.Cox.ph.check.p.of.min.effect.size=NA,
                          ##-------:
                          multivariate.Cox.FS.FDR=NA,
                          ph.check.random.global.p=NA,
                          multivariate.Cox.ph.random.violation.rate=NA,
                          univariate.Cox.FS.FDR=NA,
                          univariate.Cox.ph.random.violation.rate=NA,
                          ##--------:
                          univariate.logistic.regression.FS.success.rate=NA,
                          univariate.logistic.regression.effect.size.ranking.success.rate=NA,
                          univariate.logistic.regression.FDR=NA,
                          ##--------:
                          univariate.gaussian.regression.FS.success.rate=NA,
                          univariate.gaussian.regression.effect.size.ranking.success.rate=NA,
                          univariate.gaussian.regression.p.values=NA,
                          univariate.gaussian.regression.FDR=NA,
                          ##--------:
                          multivariate.gaussian.regression.effect.size.ranking.success.rate=NA,
                          multivariate.gaussian.regression.p.values=NA,
                          ##--------:
                          glmnet.Cox.lambda.min.FS.success.rate=NA,
                          glmnet.Cox.lambda.min.effect.size.ranking.success.rate=NA,
                          glmnet.Cox.lambda.min.FDR=NA,
                          ##--------:
                          glmnet.Cox.lambda.1se.FS.success.rate=NA,
                          glmnet.Cox.lambda.1se.effect.size.ranking.success.rate=NA,
                          glmnet.Cox.lambda.1se.beta=NA,
                          glmnet.Cox.lambda.1se.FDR=NA,
                          ##--------:
                          glmnet.Gaussian.FS.success.rate=NA,
                          glmnet.Gaussian.effect.size.ranking.success.rate=NA,
                          glmnet.Gaussian.FDR=NA)

for(case.id in 1:sim.total.n){
  ## --------- baseline hazard parameters:
  lambda <- sample(c(10:200)/5,1) ## this is to make sure all lambda>1
  alpha <- sample(c(10:200)/5,1)
  
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
  z <- mvrnorm(n,mu = mu.vector, Sigma = cov.matrix) ## --------  fake x
  colnames(z) <- paste("Z",1:p,sep="")
  xb <- x%*%b
  xb <- as.vector(xb)
  ST = ((1/lambda)*( - log(U)/exp(xb) ))^{1/alpha} 
  data <- data.frame(ST, x,z)
  
  ## --------------- true parameters:
  para.truth <- list(lambda,
                     alpha,
                     b,
                     mu.vector,
                     cov.matrix)
  names(para.truth) <- c("lambda",
                         "alpha",
                         "b",
                         "mu.vector",
                         "cov.matrix")
  # para.truth[1:4]
  
  ## --------------- try fitting Cox regression with no censored data ...
  #  head(data)
  data$OS.status <- ifelse(runif(n)>0.9,  0,1)   ## 10% right-censored ## 1   ## no censored data ## 
  data$OS.time <- data$ST
  data$OS.time[which(data$OS.status==0)] <- data$OS.time[which(data$OS.status==0)]*runif(length(which(data$OS.status==0)))
  data$OS.SurvObj <- with(data, Surv(OS.time, OS.status))
  
  # paste(names(data), collapse = " + ")
  #tmp.fit <- survfit(OS.SurvObj ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data, conf.type = "log-log")
  start_time <- Sys.time()
  #tmp.fit <- coxph(OS.SurvObj ~ X1 + X2 + X3 + X4 + X5, data)
  pp <- sample(1:p, p/2, replace = F)
  tmp.fit <- coxph(as.formula(paste("OS.SurvObj ~ ", paste(paste(rep("X",p/2),pp,sep=""), collapse = " + "),collapse = "")), data)
  end_time <- Sys.time()
  
  run.time <- end_time - start_time
  fit.results <- summary(tmp.fit)
  # print(fit.results)
  fit.results <- fit.results$coefficients
  feature.selection.success.rate <- length(which(fit.results[,5]<=0.05))/length(fit.results[,5])
  # print(feature.selection.success.rate)
  #tmp.fit
  
  # ## -----------
  # simulation.case <- list(para.truth,
  #                         simulation.quality,
  #                         fit.results,
  #                         run.time,
  #                         feature.selection.success.rate)
  # names(simulation.case) <- c("para.truth",
  #                             "simulation.quality",
  #                             "fit.results",
  #                             "fitting.time",
  #                             "feature.selection.success.rate")
  # 
  # save(simulation.case, file=paste("Explore_alpha_lambda_sim5_08032022_case_", case.id, ".RData", sep=""))
  # ## -----------
  
  #sim.summary$max.ST[case.id] <- max(ST)
  sim.summary$lambda.value[case.id] <- lambda
  sim.summary$alpha.value[case.id] <- alpha
  sim.summary[case.id, 4:22] <- quantile(ST, c(1:19)/20)

  sim.summary$multivariate.Cox.FS.success.rate[case.id] <- feature.selection.success.rate
  ph.check <- cox.zph(tmp.fit)
  #print(ph.check)
  ph.check.global.index <- dim(ph.check$table)[1]
  sim.summary$ph.check.global.p[case.id] <- ph.check$table[,"p"][ph.check.global.index]
  sim.summary$multivariate.Cox.ph.false.violation.rate[case.id] <- length(which(ph.check$table[,"p"][-ph.check.global.index]<0.05))/(ph.check.global.index-1)
  #print(sim.summary);
  
  #table(fit.results[,5]<=0.05, ph.check$table[,"p"][-ph.check.global.index]<0.05)
  #table(rank(fit.results[,5]), rank(ph.check$table[,"p"][-ph.check.global.index]))
  # table(fit.results[,5], pp)
  # table(rank(fit.results[,5]) == rank(pp*(-1)))
  sim.summary$multivariate.Cox.effect.size.ranking.success.rate[case.id]  <- length(which((rank(fit.results[,5]) == rank(pp*(-1)))))/length(pp)
  #readline("pause...")
  sim.summary$multivariate.Cox.ph.check.p.of.max.effect.size[case.id] <- fit.results[,5][which(pp==max(pp))]
  sim.summary$multivariate.Cox.max.effect.size[case.id] <- max(pp)
  sim.summary$multivariate.Cox.ph.check.p.of.min.effect.size[case.id] <- fit.results[,5][which(pp==min(pp))]
  sim.summary$multivariate.Cox.min.effect.size[case.id] <- min(pp)

  
  univariate.Cox.summary <- data.frame(var.name=paste("X",pp,sep=""), true.effect.size = pp, p.value=NA, ph.check.p=NA)
  for(j in pp){ ## ------------- fitting univariate Cox
    tmp.fit <- coxph(as.formula(paste("OS.SurvObj ~ ", paste("X",j,sep=""),collapse = "")), data)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    univariate.Cox.summary$p.value[which(univariate.Cox.summary$true.effect.size==j)] <- fit.results[,5]
    ph.check <- cox.zph(tmp.fit)
    univariate.Cox.summary$ph.check.p[which(univariate.Cox.summary$true.effect.size==j)] <- ph.check$table[,"p"][1]
  }
  
  sim.summary$univariate.Cox.FS.success.rate[case.id]  <- length(which(univariate.Cox.summary$p.value<0.05))/length(pp)
  sim.summary$univariate.Cox.effect.size.ranking.success.rate[case.id]  <- length(which((rank(univariate.Cox.summary$p.value) == rank(pp*(-1)))))/length(pp)
  
  sim.summary$univariate.Cox.ph.false.violation.rate[case.id]  <- length(which(univariate.Cox.summary$ph.check.p<0.05))/length(pp)
  sim.summary$univariate.Cox.ph.check.p.of.max.effect.size[case.id] <- univariate.Cox.summary$ph.check.p[which(univariate.Cox.summary$true.effect.size==max(univariate.Cox.summary$true.effect.size))]
  sim.summary$univariate.Cox.ph.check.p.of.min.effect.size[case.id] <- univariate.Cox.summary$ph.check.p[which(univariate.Cox.summary$true.effect.size==min(univariate.Cox.summary$true.effect.size))]

  ## ------------------------------ check false discovery rate of multivariate and univariate Cox regression:
  tmp.fit <- coxph(as.formula(paste("OS.SurvObj ~ ", paste(paste(rep("Z",p/2),pp,sep=""), collapse = " + "),collapse = "")), data)
  fit.results <- summary(tmp.fit)
  fit.results <- fit.results$coefficients
  false.discovery.rate <- length(which(fit.results[,5]<=0.05))/length(fit.results[,5])
  sim.summary$multivariate.Cox.FS.FDR[case.id] <- false.discovery.rate
  
  ph.check <- cox.zph(tmp.fit)
  ph.check.global.index <- dim(ph.check$table)[1]
  sim.summary$ph.check.random.global.p[case.id] <- ph.check$table[,"p"][ph.check.global.index]
  sim.summary$multivariate.Cox.ph.random.violation.rate[case.id] <- length(which(ph.check$table[,"p"][-ph.check.global.index]<0.05))/(ph.check.global.index-1)
  
  univariate.Cox.summary <- data.frame(var.name=paste("Z",pp,sep=""), var.index = pp, p.value=NA, ph.check.p=NA)
  for(j in pp){ ## ------------- fitting univariate Cox with random noise:
    tmp.fit <- coxph(as.formula(paste("OS.SurvObj ~ ", paste("Z",j,sep=""),collapse = "")), data)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    univariate.Cox.summary$p.value[which(univariate.Cox.summary$var.index==j)] <- fit.results[,5]
    ph.check <- cox.zph(tmp.fit)
    univariate.Cox.summary$ph.check.p[which(univariate.Cox.summary$var.index==j)] <- ph.check$table[,"p"][1]
  }
  
  sim.summary$univariate.Cox.FS.FDR[case.id]  <- length(which(univariate.Cox.summary$p.value<0.05))/length(pp)
  sim.summary$univariate.Cox.ph.random.violation.rate[case.id]  <- length(which(univariate.Cox.summary$ph.check.p<0.05))/length(pp)

  ## ----------------------:
  univariate.logistic.regression.summary <- data.frame(var.name=paste("X",pp,sep=""), true.effect.size = pp, p.value=NA)
  for(j in pp){ ## ------------- fitting univariate logistic.regression
    tmp.fit <- glm(as.formula(paste("OS.status ~ OS.time + ", paste("X",j,sep=""),collapse = "")), data, family=binomial)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    univariate.logistic.regression.summary$p.value[which(univariate.logistic.regression.summary$true.effect.size==j)] <- fit.results[3,4]
  }
  
  sim.summary$univariate.logistic.regression.FS.success.rate[case.id]  <- length(which(univariate.logistic.regression.summary$p.value<0.05))/length(pp)
  sim.summary$univariate.logistic.regression.effect.size.ranking.success.rate[case.id]  <- length(which((rank(univariate.logistic.regression.summary$p.value) == rank(pp*(-1)))))/length(pp)
  
  ## ----------------------:
  univariate.logistic.regression.summary <- data.frame(var.name=paste("Z",pp,sep=""), var.index = pp, p.value=NA)
  for(j in pp){ ## ------------- fitting logistic.regression with random noise
    tmp.fit <- glm(as.formula(paste("OS.status ~ OS.time + ", paste("Z",j,sep=""),collapse = "")), data, family=binomial)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    univariate.logistic.regression.summary$p.value[which(univariate.logistic.regression.summary$var.index==j)] <- fit.results[3,4]
  }
  
  sim.summary$univariate.logistic.regression.FDR[case.id]  <- length(which(univariate.logistic.regression.summary$p.value<0.05))/length(pp)
  
  ## ----------------------:
  univariate.gaussian.regression.summary <- data.frame(var.name=paste("X",pp,sep=""), true.effect.size = pp, p.value=NA)
  for(j in pp){ ## ------------- fitting gaussian.regression
    tmp.fit <- lm(as.formula(paste("log(OS.time) ~ OS.status + ", paste("X",j,sep=""),collapse = "")), data)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    univariate.gaussian.regression.summary$p.value[which(univariate.gaussian.regression.summary$true.effect.size==j)] <- fit.results[3,4]
  }
  
  sim.summary$univariate.gaussian.regression.FS.success.rate[case.id]  <- length(which(univariate.gaussian.regression.summary$p.value<0.05))/length(pp)
  sim.summary$univariate.gaussian.regression.effect.size.ranking.success.rate[case.id]  <- length(which((rank(univariate.gaussian.regression.summary$p.value) == rank(pp*(-1)))))/length(pp)
  sim.summary$univariate.gaussian.regression.p.values[case.id] <- paste(paste(univariate.gaussian.regression.summary$var.name,
                                                                        ": ",
                                                                        univariate.gaussian.regression.summary$p.value,
                                                                        sep=""), collapse=";  " )
  
  if(length(which(univariate.gaussian.regression.summary$p.value<0.05))>0){
    selected.feature <- as.character(univariate.gaussian.regression.summary$var.name[which(univariate.gaussian.regression.summary$p.value<0.05)])
  }else{selected.feature <- c("No X was selected.")}

  ## ----------------------:
  univariate.gaussian.regression.summary <- data.frame(var.name=paste("Z",pp,sep=""), var.index = pp, p.value=NA)
  for(j in pp){ ## ------------- fitting gaussian.regression with random noise
    tmp.fit <- lm(as.formula(paste("log(OS.time) ~ OS.status + ", paste("Z",j,sep=""),collapse = "")), data)
    fit.results <- summary(tmp.fit)
    fit.results <- fit.results$coefficients
    univariate.gaussian.regression.summary$p.value[which(univariate.gaussian.regression.summary$var.index==j)] <- fit.results[3,4]
  }
  
  sim.summary$univariate.gaussian.regression.FDR[case.id]  <- length(which(univariate.gaussian.regression.summary$p.value<0.05))/length(pp)
  
  if(length(which(univariate.gaussian.regression.summary$p.value<0.05))>0){
    selected.feature <- c(selected.feature,
                          as.character(univariate.gaussian.regression.summary$var.name[which(univariate.gaussian.regression.summary$p.value<0.05)])
    )
  }else{selected.feature <- c(selected.feature, "No Z was selected.")}
  
  ## ------------- fitting multivariate gaussian.regression with features selected by Gaussian model with 2 covariates:
  selected.feature <- data.frame(var.name=names(data)[which(names(data)%in%selected.feature)])
  selected.feature$var.category <- substr(selected.feature$var.name, 1,1)
  selected.feature$true.effect.size <- as.integer(substring(selected.feature$var.name, 2))
  if(length(which(selected.feature$var.category=="Z"))>0){
    selected.feature$true.effect.size[which(selected.feature$var.category=="Z")] <- 0
  }

  tmp.fit <- lm(as.formula(paste("log(OS.time) ~ OS.status + ", paste(selected.feature$var.name,collapse = " + "),collapse = "")), data)
  fit.results <- summary(tmp.fit)
  fit.results <- fit.results$coefficients
  selected.feature$multivariate.Gaussian.p.value <- fit.results[-c(1:2),4]
  
  multivariate.gaussian.regression.summary <- data.frame(var.name=paste("X",pp,sep=""), true.effect.size = pp)
  multivariate.gaussian.regression.summary <- merge(multivariate.gaussian.regression.summary,
                                                    selected.feature[,c(1,4)],
                                                    by="var.name", all.x = T)
  if(length(which(is.na(multivariate.gaussian.regression.summary$multivariate.Gaussian.p.value)))>0){
    multivariate.gaussian.regression.summary$multivariate.Gaussian.p.value[which(is.na(multivariate.gaussian.regression.summary$multivariate.Gaussian.p.value))] <- 0.05
  }
  
  sim.summary$multivariate.gaussian.regression.effect.size.ranking.success.rate[case.id] <- length(which((rank(multivariate.gaussian.regression.summary$multivariate.Gaussian.p.value) == rank(multivariate.gaussian.regression.summary$true.effect.size*(-1)))))/5
  sim.summary$multivariate.gaussian.regression.p.values[case.id] <- paste(paste(row.names(fit.results),
                                                                                ": ",
                                                                                fit.results[,4],
                                                                                sep=""), collapse=";  " )
  
  ## ------------------------ fit lasso regularized Cox:
  cox.path.X <- as.matrix(data[,which(names(data)%in%c(paste("X", pp, sep=""), paste("Z", pp, sep="")))])
  y = cbind(time = data$OS.time, status = data$OS.status)
  foldid = sample(rep(seq(10), length = n))
  fit1_cv = cv.glmnet(cox.path.X, y, family = "cox", foldid = foldid)
  lambda.min.v <- fit1_cv$lambda.min  ## value of lambda that gives minimum "The mean cross-validated error".
  lambda.1se.v <- fit1_cv$lambda.1se  ## largest value of lambda such that error is within 1 standard error of the mini- mum.
  # plot(fit1_cv)
  # title("Cox Family", line = 2.5)
  cox.path.jsurv <- data$OS.SurvObj
  tmp.fit.1 <- glmnet:::cox.path(cox.path.X, cox.path.jsurv, lambda = lambda.min.v)
  tmp.fit.2 <- glmnet:::cox.path(cox.path.X, cox.path.jsurv, lambda = lambda.1se.v)
  tmp.fit.1$beta
  tmp.fit.2$beta
  
  sim.summary$glmnet.Cox.lambda.min.FS.success.rate[case.id]  <- length(which(tmp.fit.1$beta[1:5]!=0))/5
  sim.summary$glmnet.Cox.lambda.min.effect.size.ranking.success.rate[case.id]  <-  length(which(c(rank(abs(tmp.fit.1$beta[5:1])) == 5:1)))/5
  sim.summary$glmnet.Cox.lambda.min.FDR[case.id] <- length(which(tmp.fit.1$beta[6:10]!=0))/5
  
  sim.summary$glmnet.Cox.lambda.1se.FS.success.rate[case.id]  <- length(which(tmp.fit.2$beta[1:5]!=0))/5
  sim.summary$glmnet.Cox.lambda.1se.effect.size.ranking.success.rate[case.id]  <- length(which(c(rank(abs(tmp.fit.2$beta[5:1])) == 5:1)))/5
  sim.summary$glmnet.Cox.lambda.1se.beta[case.id] <- paste(paste(rownames(tmp.fit.2$beta)[1:5],
                                                                 ": ",
                                                                 tmp.fit.2$beta[1:5],
                                                                 sep=""), collapse=";  " )
  sim.summary$glmnet.Cox.lambda.1se.FDR[case.id] <- length(which(tmp.fit.2$beta[6:10]!=0))/5
  
  ## ------------------------ fit Gaussian model with elastic net regularization:
  #cv.glmnet(cox.path.X, y)
  cox.path.X <- as.matrix(data[,which(names(data)%in%c(paste("X", pp, sep=""), paste("Z", pp, sep="")))])
  cox.path.X <- cbind(cox.path.X, data$OS.status)
  y <- data$OS.time
  tmp.fit <- glmnet(cox.path.X, log(y))
  tmp.beta <- coef(tmp.fit, s = 0.05)
  tmp.beta <- tmp.beta[-1]
  
  sim.summary$glmnet.Gaussian.FS.success.rate[case.id]  <- length(which(tmp.beta[1:5]!=0))/5
  sim.summary$glmnet.Gaussian.effect.size.ranking.success.rate[case.id]  <- length(which(c(rank(abs(tmp.beta[5:1])) == 5:1)))/5
  sim.summary$glmnet.Gaussian.FDR[case.id] <- length(which(tmp.beta[6:10]!=0))/5
}

save(sim.summary, file="Cox_regression_simulation_n1000_sim_7_09172022summary.RData")



## --------------
multivariate.Cox.box.plot.data.tmp <- sim.summary[,c("multivariate.Cox.max.effect.size","multivariate.Cox.ph.check.p.of.max.effect.size")]
names(multivariate.Cox.box.plot.data.tmp) <- c("true.effect.size", "multivariate.Cox.ph.check.p.value")
#head(multivariate.Cox.box.plot.data.tmp)
multivariate.Cox.box.plot.data <- multivariate.Cox.box.plot.data.tmp
multivariate.Cox.box.plot.data.tmp <- sim.summary[,c("multivariate.Cox.min.effect.size","multivariate.Cox.ph.check.p.of.min.effect.size")]
names(multivariate.Cox.box.plot.data.tmp) <- c("true.effect.size", "multivariate.Cox.ph.check.p.value")
#head(multivariate.Cox.box.plot.data.tmp)
multivariate.Cox.box.plot.data <- rbind(multivariate.Cox.box.plot.data, multivariate.Cox.box.plot.data.tmp)
#table(multivariate.Cox.box.plot.data$true.effect.size)

univariate.Cox.box.plot.data.tmp <- sim.summary[,c("multivariate.Cox.max.effect.size","univariate.Cox.ph.check.p.of.max.effect.size")]
names(univariate.Cox.box.plot.data.tmp) <- c("true.effect.size", "univariate.Cox.ph.check.p.value")
#head(univariate.Cox.box.plot.data.tmp)
univariate.Cox.box.plot.data <- univariate.Cox.box.plot.data.tmp
univariate.Cox.box.plot.data.tmp <- sim.summary[,c("multivariate.Cox.min.effect.size","univariate.Cox.ph.check.p.of.min.effect.size")]
names(univariate.Cox.box.plot.data.tmp) <- c("true.effect.size", "univariate.Cox.ph.check.p.value")
#head(univariate.Cox.box.plot.data.tmp)
univariate.Cox.box.plot.data <- rbind(univariate.Cox.box.plot.data, univariate.Cox.box.plot.data.tmp)
#table(univariate.Cox.box.plot.data$true.effect.size)

## --------------
sink("sim_7_09172022_n1000_results_summary.txt")
head(sim.summary)
cat("\n");cat("----------------------- Sensitivity ---------------------------------\n");
print("table(sim.summary$multivariate.Cox.FS.success.rate)");
table(sim.summary$multivariate.Cox.FS.success.rate)
prop.table(table(sim.summary$multivariate.Cox.FS.success.rate))
cat("\n");cat("\n");print("table(sim.summary$univariate.Cox.FS.success.rate)");
table(sim.summary$univariate.Cox.FS.success.rate)
prop.table(table(sim.summary$univariate.Cox.FS.success.rate))
cat("\n");cat("\n");print("table(sim.summary$univariate.logistic.regression.FS.success.rate)");
table(sim.summary$univariate.logistic.regression.FS.success.rate)
prop.table(table(sim.summary$univariate.logistic.regression.FS.success.rate))
cat("\n");cat("\n");print("table(sim.summary$univariate.gaussian.regression.FS.success.rate)");
table(sim.summary$univariate.gaussian.regression.FS.success.rate)
prop.table(table(sim.summary$univariate.gaussian.regression.FS.success.rate))
cat("\n");cat("\n");print("table(sim.summary$glmnet.Cox.lambda.min.FS.success.rate)");
table(sim.summary$glmnet.Cox.lambda.min.FS.success.rate)
prop.table(table(sim.summary$glmnet.Cox.lambda.min.FS.success.rate))
cat("\n");cat("\n");print("table(sim.summary$glmnet.Cox.lambda.1se.FS.success.rate)");
table(sim.summary$glmnet.Cox.lambda.1se.FS.success.rate)
prop.table(table(sim.summary$glmnet.Cox.lambda.1se.FS.success.rate))
## -------:
cat("\n");cat("--------------------------Effect size ranking accuracy------------------------------\n");
print("prop.table(table(sim.summary$multivariate.Cox.effect.size.ranking.success.rate))");
table(sim.summary$multivariate.Cox.effect.size.ranking.success.rate)
prop.table(table(sim.summary$multivariate.Cox.effect.size.ranking.success.rate))
cat("\n");print("prop.table(table(sim.summary$univariate.Cox.effect.size.ranking.success.rate))");
table(sim.summary$univariate.Cox.effect.size.ranking.success.rate)
prop.table(table(sim.summary$univariate.Cox.effect.size.ranking.success.rate))
cat("\n");print("prop.table(table(sim.summary$univariate.logistic.regression.effect.size.ranking.success.rate))");
table(sim.summary$univariate.logistic.regression.effect.size.ranking.success.rate)
prop.table(table(sim.summary$univariate.logistic.regression.effect.size.ranking.success.rate))
cat("\n");print("prop.table(table(sim.summary$univariate.gaussian.regression.effect.size.ranking.success.rate))");
table(sim.summary$univariate.gaussian.regression.effect.size.ranking.success.rate)
prop.table(table(sim.summary$univariate.gaussian.regression.effect.size.ranking.success.rate))
cat("\n");print("prop.table(table(sim.summary$multivariate.gaussian.regression.effect.size.ranking.success.rate))");
table(sim.summary$multivariate.gaussian.regression.effect.size.ranking.success.rate)
prop.table(table(sim.summary$multivariate.gaussian.regression.effect.size.ranking.success.rate))
cat("\n");print("prop.table(table(sim.summary$glmnet.Cox.lambda.min.effect.size.ranking.success.rate))");
table(sim.summary$glmnet.Cox.lambda.min.effect.size.ranking.success.rate)
prop.table(table(sim.summary$glmnet.Cox.lambda.min.effect.size.ranking.success.rate))
cat("\n");print("prop.table(table(sim.summary$glmnet.Cox.lambda.1se.effect.size.ranking.success.rate))");
table(sim.summary$glmnet.Cox.lambda.1se.effect.size.ranking.success.rate)
prop.table(table(sim.summary$glmnet.Cox.lambda.1se.effect.size.ranking.success.rate))
## -------:
cat("\n");cat("--------------------------Specificity------------------------------\n");
print("table(sim.summary$multivariate.Cox.FS.FDR)");
table(sim.summary$multivariate.Cox.FS.FDR)
prop.table(table(sim.summary$multivariate.Cox.FS.FDR))
cat("\n");cat("\n");print("table(sim.summary$univariate.Cox.FS.FDR)");
table(sim.summary$univariate.Cox.FS.FDR)
prop.table(table(sim.summary$univariate.Cox.FS.FDR))
cat("\n");cat("\n");print("table(sim.summary$univariate.logistic.regression.FDR)");
table(sim.summary$univariate.logistic.regression.FDR)
prop.table(table(sim.summary$univariate.logistic.regression.FDR))
cat("\n");cat("\n");print("table(sim.summary$univariate.gaussian.regression.FDR)");
table(sim.summary$univariate.gaussian.regression.FDR)
prop.table(table(sim.summary$univariate.gaussian.regression.FDR))
cat("\n");cat("\n");print("table(sim.summary$glmnet.Cox.lambda.min.FDR)");
table(sim.summary$glmnet.Cox.lambda.min.FDR)
prop.table(table(sim.summary$glmnet.Cox.lambda.min.FDR))
cat("\n");cat("\n");print("table(sim.summary$glmnet.Cox.lambda.1se.FDR)");
table(sim.summary$glmnet.Cox.lambda.1se.FDR)
prop.table(table(sim.summary$glmnet.Cox.lambda.1se.FDR))
## -------------------------------------------------------------------------------------:
cat("\n");cat("\n");
sapply(sim.summary[,-c(1:2)], summary)
## -------:
cat("\n");cat("--------------------------------------------------------\n");
print("table(sim.summary$multivariate.Cox.ph.false.violation.rate)");
table(sim.summary$multivariate.Cox.ph.false.violation.rate)
prop.table(table(sim.summary$multivariate.Cox.ph.false.violation.rate))
cat("\n");cat("\n");print("table(sim.summary$univariate.Cox.ph.false.violation.rate)");
table(sim.summary$univariate.Cox.ph.false.violation.rate)
prop.table(table(sim.summary$univariate.Cox.ph.false.violation.rate))
## -------:
cat("\n");cat("\n");print("prop.table(table(sim.summary$ph.check.global.p<0.05))");
table(sim.summary$ph.check.global.p<0.05)
prop.table(table(sim.summary$ph.check.global.p<0.05))
## -------:
cat("\n");cat("--------------------------------------------------------\n");
summary(glm(multivariate.Cox.ph.check.p.value ~ true.effect.size, data=multivariate.Cox.box.plot.data))
summary(glm(univariate.Cox.ph.check.p.value ~ true.effect.size, data=univariate.Cox.box.plot.data))
## -------:
multivariate.Cox.box.plot.data$true.effect.size <- as.factor(multivariate.Cox.box.plot.data$true.effect.size)
univariate.Cox.box.plot.data$true.effect.size <- as.factor(univariate.Cox.box.plot.data$true.effect.size)
cat("\n");print("tapply(multivariate.Cox.box.plot.data$multivariate.Cox.ph.check.p.value, multivariate.Cox.box.plot.data$true.effect.size, summary)");
tapply(multivariate.Cox.box.plot.data$multivariate.Cox.ph.check.p.value, multivariate.Cox.box.plot.data$true.effect.size, summary)
cat("\n");print("tapply(univariate.Cox.box.plot.data$univariate.Cox.ph.check.p.value, univariate.Cox.box.plot.data$true.effect.size, summary)");
tapply(univariate.Cox.box.plot.data$univariate.Cox.ph.check.p.value, univariate.Cox.box.plot.data$true.effect.size, summary)
## -------:
cat("\n");cat("--------------------------------------------------------\n");
print("table(sim.summary$multivariate.Cox.ph.random.violation.rate)");
table(sim.summary$multivariate.Cox.ph.random.violation.rate)
prop.table(table(sim.summary$multivariate.Cox.ph.random.violation.rate))
cat("\n");cat("\n");print("table(sim.summary$univariate.Cox.ph.random.violation.rate)");
table(sim.summary$univariate.Cox.ph.random.violation.rate)
prop.table(table(sim.summary$univariate.Cox.ph.random.violation.rate))
## -------:
cat("\n");print("prop.table(table(sim.summary$ph.check.random.global.p<0.05))");
table(sim.summary$ph.check.random.global.p<0.05)
prop.table(table(sim.summary$ph.check.random.global.p<0.05))
## -------:
sink()

boxplot(multivariate.Cox.ph.check.p.value ~ true.effect.size, data=multivariate.Cox.box.plot.data)
boxplot(univariate.Cox.ph.check.p.value ~ true.effect.size, data=univariate.Cox.box.plot.data)



