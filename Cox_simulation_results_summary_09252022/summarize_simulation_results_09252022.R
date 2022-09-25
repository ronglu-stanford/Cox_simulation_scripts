rm(list=ls())
## -----------------------------
setwd("./Cox_simulation_scripts_09192022_results")
all.files <- list.files()
all.RData.files <- all.files[grep("RData",  all.files)]
all.RData.files <- all.RData.files[c(1,5,9,13,
                                     3,7,11,15,
                                     2,6,10,14,
                                     4,8,12,16)]

all.RData.files
## -----------------------------
load(all.RData.files[1])
sim.summary.tb <- data.frame(model.name=c("Multivariate Cox",
                                          "Univariate Cox",
                                          "Logistic Regression",
                                          "Gaussian Regression (1 feature)",
                                          "glmnet Cox (lambda min)",
                                          "glmnet Cox (lambda 1se)",
                                          "glmnet Gaussian (lambda=0.05)",
                                          "Gaussian Regression (selected features)"#,
))
sim.summary.tb$Sensitivity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.success.rate",
                                                      "univariate.Cox.FS.success.rate",
                                                      "univariate.logistic.regression.FS.success.rate",
                                                      "univariate.gaussian.regression.FS.success.rate",
                                                      "glmnet.Cox.lambda.min.FS.success.rate",
                                                      "glmnet.Cox.lambda.1se.FS.success.rate",
                                                      "glmnet.Gaussian.FS.success.rate")], 
                                       function(x){
                                         tmp.tb <- table(x)/10000
                                         if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
                                       }),
                                NA)
sim.summary.tb$Specificity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.FDR",
                                                      "univariate.Cox.FS.FDR",
                                                      "univariate.logistic.regression.FDR",
                                                      "univariate.gaussian.regression.FDR",
                                                      "glmnet.Cox.lambda.min.FDR",
                                                      "glmnet.Cox.lambda.1se.FDR",
                                                      "glmnet.Gaussian.FDR")], 
                                       function(x){
                                         tmp.tb <- table(x)/10000
                                         if("0" %in% names(tmp.tb)){return(tmp.tb["0"])}else(return(0))
                                       }),
                                NA)
sim.summary.tb$Balanced.Accuracy <- (sim.summary.tb$Sensitivity + sim.summary.tb$Specificity)/2

sim.summary.tb$Accuracy.of.effect.size.ranking <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                                        "univariate.Cox.effect.size.ranking.success.rate",
                                                                        "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                                        "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                                        "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                                        "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                                        "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                                        "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
)], 
function(x){
  tmp.tb <- table(x)/10000
  if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
})
sim.summary.tb$fitting.error.rate <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                           "univariate.Cox.effect.size.ranking.success.rate",
                                                           "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                           "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                           "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                           "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                           "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                           "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
)], 
function(x){
  tmp.tb <- length(which(is.na(x)))/10000
  return(tmp.tb)
})
sim.summary.tb$file.name <- all.RData.files[1]
sim.summary.tb
saved.sim.summary.tb <- sim.summary.tb

## --------------------------------:
for(i in 2:length(all.RData.files)){
  load(all.RData.files[i])
  sim.summary.tb <- data.frame(model.name=c("Multivariate Cox",
                                            "Univariate Cox",
                                            "Logistic Regression",
                                            "Gaussian Regression (1 feature)",
                                            "glmnet Cox (lambda min)",
                                            "glmnet Cox (lambda 1se)",
                                            "glmnet Gaussian (lambda=0.05)",
                                            "Gaussian Regression (selected features)"#,
  ))
  sim.summary.tb$Sensitivity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.success.rate",
                                                        "univariate.Cox.FS.success.rate",
                                                        "univariate.logistic.regression.FS.success.rate",
                                                        "univariate.gaussian.regression.FS.success.rate",
                                                        "glmnet.Cox.lambda.min.FS.success.rate",
                                                        "glmnet.Cox.lambda.1se.FS.success.rate",
                                                        "glmnet.Gaussian.FS.success.rate")], 
                                         function(x){
                                           tmp.tb <- table(x)/10000
                                           if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
                                         }),
                                  NA)
  sim.summary.tb$Specificity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.FDR",
                                                        "univariate.Cox.FS.FDR",
                                                        "univariate.logistic.regression.FDR",
                                                        "univariate.gaussian.regression.FDR",
                                                        "glmnet.Cox.lambda.min.FDR",
                                                        "glmnet.Cox.lambda.1se.FDR",
                                                        "glmnet.Gaussian.FDR")], 
                                         function(x){
                                           tmp.tb <- table(x)/10000
                                           if("0" %in% names(tmp.tb)){return(tmp.tb["0"])}else(return(0))
                                         }),
                                  NA)
  sim.summary.tb$Balanced.Accuracy <- (sim.summary.tb$Sensitivity + sim.summary.tb$Specificity)/2
  
  sim.summary.tb$Accuracy.of.effect.size.ranking <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                                          "univariate.Cox.effect.size.ranking.success.rate",
                                                                          "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                                          "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                                          "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                                          "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                                          "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                                          "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
  )], 
  function(x){
    tmp.tb <- table(x)/10000
    if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
  })
  sim.summary.tb$fitting.error.rate <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                             "univariate.Cox.effect.size.ranking.success.rate",
                                                             "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                             "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                             "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                             "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                             "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                             "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
  )], 
  function(x){
    tmp.tb <- length(which(is.na(x)))/10000
    return(tmp.tb)
  })
  
  sim.summary.tb$file.name <- all.RData.files[i]
  saved.sim.summary.tb <- rbind(saved.sim.summary.tb, sim.summary.tb)
}

## the way that the simulation scripts fit the generic multivariate Cox regression is not appropriate for evaluating sensitivity and specificity:
saved.sim.summary.tb$Sensitivity[which(saved.sim.summary.tb$model.name=="Multivariate Cox")] <- NA
saved.sim.summary.tb$Specificity[which(saved.sim.summary.tb$model.name=="Multivariate Cox")] <- NA
saved.sim.summary.tb$Balanced.Accuracy[which(saved.sim.summary.tb$model.name=="Multivariate Cox")] <- NA
saved.sim.summary.tb[1:10,]
## -----------------------------
saved.sim.summary.tb.part1 <- saved.sim.summary.tb

# ----------------------------------------------------
setwd("./Cox_simulation_scripts_09172022_results")
all.files <- list.files()
all.RData.files <- all.files[grep("RData",  all.files)]
all.RData.files <- all.RData.files[c(9,1,5,13,
                                     11,3,7,15,
                                     10,2,6,14,
                                     12,4,8,16)]
all.RData.files
## -----------------------------
load(all.RData.files[1])
sim.summary.tb <- data.frame(model.name=c("Multivariate Cox",
                                          "Univariate Cox",
                                          "Logistic Regression",
                                          "Gaussian Regression (1 feature)",
                                          "glmnet Cox (lambda min)",
                                          "glmnet Cox (lambda 1se)",
                                          "glmnet Gaussian (lambda=0.05)",
                                          "Gaussian Regression (selected features)"#,
                                          ))
sim.summary.tb$Sensitivity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.success.rate",
                                                      "univariate.Cox.FS.success.rate",
                                                    "univariate.logistic.regression.FS.success.rate",
                                                    "univariate.gaussian.regression.FS.success.rate",
                                                    "glmnet.Cox.lambda.min.FS.success.rate",
                                                    "glmnet.Cox.lambda.1se.FS.success.rate",
                                                    "glmnet.Gaussian.FS.success.rate")], 
                                     function(x){
                                       tmp.tb <- table(x)/10000
                                       if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
                                     }),
                                NA)
sim.summary.tb$Specificity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.FDR",
                                                      "univariate.Cox.FS.FDR",
                                                      "univariate.logistic.regression.FDR",
                                                      "univariate.gaussian.regression.FDR",
                                                      "glmnet.Cox.lambda.min.FDR",
                                                      "glmnet.Cox.lambda.1se.FDR",
                                                      "glmnet.Gaussian.FDR")], 
                                       function(x){
                                         tmp.tb <- table(x)/10000
                                         if("0" %in% names(tmp.tb)){return(tmp.tb["0"])}else(return(0))
                                       }),
                                NA)
sim.summary.tb$Balanced.Accuracy <- (sim.summary.tb$Sensitivity + sim.summary.tb$Specificity)/2

sim.summary.tb$Accuracy.of.effect.size.ranking <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                                        "univariate.Cox.effect.size.ranking.success.rate",
                                                                        "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                                        "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                                        "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                                        "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                                        "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                                        "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
                                                                        )], 
                                                         function(x){
                                                           tmp.tb <- table(x)/10000
                                                           if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
                                                         })
sim.summary.tb$fitting.error.rate <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                           "univariate.Cox.effect.size.ranking.success.rate",
                                                           "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                           "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                           "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                           "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                           "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                           "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
                                                           )], 
                                            function(x){
                                              tmp.tb <- length(which(is.na(x)))/10000
                                              return(tmp.tb)
                                            })
sim.summary.tb$file.name <- all.RData.files[1]
sim.summary.tb
saved.sim.summary.tb <- sim.summary.tb

## --------------------------------:
for(i in 2:length(all.RData.files)){
  load(all.RData.files[i])
  sim.summary.tb <- data.frame(model.name=c("Multivariate Cox",
                                            "Univariate Cox",
                                            "Logistic Regression",
                                            "Gaussian Regression (1 feature)",
                                            "glmnet Cox (lambda min)",
                                            "glmnet Cox (lambda 1se)",
                                            "glmnet Gaussian (lambda=0.05)",
                                            "Gaussian Regression (selected features)"#,
  ))
  sim.summary.tb$Sensitivity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.success.rate",
                                                        "univariate.Cox.FS.success.rate",
                                                        "univariate.logistic.regression.FS.success.rate",
                                                        "univariate.gaussian.regression.FS.success.rate",
                                                        "glmnet.Cox.lambda.min.FS.success.rate",
                                                        "glmnet.Cox.lambda.1se.FS.success.rate",
                                                        "glmnet.Gaussian.FS.success.rate")], 
                                         function(x){
                                           tmp.tb <- table(x)/10000
                                           if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
                                         }),
                                  NA)
  sim.summary.tb$Specificity <- c(sapply(sim.summary[,c("multivariate.Cox.FS.FDR",
                                                        "univariate.Cox.FS.FDR",
                                                        "univariate.logistic.regression.FDR",
                                                        "univariate.gaussian.regression.FDR",
                                                        "glmnet.Cox.lambda.min.FDR",
                                                        "glmnet.Cox.lambda.1se.FDR",
                                                        "glmnet.Gaussian.FDR")], 
                                         function(x){
                                           tmp.tb <- table(x)/10000
                                           if("0" %in% names(tmp.tb)){return(tmp.tb["0"])}else(return(0))
                                         }),
                                  NA)
  sim.summary.tb$Balanced.Accuracy <- (sim.summary.tb$Sensitivity + sim.summary.tb$Specificity)/2
  
  sim.summary.tb$Accuracy.of.effect.size.ranking <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                                          "univariate.Cox.effect.size.ranking.success.rate",
                                                                          "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                                          "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                                          "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                                          "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                                          "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                                          "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
  )], 
  function(x){
    tmp.tb <- table(x)/10000
    if("1" %in% names(tmp.tb)){return(tmp.tb["1"])}else(return(0))
  })
  sim.summary.tb$fitting.error.rate <- sapply(sim.summary[,c("multivariate.Cox.effect.size.ranking.success.rate",
                                                             "univariate.Cox.effect.size.ranking.success.rate",
                                                             "univariate.logistic.regression.effect.size.ranking.success.rate",
                                                             "univariate.gaussian.regression.effect.size.ranking.success.rate",
                                                             "glmnet.Cox.lambda.min.effect.size.ranking.success.rate",
                                                             "glmnet.Cox.lambda.1se.effect.size.ranking.success.rate",
                                                             "glmnet.Gaussian.effect.size.ranking.success.rate",
                                                             "multivariate.gaussian.regression.effect.size.ranking.success.rate"#,
                                                             )], 
  function(x){
    tmp.tb <- length(which(is.na(x)))/10000
    return(tmp.tb)
  })
  
  sim.summary.tb$file.name <- all.RData.files[i]
  saved.sim.summary.tb <- rbind(saved.sim.summary.tb, sim.summary.tb)
}

## the way that the simulation scripts fit the generic multivariate Cox regression is not appropriate for evaluating sensitivity and specificity:
saved.sim.summary.tb$Sensitivity[which(saved.sim.summary.tb$model.name=="Multivariate Cox")] <- NA
saved.sim.summary.tb$Specificity[which(saved.sim.summary.tb$model.name=="Multivariate Cox")] <- NA
saved.sim.summary.tb$Balanced.Accuracy[which(saved.sim.summary.tb$model.name=="Multivariate Cox")] <- NA
saved.sim.summary.tb[1:10,]
## -----------------------------
saved.sim.summary.tb.part2 <- saved.sim.summary.tb

saved.sim.summary.tb <- rbind(saved.sim.summary.tb.part1, saved.sim.summary.tb.part2)
write.csv(saved.sim.summary.tb, file="saved_simulation_results_summary_09252022.csv", na="--", row.names = F)
