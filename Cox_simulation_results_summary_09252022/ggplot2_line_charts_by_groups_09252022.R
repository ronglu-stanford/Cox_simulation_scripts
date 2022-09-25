rm(list=ls())

library(ggplot2)

data <- read.csv("saved_simulation_results_summary_09252022v2.csv", stringsAsFactors = F, na.strings = c("--"))
dim(data)
# [1] 256  10

head(data)
data$event.sample.size.factor <- factor(as.character(data$event.sample.size), levels=as.character(sort(unique(data$event.sample.size))))
data$model <- data$model.name
data$model[which(data$model=="glmnet Gaussian (lambda=0.05)")] <- "Gaussian Elastic Net (lambda=0.05)" 
data$model.v2 <- data$model

tmp.data <- data[which(data$feature.corr==0),]
plot.1 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Accuracy.of.effect.size.ranking, group=model, color=model)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nProbability of correctly ranking the effects of all true features\n") + 
  ggtitle("\n\nAll features are independent") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "independent_features_effect_size_ranking_09252022.pdf" (10.65 X 6.10)

tmp.data <- data[which(data$feature.corr==0.8),]
plot.2 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Accuracy.of.effect.size.ranking, group=model, color=model)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nProbability of correctly ranking the effects of all true features\n") + 
  ggtitle("\n\nAll true features are highly correlated") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "correlated_features_effect_size_ranking_09252022.pdf" (10.65 X 6.10)

##  --------------
## rm "Gaussian Regression (selected features)" from "balanced accuracy" plots
data$model.v2 <- factor(data$model, levels = c("Gaussian Elastic Net (lambda=0.05)",
                                                  "Gaussian Regression (1 feature)",
                                                  #"Gaussian Regression (selected features)",
                                                  "glmnet Cox (lambda 1se)",
                                                  "glmnet Cox (lambda min)",
                                                  "Logistic Regression",
                                                  #"Multivariate Cox",
                                                  "Univariate Cox"))
data$model.v2 <- as.character(data$model.v2)

tmp.data <- data[which(data$feature.corr==0 & !is.na(data$model.v2)),]
plot.3 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Balanced.Accuracy, group=model.v2, color=model.v2)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nFeature Selection (Sensitivity + Specificity)/2\n") + 
  ggtitle("\n\nAll features are independent") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "independent_features_average_sensitivity_specificity_09252022.pdf" (10.65 X 6.10)

tmp.data <- data[which(data$feature.corr==0.8 & !is.na(data$model.v2)),]
plot.4 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Balanced.Accuracy, group=model.v2, color=model.v2)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nFeature Selection (Sensitivity + Specificity)/2\n") + 
  ggtitle("\n\nAll true features are highly correlated") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "correlated_features_average_sensitivity_specificity_09252022.pdf" (10.65 X 6.10)

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pdf("Four_panel_Summary_Plots_09252022.pdf",width=10.65*2 ,height=6.10*2)   ## Figure 1
multiplot(plotlist =list( plot.3, plot.4,plot.1, plot.2), cols=2)
dev.off()


## ------------------------------------------------------------------


tmp.data <- data[which(data$feature.corr==0 & !is.na(data$model.v2)),]
plot.5 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Sensitivity, group=model.v2, color=model.v2)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nFeature Selection Sensitivity\n") + 
  ggtitle("\n\nAll features are independent") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "independent_features_sensitivity_09252022.pdf" (10.65 X 6.10)

tmp.data <- data[which(data$feature.corr==0.8 & !is.na(data$model.v2)),]
plot.6 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Sensitivity, group=model.v2, color=model.v2)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nFeature Selection Sensitivity\n") + 
  ggtitle("\n\nAll true features are highly correlated") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "correlated_features_sensitivity_09252022.pdf" (10.65 X 6.10)

## -------------------------------------------------------

tmp.data <- data[which(data$feature.corr==0 & !is.na(data$model.v2)),]
plot.7 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Specificity, group=model.v2, color=model.v2)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nFeature Selection Specificity\n") + 
  ggtitle("\n\nAll features are independent") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "independent_features_specificity_09252022.pdf" (10.65 X 6.10)

tmp.data <- data[which(data$feature.corr==0.8 & !is.na(data$model.v2)),]
plot.8 <- ggplot(data=tmp.data, aes(x=event.sample.size.factor, y=Specificity, group=model.v2, color=model.v2)) +
  geom_line() + geom_point() + 
  xlab("\nNumber of Events\n\n") + ylab("\n\nFeature Selection Specificity\n") + 
  ggtitle("\n\nAll true features are highly correlated") + 
  #labs(subtitle = expression(paste("\n\n[Assume 5 out of 10 true features (",rho,"=0) had data available]")) ) + 
  labs(subtitle ="(Assume 5 out of 10 true features had data available)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),legend.title= element_blank())
## saved as "correlated_features_specificity_09252022.pdf" (10.65 X 6.10)


pdf("Sensitivity_Specificity_Summary_Plots_09252022.pdf",width=10.65*2 ,height=6.10*2)   ## Sup Figure 1
multiplot(plotlist =list( plot.5, plot.6,plot.7, plot.8), cols=2)
dev.off()