### Loading the workspace
load("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Pvalue_combine.RData")
library(ranger)
sim_linear_RF <- function(p=50, n=1000, k=10, rho = 0.25, SNR = 1){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = 1:k
  beta = 1 * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n, sd = 1/SNR)

  result <- stable_MSE_Test.default(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250,
                                    n_val = 50, base.learner = "rtree", ranger = T, p = 0.5)
  print(result[["Agg_Pvalue"]])
  return(result)
}

### Flattened sine function applied to l2norm of first k features
sim_flattened_sine_RF <- function(p=50, n=1000, k=10, rho = 0.25, SNR = 1){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = 1:k
  flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))
  
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(n, sd = 1)

  result <- stable_MSE_Test.default(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, 
                                    n_val = 50, base.learner = "rtree", ranger = T, p = 0.5)
  print(result[["Agg_Pvalue"]])
  return(result)
}

###############################################################
# Fish covariate design
###############################################################
qsar_fish_toxicity <- read.csv("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Data/qsar_fish_toxicity.csv", header=FALSE, sep=";")
names(qsar_fish_toxicity)[7] <- "LC50"

### Function for simulating a single knock-off on a linear model
sim_linear_fish_RF <-  function(SNR = 1, k = 2, p = ncol(qsar_fish_toxicity) - 1, n = nrow(qsar_fish_toxicity)){
  X <- as.matrix(qsar_fish_toxicity[sample(nrow(qsar_fish_toxicity), n, F),-7])
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X)), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  nonzero = 1:k
  beta = 1 * (1:p %in% nonzero)
  y = X %*% beta + rnorm(nrow(X), sd = 1/SNR)
  X <- data.frame(X)
  names(X) <- paste("X", 1:p, sep = "")

  result <- stable_MSE_Test.default(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, 
                                    n_val = 50, base.learner = "rtree", ranger = T, p = 0.5)
  print(result[["Agg_Pvalue"]])
  return(result)
}


### Flattened sine function applied to l2norm of first k features
sim_flattened_sine_fish_RF <- function(SNR = 1, k = 2, p = ncol(qsar_fish_toxicity) - 1, n = nrow(qsar_fish_toxicity)){
  X <- as.matrix(qsar_fish_toxicity[sample(nrow(qsar_fish_toxicity), n, F),-7])
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X)), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  #print(round(head(X, 5), 3))
  knockoffs = function(X) create.second_order(X)
  flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(nrow(X), sd = 1/SNR)
  X <- data.frame(X)
  names(X) <- paste("X", 1:p, sep = "")
  
  result <- stable_MSE_Test.default(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, 
                                    n_val = 50, base.learner = "rtree", ranger = T, p = 0.5)
  print(result[["Agg_Pvalue"]])
  return(result)
}


run <- F
if(run){
  set.seed(367114)
  sim_params <- expand.grid(SNR = seq(0.5, 5, length.out = 10), k = seq(2, 25, length.out = 6))
  rho <- 0.25
  p <- 25
  n <- nrow(qsar_fish_toxicity)
  nsim <- 2 # for educational purposes
  result_RF <- vector(mode = "list", length = nrow(sim_params))
  for(i in 1:nrow(sim_params)){
    result_RF[[i]][["parameters"]] <- sim_params[i,]
    suppressWarnings(result_RF[[i]][["linear_gaussian"]] <- sapply(1:nsim, FUN = function(x) sim_linear_RF(n = n, p = p, rho = rho, SNR = sim_params$SNR[i], k = sim_params$k[i])))
    suppressWarnings(result_RF[[i]][["fs_gaussian"]] <- sapply(1:nsim, FUN = function(x) sim_flattened_sine_RF(n = n, p = p, rho = rho, SNR = sim_params$SNR[i], k = sim_params$k[i])))
    suppressWarnings(result_RF[[i]][["linear_fish"]] <- sapply(1:nsim, FUN = function(x) sim_linear_fish_RF(SNR = sim_params$SNR[i], k = sim_params$k[i], p = p)))
    suppressWarnings(result_RF[[i]][["fs_fish"]] <- sapply(1:nsim, FUN = function(x) sim_flattened_sine_fish_RF(SNR = sim_params$SNR[i], k = sim_params$k[i], p = p)))
  }
  #save(result_RF, file = paste("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/RFStab_nsim=2_n=", n, "p=", p, "rho=", rho, ".RData", sep = ""))
}


plot <- T
if(plot){
  library(dplyr)
  result_summarize <- function(obj){
    obj_no_params <- obj[-which(names(obj) == c("parameters"))]
    summaries <- function(out){
      median_pval <- median(unlist(out[1,]))
      power_10 <- mean(unlist(out[1,]) < 0.1)
      power_05 <- mean(unlist(out[1,]) < 0.05)
      mean_pval <- mean(unlist(out[1,]))
      min_pval <- min(unlist(out[1,]))
      
      summaried <- data.frame(SNR =obj[["parameters"]][["SNR"]], k = obj[["parameters"]][["k"]], 
                              Metric = c("median_pval", "power_10", "power_05", "mean_pval", "min_pval"),
                              Value = c(median_pval, power_10, power_05, mean_pval, min_pval))
      return(summaried)
    }
    
    sim_summaries <- lapply(obj_no_params, summaries)
    for(i in 1:length(sim_summaries)){
      sim_type <- names(sim_summaries)[i]
      split_str <- unlist(strsplit(sim_type, split = "_"))
      model_Y <- split_str[1]
      model_X <- split_str[2]
      sim_summaries[[i]][["model_X"]] <- model_X
      sim_summaries[[i]][["model_Y"]] <- model_Y
    }
    out <- do.call(rbind, c(sim_summaries, make.row.names = F))
    return(out)
  }
  clean_results_RF <- do.call(rbind, c(lapply(result_RF, result_summarize), make.row.names = F))
  
  library(ggplot2)
  plot_themes <- theme_bw() + theme(strip.text = element_text(size = 13), axis.title = element_text(size = 13))
  power_k_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "power_10"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("Power") + plot_themes
  power_SNR_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "power_10"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("Power") + plot_themes
  
  min_pval_k_plot <-  ggplot(data = clean_results_RF %>% filter(Metric == "min_pval"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("min_pval") + plot_themes
  min_pval_SNR_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "min_pval"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("min_pval") + plot_themes
}


