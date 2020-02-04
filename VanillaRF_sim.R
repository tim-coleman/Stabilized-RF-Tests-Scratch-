### Loading the workspace
load("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Pval_wkspace.RData")
library(ranger)
sim_linear_vRF <- function(p=50, n=1000, k=10, rho = 0.25, SNR = 1){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = 1:k
  beta = 1 * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n, sd = 1/SNR)
  
  result <- MSE_Test(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, NTest = 10,
                     base.learner = "rtree", ranger = T, p = 0.5)
  return(result[["Pvalue"]])
}

### Flattened sine function applied to l2norm of first k features
sim_flattened_sine_vRF <- function(p=50, n=1000, k=10, rho = 0.25, SNR = 1){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = 1:k
  flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))
  
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(n, sd = 1)
  pairs(cbind(X, y))
  plot(y~ l2_features)
  
  result <- MSE_Test(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, NTest = 10,
                     base.learner = "rtree", ranger = T, p = 0.5)
  return(result[["Pvalue"]])
}

##############################################################
# Fish covariate design
###############################################################
qsar_fish_toxicity <- read.csv("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Data/qsar_fish_toxicity.csv", header=FALSE, sep=";")
names(qsar_fish_toxicity)[7] <- "LC50"

### Function for simulating a single knock-off on a linear model
sim_linear_fish_vRF <- function(SNR = 1, k = 2, p = ncol(qsar_fish_toxicity) - 1, n = nrow(qsar_fish_toxicity)){
  X <- as.matrix(qsar_fish_toxicity[sample(nrow(qsar_fish_toxicity), n, F),-7]) 
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X)), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  X <- scale(X)
  nonzero = 1:k
  beta = 1 * (1:p %in% nonzero)
  y = X %*% beta + rnorm(nrow(X), sd = 1/SNR)
  X <- data.frame(X)
  names(X) <- paste("X", 1:p, sep = "")
  
  result <- MSE_Test(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, NTest = 10,
                     base.learner = "rtree", ranger = T, p = 0.5)
  return(result[["Pvalue"]])
}


### Flattened sine function applied to l2norm of first k features
sim_flattened_sine_fish_vRF <-  function(SNR = 1, k = 2, p = ncol(qsar_fish_toxicity) - 1, n = nrow(qsar_fish_toxicity)){
  X <- as.matrix(qsar_fish_toxicity[sample(nrow(qsar_fish_toxicity), n, F),-7])
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X)), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  X <- scale(X)
  knockoffs = function(X) create.second_order(X)
  flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(nrow(X), sd = 1/SNR)
  X <- data.frame(X)
  names(X) <- paste("X", 1:p, sep = "")
  pairs(cbind(X, y))
  
  result <- MSE_Test(X = data.frame(X), y = y, var = c("X1"), NTree = 125, B = 250, NTest = 10,
                     base.learner = "rtree", ranger = T, p = 0.5)
  return(result[["Pvalue"]])
}


run <- T
if(run){
  set.seed(367114)
  sim_params <- expand.grid(SNR = seq(0.5, 5, length.out = 10), k = seq(2, 25, length.out = 6))
  #sim_params <- expand.grid(SNR = seq(0.5, 5, length.out = 3), k = seq(2, 25, length.out = 4))
  rho <- 0.25
  p <- 25
  n <- 100
  nsim <- 2 #for educational purposes
  result_vRF <- vector(mode = "list", length = nrow(sim_params))
  for(i in 1:nrow(sim_params)){
    print(sim_params[i,])
    result_vRF[[i]][["parameters"]] <- sim_params[i,]
    suppressWarnings(result_vRF[[i]][["linear_gaussian"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with linear gaussian...\n")}; 
        sim_linear_vRF(n = n, p = p, rho = rho, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
    suppressWarnings(result_vRF[[i]][["fs_gaussian"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with flattened sine gaussian...\n")}; 
      sim_flattened_sine_vRF(n = n, p = p, rho = rho, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
    suppressWarnings(result_vRF[[i]][["linear_fish"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with linear fish...\n")}; 
      sim_linear_fish_vRF(n = n, p = p, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
    suppressWarnings(result_vRF[[i]][["fs_fish"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with flattened sine fish ...\n")}; 
      sim_flattened_sine_fish_vRF(n = n, p = p, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
  }
  save(result_vRF, file = paste("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/vRFStab_nsim=100_n=", n, "p=", p, "rho=", rho, ".RData", sep = ""))
}

run <- T
if(run){
  set.seed(367114)
  sim_params <- expand.grid(SNR = seq(0.5, 5, length.out = 5), k = seq(2, 25, length.out = 6))
  #sim_params <- expand.grid(SNR = seq(0.5, 5, length.out = 3), k = seq(2, 25, length.out = 4))
  rho <- 0.25
  p <- 50
  n <- 190
  nsim <- 10 # for educational purposes
  result_vRF <- vector(mode = "list", length = nrow(sim_params))
  for(i in 1:nrow(sim_params)){
    print(sim_params[i,])
    result_vRF[[i]][["parameters"]] <- sim_params[i,]
    suppressWarnings(result_vRF[[i]][["linear_gaussian"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with linear gaussian...\n")}; 
      sim_linear_vRF(n = n, p = p, rho = rho, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
    suppressWarnings(result_vRF[[i]][["fs_gaussian"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with flattened sine gaussian...\n")}; 
      sim_flattened_sine_vRF(n = n, p = p, rho = rho, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
    suppressWarnings(result_vRF[[i]][["linear_fish"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with linear fish...\n")}; 
      sim_linear_fish_vRF(n = n, p = p, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
    suppressWarnings(result_vRF[[i]][["fs_fish"]] <- sapply(1:nsim, FUN = function(x){ if(x == nsim/2){cat("Halfway done with flattened sine fish ...\n")}; 
      sim_flattened_sine_fish_vRF(n = n, p = p, SNR = sim_params$SNR[i], k = sim_params$k[i])}))
  }
  save(result_vRF, file = paste("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/vRFStab_nsim=100_n=", n, "p=", p, "rho=", rho, ".RData", sep = ""))
}

###################################
# Loading from cluster
###################################

cluster_load <- T
if(cluster_load){
load("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Results/vRFStab_nsim=100_n=908p=25rho=0.25.RData",
       verbose = T)
result_vRF <- result_RF
}

plot <- T
if(plot){
  library(dplyr)
  library(ggplot2)
  result_summarize <- function(obj){
    obj_no_params <- obj[-which(names(obj) == c("parameters"))]
    summaries <- function(out){
      median_pval <- median(unlist(out))
      power_10 <- mean(unlist(out) < 0.1)
      power_05 <- mean(unlist(out) < 0.05)
      mean_pval <- mean(unlist(out))
      min_pval <- min(unlist(out))
      pval_variance <- var(unlist(out))
      
      summaried <- data.frame(SNR =obj[["parameters"]][["SNR"]], k = obj[["parameters"]][["k"]], 
                              Metric = c("median_pval", "power_10", "power_05", "mean_pval", "min_pval", "pval_variance"),
                              Value = c(median_pval, power_10, power_05, mean_pval, min_pval, pval_variance))

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
  clean_results_vRF <- do.call(rbind, c(lapply(result_vRF, result_summarize), make.row.names = F))
  
  library(ggplot2)
  plot_themes <- theme_bw() + theme(strip.text = element_text(size = 13), axis.title = element_text(size = 13))
  vRF_power_k_plot <- ggplot(data = clean_results_vRF %>% filter(Metric == "power_10"), aes(x = k/25, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab("Proportion of Important Variables (s/25)") + ylab("Power") + plot_themes
  vRF_power_SNR_plot <- ggplot(data = clean_results_vRF %>% filter(Metric == "power_10"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("Power") + plot_themes
  
  vRF_pvalvar_k_plot <- ggplot(data = clean_results_vRF %>% filter(Metric == "pval_variance"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab("Proportion of Important Variables") + ylab("P-value Variance") + plot_themes
  vRF_pvalvar_SNR_plot <- ggplot(data = clean_results_vRF %>% filter(Metric == "pval_variance"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("P-value Variance") + plot_themes
  
  
  vRF_min_pval_k_plot <-  ggplot(data = clean_results_vRF %>% filter(Metric == "min_pval"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab("Proportion of Important Variables") + ylab("min_pval") + plot_themes
  vRF_min_pval_SNR_plot <- ggplot(data = clean_results_vRF %>% filter(Metric == "min_pval"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("min_pval") + plot_themes
}

save_plots <- T
if(save_plots){
  plot_sizes <- list(width = 6, height = 3.8)
  
  vRF_power_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/vRF_power_k_plot.pdf'))
  do.call(pdf, args = vRF_power_k_save)
  vRF_power_k_plot
  dev.off()
  
  vRF_pvalvar_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/vRF_pvalvar_k_plot.pdf'))
  do.call(pdf, args = vRF_pvalvar_k_save)
  vRF_pvalvar_k_plot
  dev.off()
  
  vRF_power_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/vRF_power_SNR_plot.pdf'))
  do.call(pdf, args = vRF_power_SNR_save)
  vRF_power_SNR_plot
  dev.off()
  
}



