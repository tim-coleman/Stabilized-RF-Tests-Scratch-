### Basic LASSO example
library(knockoff)
library(dplyr)

###############################################################
# Normal covariate design, should be easy for knockoffs
###############################################################

### Function for simulating a single knock-off on a linear model
sim_linear <- function(p=50, n=1000, k=10, rho = 0.35, SNR = 1, stat = stat.glmnet_coefdiff){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = 1:k
  beta = 1 * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n, sd = 1/SNR)

  knockoffs = function(X) create.gaussian(X, mu, Sigma)
  result = knockoff.filter(X, y, knockoffs=knockoffs, fdr = 0.1, statistic = stat)
  # print(result$selected)
  # print(result$threshold)
  #plot(result$statistic, type = 'h')
  #abline(h = result$threshold, col = 'lightskyblue', lwd = 1.5, lty = 'dashed')
  result <- result[c("threshold", "statistic", "selected")]
  return(result)
}

### Flattened sine function applied to l2norm of first k features
sim_flattened_sine <- function(p=50, n=1000, k=10, rho = 0.35, SNR = 1, stat = stat.glmnet_coefdiff){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  nonzero = 1:k
  flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))

  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(n, sd = 1/SNR)
  #par(mfrow = c(1,2))
  #plot(y~l2_features)
  
  knockoffs = function(X) create.gaussian(X, mu, Sigma)
  result = knockoff.filter(X, y, knockoffs=knockoffs, fdr = 0.1, statistic = stat)
  # print(result$selected)
  # print(result$threshold)
  #plot(result$statistic, type = 'h')
  #abline(h = result$threshold, col = 'lightskyblue', lwd = 1.5, lty = 'dashed')
  #par(mfrow = c(1,1))
  result <- result[c("threshold", "statistic", "selected")] 
  return(result)
}


###############################################################
# Fish covariate design, should be hard for knockoffs
###############################################################
qsar_fish_toxicity <- read.csv("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Data/qsar_fish_toxicity.csv", header=FALSE, sep=";")
names(qsar_fish_toxicity)[7] <- "LC50"

### Function for simulating a single knock-off on a linear model
sim_linear_fish <- function(SNR = 1, k = 2, stat = stat.glmnet_coefdiff, p = ncol(qsar_fish_toxicity) - 1, n = nrow(qsar_fish_toxicity)){
  X <- as.matrix(qsar_fish_toxicity[sample(nrow(qsar_fish_toxicity), n, F),-7])
  if(p - ncol(X) > 0){
  X_extra <- X[sample(nrow(X)), sample(ncol(X), p-ncol(X), replace = T)]
  X <- cbind(X, X_extra)
  }
  X <- scale(X)
  nonzero = 1:k
  beta = 1 * (1:p %in% nonzero)
  y = X %*% beta + rnorm(nrow(X), sd = 1/SNR)
  
  knockoffs = function(X) create.second_order(X)
  result = knockoff.filter(X, y, knockoffs=knockoffs, fdr = 0.1, statistic = stat)
  # print(result$selected)
  # print(result$threshold)
  #plot(result$statistic, type = 'h')
  #abline(h = result$threshold, col = 'lightskyblue', lwd = 1.5, lty = 'dashed')
  result <- result[c("threshold", "statistic", "selected")]
  return(result)
}

### Flattened sine function applied to l2norm of first k features
sim_flattened_sine_fish <-  function(SNR = 1, k = 2, stat = stat.glmnet_coefdiff, p = ncol(qsar_fish_toxicity) - 1, n = nrow(qsar_fish_toxicity)){
  X <- as.matrix(qsar_fish_toxicity[sample(nrow(qsar_fish_toxicity), n, F),-7])
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X)), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  X <- scale(X)
  #print(round(head(X, 5), 3))
  knockoffs = function(X) create.second_order(X)
  flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(nrow(X), sd = 1/SNR)
  
  #par(mfrow = c(1,2))
  #lot(y~l2_features)
  
  result = knockoff.filter(X, y, knockoffs=knockoffs, fdr = 0.1, statistic = stat)
  # print(result$selected)
  # print(result$threshold)
  #plot(result$statistic, type = 'h')
  #abline(h = result$threshold, col = 'lightskyblue', lwd = 1.5, lty = 'dashed')
  result <- result[c("threshold", "statistic", "selected")]
  return(result)
}


run <- T
if(run){
  library(pbapply)
  set.seed(367116)
  sim_params <- expand.grid(SNR = seq(0.5, 5, length.out = 10), k = seq(2, 25, length.out = 6))
  rho <- 0.25
  stat <- stat.random_forest
  p <- 25
  n <- nrow(qsar_fish_toxicity)
  nsim <- 100
  result <- vector(mode = "list", length = nrow(sim_params))
  for(i in 1:nrow(sim_params)){
    cat("Iteration:", i, "of", nrow(sim_params), "\n")
    result[[i]][["parameters"]] <- sim_params[i,]
    suppressWarnings(result[[i]][["linear_gaussian"]] <- pbsapply(1:nsim, FUN = function(x) sim_linear(n = n, p = p, rho = rho, stat = stat, SNR = sim_params$SNR[i], k = sim_params$k[i])))
    suppressWarnings(result[[i]][["fs_gaussian"]] <- pbsapply(1:nsim, FUN = function(x) sim_flattened_sine(n = n, p = p, rho = rho, stat = stat, SNR = sim_params$SNR[i], k = sim_params$k[i])))
    suppressWarnings(result[[i]][["linear_fish"]] <- pbsapply(1:nsim, FUN = function(x) sim_linear_fish(SNR = sim_params$SNR[i], stat = stat, k = sim_params$k[i], p = p)))
    suppressWarnings(result[[i]][["fs_fish"]] <- pbsapply(1:nsim, FUN = function(x) sim_flattened_sine_fish(SNR = sim_params$SNR[i], k = sim_params$k[i], p = p)))
  }
  save(result, file = paste("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/knockoff_RFstat_sim_n=", n, "p=", p, "rho=", rho, ".RData", sep = ""))
}


process_results <- T
if(process_results){
  load("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/knockoff_RFstat_sim_n=908p=25rho=0.25.RData")
  result_summarize <- function(obj){
    important <- 1:obj[["parameters"]][["k"]]; p <- length(obj[["linear_gaussian"]][,1][["statistic"]]);
    not_important <- (1:p)[-c(important)]
    obj_no_params <- obj[-which(names(obj) == c("parameters"))]
    summaries <- function(out){
      num_selected <- mean(sapply(out[3,], length))
      mean_thresh <- mean(unlist(out[1,]))
      prop_finite <- mean(is.finite(unlist(out[1,])))
      
      # We evaluate power by seeing if variable 1 is in the model or not
      # We evaluate Type I error rate by counting the # of unimportant variables selected
      power <- mean(sapply(out[3,], FUN = function(x) 1 %in% x))
      Type_I <- mean(sapply(out[3,], FUN = function(x) mean(not_important %in% x)))
      FDR <- mean(sapply(out[3,], FUN = function(x){
                                                    if(length(x) == 0){
                                                      FDR <- 0
                                                    }else{
                                                      FDR <- sum(not_important %in% x)/length(x)
                                                    }         
                                                      return(FDR)}))
                                                            
      summaried <- data.frame(SNR =obj[["parameters"]][["SNR"]], k = obj[["parameters"]][["k"]], 
                              Metric = c("num_selected", "mean_thresh", "prop_finite", "power", "Type_I", "FDR"),
                              Value = c(num_selected, mean_thresh, prop_finite, power, Type_I, FDR))
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
  clean_results <- do.call(rbind, c(lapply(result, result_summarize), make.row.names = F))
}


plot_results <- T
if(plot_results){
  library(ggplot2)
  plot_themes <- theme_bw() + theme(strip.text = element_text(size = 13), axis.title = element_text(size = 13))
  FDR_k_plot <- ggplot(data = clean_results %>% filter(Metric == "FDR"), aes(x = k/25, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("Proportion of Important Variables") + ylab("FDR") + geom_hline(aes(yintercept = 0.10), lty = 'dashed', col = 'tomato', size = 1.2) + plot_themes
  FDR_SNR_plot <- ggplot(data = clean_results %>% filter(Metric == "FDR"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("SNR") + ylab("FDR") + geom_hline(aes(yintercept = 0.10), lty = 'dashed', col = 'tomato', size = 1.2) + plot_themes
  Type_I_k_plot <- ggplot(data = clean_results %>% filter(Metric == "Type_I"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("Type I Error") + plot_themes
  Type_I_SNR_plot <- ggplot(data = clean_results %>% filter(Metric == "Type_I"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("SNR") + ylab("Type I Error") + plot_themes  
  power_k_plot <- ggplot(data = clean_results %>% filter(Metric == "power"), aes(x = k/25, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("Proportion of Important Variables (s/25)") + ylab("Power") + plot_themes
  power_SNR_plot <- ggplot(data = clean_results %>% filter(Metric == "power"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("Power") + plot_themes
  num_selected_k_plot <- ggplot(data = clean_results %>% filter(Metric == "num_selected"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("Number Variables Selected") + plot_themes
  num_selected_SNR_plot <- ggplot(data = clean_results %>% filter(Metric == "num_selected"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("SNR") + ylab("Number Variables Selected") + plot_themes
}

save_plots <- T
if(save_plots){
  plot_sizes <- list(width = 6, height = 3.8)
  
  FDR_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_FDR_k_plot.pdf'))
  do.call(pdf, args = FDR_k_save)
  FDR_k_plot
  dev.off()
  
  FDR_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_FDR_SNR_plot.pdf'))
  do.call(pdf, args = FDR_SNR_save)
  FDR_SNR_plot
  dev.off()
  
  num_selected_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_num_selected_k_plot.pdf'))
  do.call(pdf, args = num_selected_k_save)
  num_selected_k_plot
  dev.off()
  
  num_selected_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_num_selected_SNR_plot.pdf'))
  do.call(pdf, args = num_selected_SNR_save)
  num_selected_SNR_plot
  dev.off()
  
  power_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_power_k_plot.pdf'))
  do.call(pdf, args = power_k_save)
  power_k_plot
  dev.off()
  
  power_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_power_SNR_plot.pdf'))
  do.call(pdf, args = power_SNR_save)
  power_SNR_plot
  dev.off()
  
  Type_I_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_Type_I_k_plot.pdf'))
  do.call(pdf, args = Type_I_k_save)
  Type_I_k_plot
  dev.off()
  
  Type_I_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/Knockoff_Type_I_SNR_plot.pdf'))
  do.call(pdf, args = Type_I_SNR_save)
  Type_I_SNR_plot
  dev.off()
  }

