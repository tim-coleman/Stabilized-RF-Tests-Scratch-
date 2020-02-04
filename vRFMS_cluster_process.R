#######################################
# Loading the cluster results 
#######################################

### null case
load("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Results/vRF_Multisplit_h0_nsim=100_n=908p=25rho=0.25.RData")
result_RF_h0 <- result_RF
rm(result_RF)
library(dplyr)
library(ggplot2)
######################################
# Processing
######################################

plot_H0 <- T
if(plot_H0){
  library(dplyr)
  result_summarize <- function(obj){
    obj_no_params <- obj[-which(names(obj) == c("parameters"))]
    summaries <- function(out){
      median_pval <- median(unlist(out[1,]))
      power_10 <- mean(unlist(out[1,]) < 0.1)
      power_05 <- mean(unlist(out[1,]) < 0.05)
      mean_pval <- mean(unlist(out[1,]))
      min_pval <- min(unlist(out[1,]))
      pval_variance <- var(unlist(out[1,]))
      
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
  
  clean_results_RF_h0 <- do.call(rbind, c(lapply(result_RF_h0, result_summarize), make.row.names = F))
  
  library(ggplot2)
  plot_themes <- theme_bw() + theme(strip.text = element_text(size = 13), axis.title = element_text(size = 13))
  RF_power_k_plot <- ggplot(data = clean_results_RF_h0%>% filter(Metric == "power_10"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("Power") + plot_themes
  RF_power_SNR_plot <- ggplot(data = clean_results_RF_h0 %>% filter(Metric == "power_10"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("Power") + plot_themes
  
  RF_pvalvar_k_plot <- ggplot(data = clean_results_RF_h0 %>% filter(Metric == "pval_variance"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("P-value Variance") + plot_themes
  RF_pvalvar_SNR_plot <- ggplot(data = clean_results_RF_h0 %>% filter(Metric == "pval_variance"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("P-value Variance") + plot_themes
  
  RF_min_pval_k_plot <-  ggplot(data = clean_results_RF_h0 %>% filter(Metric == "min_pval"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("min_pval") + plot_themes
  RF_min_pval_SNR_plot <- ggplot(data = clean_results_RF_h0 %>% filter(Metric == "min_pval"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("min_pval") + plot_themes

  plot_sizes <- list(width = 7, height = 4.5)
  
  RF_power_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/vRFMS_h0_k_plot.pdf'))
  do.call(pdf, args = RF_power_k_save)
  RF_power_k_plot
  dev.off()
  
  RF_power_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/vRFMS_h0_SNR_plot.pdf'))
  do.call(pdf, args = RF_power_SNR_save)
  RF_power_SNR_plot
  dev.off()
  
}


#######################################
# Loading the high dimensional cluster results 
#######################################
load("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Results/RFStab_highdimsim_nsim100_n=100p=100rho=0.25.RData")
library(dplyr)
library(ggplot2)
######################################
# Processing
######################################

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
      pval_variance <- var(unlist(out[1,]))
      power_var <- var(unlist(out[1,] < 0.05))
      
      summaried <- data.frame(SNR =obj[["parameters"]][["SNR"]], k = obj[["parameters"]][["k"]], 
                              Metric = c("median_pval", "power_10", "power_05", "mean_pval", "min_pval", "pval_variance", "power_var"),
                              Value = c(median_pval, power_10, power_05, mean_pval, min_pval, pval_variance, power_var))
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
  RF_power_k_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "power_10"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("Power") + plot_themes
  RF_power_SNR_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "power_10"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("Power") + plot_themes
  
  RF_pvalvar_k_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "pval_variance"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("P-value Variance") + plot_themes
  RF_pvalvar_SNR_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "pval_variance"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("P-value Variance") + plot_themes
  
  RF_min_pval_k_plot <-  ggplot(data = clean_results_RF %>% filter(Metric == "min_pval"), aes(x = k, y = Value, group = SNR)) + 
    geom_line(aes(col = SNR)) + geom_point(aes(col = SNR)) + facet_wrap(c("model_X", "model_Y")) + 
    xlab ("k") + ylab("min_pval") + plot_themes
  RF_min_pval_SNR_plot <- ggplot(data = clean_results_RF %>% filter(Metric == "min_pval"), aes(x = SNR, y = Value, group = k)) + 
    geom_line(aes(col = k)) + geom_point(aes(col = k)) + facet_wrap(c("model_X", "model_Y")) +
    xlab ("SNR") + ylab("min_pval") + plot_themes
}

save_plots <- T
if(save_plots){
  plot_sizes <- list(width = 7, height = 4.5)
  
  RF_power_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/RF_power_k_plot.pdf'))
  do.call(pdf, args = RF_power_k_save)
  RF_power_k_plot
  dev.off()
  
  RF_pvalvar_k_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/RF_pvalvar_k_plot.pdf'))
  do.call(pdf, args = RF_pvalvar_k_save)
  RF_pvalvar_k_plot
  dev.off()
  
  RF_power_SNR_save <- c(plot_sizes, file = c('C:/Users/drain/Box Sync/MeSntch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Plots/RF_power_SNR_plot.pdf'))
  do.call(pdf, args = RF_power_SNR_save)
  RF_power_SNR_plot
  dev.off()
  
}

