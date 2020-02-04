### This is an adaptation of the P-value procedure from Meinshausen et al. (2008)
### "p-Values for High-Dimensional Regression"

source('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/R Stuff/MSE_Test_File.R')

cluster <- T
if(!cluster){
  # using pbsapply for debugging purposes locally, changes to pbsapply
  library(pbapply)
  sapply <- pbsapply
  pboptions(char = "=")
}

bag.s <- function(X, y, base.learner = "rpart", ntree, k, mtry = ncol(X),
                  form = as.formula("y~."), return_ind = F,
                  alpha = if(base.learner == "lm") 1 else NULL ,
                  glm_cv = if(base.learner == "lm") "external" else "none", 
                  lambda = if(glm_cv == "none" & base.learner == "lm") 1 else NULL,
                  ranger = F){
  
  y <- unlist(y) # vectorizing y
  #resp_name <- strsplit(as.character(form), split = "")[[1]][1]
  resp_name <- as.character(form)[2]
  D <- data.frame(X, "y" = y)
  names(D)[ncol(X) + 1] <- resp_name
  
  if(!(base.learner %in% c("rpart", "ctree", "rtree", "lm"))){
    stop("Base learner should be one of 'rpart', 'ctree', 'lm' or 'rtree' ")
  }
  
  if(nrow(D) == 0) stop("data (D) has 0 rows")
  
  if(base.learner == "rpart"){
    c.sim <- rpart.control(minsplit=10, maxcompete=0, maxsurrogate=0, usesurrogate=0, cp = 0)
    fun = function(){
      ind <- sample(1:nrow(D), size=k, replace=FALSE)
      return(list("model" = rpart(form, data = D[ind,], control = c.sim),
                  "ind" = ind))
    }
  }
  
  if(base.learner == "ctree"){
    if(ranger){
      fun = function(){
        ind <- sample(1:nrow(D), size=k, replace=FALSE)
        return(list("model" = ranger(form, num.trees = 1, mtry = mtry, importance = 'none', data = D[ind,],
                                     min.node.size = 1, replace = F, sample.fraction = 1, splitrule = "maxstat"),
                    "ind" = ind))
      }
    }
    else{
      c.sim <- ctree_control(minsplit=1,
                             maxsurrogate=0, mtry = mtry,
                             mincriterion = .3)
      fun = function(){
        ind <- sample(1:nrow(D), size=k, replace=FALSE)
        return(list("model" = ctree(form, controls = c.sim, data = D[ind,]), 
                    "ind" = ind))
      }
    }
  }
  
  
  if(base.learner == "rtree") {
    if(ranger){
      fun =  function(){        
        ind <- sample(1:nrow(D), size=k, replace=FALSE)
        return(list("model" = ranger(form, num.trees = 1, mtry = mtry, importance = 'none', data = D[ind,],
                                     min.node.size = 1, replace = F, sample.fraction = 1),
                    "ind" = ind))
      }
    }
    else{
      fun = function(){
        ind <- sample(1:nrow(D), size = k, replace = FALSE)
        return(list("model" = randomForest(x = X[ind,], y = y[ind], ntree = 1, mtry = mtry,
                                           replace = FALSE,
                                           sampsize = floor(k)), "ind" = ind))
      }
    }
    
  }
  
  if(base.learner == "lm"){
    if(glm_cv == "external"){
      cv0 <- cv.glmnet(x = model.matrix(form, data = D), y = y, alpha = alpha)
      l_cv <- cv0[["lambda.1se"]]
      fun = function(){
        ind <- sample(1:nrow(D), size=k, replace=FALSE)
        features <- model.matrix(form, data = D[ind,])
        return(list("model" = glmnet(x = features, y = y[ind], alpha = alpha, lambda = l_cv),
                    "ind" = ind))
      }
    }
    else{
      fun = function(){
        ind <- sample(1:nrow(D), size=k, replace=FALSE)
        features <- model.matrix(form, data = D[ind,])
        if(glm_cv == "none"){
          return(list("model" = glmnet(x = features, y = y[ind], alpha = alpha, lambda = lambda),
                      "ind" = ind))
        }
        else{
          return(list("model" = cv.glmnet(x = features, y = y[ind], alpha = alpha, nfolds = 5), 
                      "ind" = ind))
        }
      }
    }
  }
  
  rfs <- lapply(1:ntree, function(x) fun())
  
  if(!return_ind){
    rfs <- lapply(rfs, function(x) x[["model"]])
  }
  
  if(ranger & base.learner %in% c("rtree", "ctree")){
    class(rfs) <- paste("bagged", base.learner, "ranger", sep = "_")
    
  }
  else{
    class(rfs) <- paste("bagged", base.learner, sep = "_")
  }
  rfs
}

MultiSplit_MSE_Test.default <- function(X, y, var, NTest, B = 100, NTree = 500, p = 1/2, n_val = nrow(X), 
                                    base.learner = "rpart", mtry = ncol(X), importance = T, 
                                    alpha = if(base.learner == "lm") 1,
                                    glm_cv = if(base.learner == "lm") "external" else "none", 
                                    lambda = if(glm_cv == "none" & base.learner == "lm") 1 else NULL, ranger = F){
  

  Pvals <- sapply(1:n_val, function(x) {
    out <- MSE_Test(X = X, y = y, NTree = NTree, B = B, p = p, NTest = NTest, var = c("X3"), base.learner = base.learner,
                    importance = importance, alpha = alpha, glm_cv = glm_cv, lambda = lambda, ranger = ranger)
    return(out[["Pvalue"]])
  })
  
  gamma_min <- 0.05
  Q_g <- function(g) min(c(1, quantile(Pvals/g, probs = g)))
  p_out <- min(1, (1-log(gamma_min)) * min(sapply(seq(gamma_min, 1, 0.01), FUN = function(g) Q_g(g))))

  return(list(Agg_Pvalue = p_out, Pvalues = Pvals))
}

### Loading in some necessary data
qsar_fish_toxicity <- read.csv("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Data/qsar_fish_toxicity.csv", header=FALSE, sep=";")
names(qsar_fish_toxicity)[7] <- "LC50"

### A nonlinear response function
flattened_sine <- function(x) sin(pi*(x-mean(x)))/((x-mean(x)))


### Simulation functions

sim_linear_RF <- function(p=50, n=1000, k=10, rho = 0.25, SNR = 1, H_0 = F){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  beta = 1 * (1:p %in% 1:k)
  y = X %*% beta + rnorm(n, sd = 1/SNR)
  X <- scale(X)
  X <- data.frame(X)
  if(H_0){
    X[["X1"]] <- sample(X[["X1"]], nrow(X), F)
    result <- MultiSplit_MSE_Test.default(X = X,  y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                      n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  } else {
    result <- MultiSplit_MSE_Test.default(X = X, y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                      n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  }
  return(result)
}

sim_flattened_sine_RF <- function(p=50, n=1000, k=10, rho = 0.25, SNR = 1, H_0 = F){
  mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  X <- scale(X)
  nonzero = 1:k
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(n, sd = 1/SNR)
  X <- data.frame(X)
  if(H_0){
    X[["X1"]] <- sample(X[["X1"]], nrow(X), F)
    result <- MultiSplit_MSE_Test.default(X = X,  y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                          n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  } else{
    result <- MultiSplit_MSE_Test.default(X = X, y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                          n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  }
  return(result)
}

sim_linear_fish_RF <- function(n=1000, k=10, rho = 0.25, SNR = 1,  H_0 = F, p = ncol(qsar_fish_toxicity) - 1){
  X <- as.matrix(qsar_fish_toxicity[,-7])
  if(n > nrow(X)){
    X <- X[sample(nrow(X), n, F),]
  } else{
    X <- X[sample(nrow(X), n, T),]
  }
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X), n, replace = T), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  X <- scale(X)
  beta = 1 * (1:p %in% 1:k)
  y = X %*% beta + rnorm(nrow(X), sd = 1/SNR)
  X <- data.frame(X)
  names(X) <- paste("X", 1:p, sep = "")
  
  if(H_0){
    X[["X1"]] <- sample(X[["X1"]], nrow(X), F)
    result <- MultiSplit_MSE_Test.default(X = X,  y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                          n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  } else {
    result <- MultiSplit_MSE_Test.default(X = X, y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                          n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  }
  return(result)
}

sim_flattened_sine_fish_RF <- function(n=1000, k=10, rho = 0.25, SNR = 1,  H_0 = F, p = ncol(qsar_fish_toxicity) - 1){
  X <- as.matrix(qsar_fish_toxicity[,-7])
  if(n > nrow(X)){
    X <- X[sample(nrow(X), n, T),]
  } else {
    X <- X[sample(nrow(X), n, T),]
  }
  if(p - ncol(X) > 0){
    X_extra <- X[sample(nrow(X), n, replace = T), sample(ncol(X), p-ncol(X), replace = T)]
    X <- cbind(X, X_extra)
  }
  X <- scale(X)
  l2_features <- apply(X[,1:k], FUN = function(x) sqrt(sum(x^2)), MARGIN = 1)
  y <- flattened_sine(l2_features) + rnorm(nrow(X), sd = 1/SNR)
  X <- data.frame(X)
  names(X) <- paste("X", 1:p, sep = "")
  if(H_0){
    X[["X1"]] <- sample(X[["X1"]], nrow(X), F)
    result <- MultiSplit_MSE_Test.default(X = X,  y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                          n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  } else{
    result <- MultiSplit_MSE_Test.default(X = X, y = y, var = c("X1"), NTree = 125, B = 250, NTest = n/10,
                                          n_val = 50, base.learner = "rtree", ranger = T, p = 0.5, importance = T)
  }
  return(result)
}


function_validate <- F
if(function_validate){
  set.seed(367116)
  SNR <- 3.5
  n <- 1000
  validated <- list(
    fs_g_h0 = sim_flattened_sine_RF(p = 10, n = n, SNR = SNR, k = 6, H_0 = T),
    fs_g_h1 = sim_flattened_sine_RF(p = 10, n = n, SNR = SNR, k = 6,  H_0 = F),
    l_g_h0 = sim_linear_RF(p = 10, n = n, SNR = 2.5, k = 6,  H_0 = T),
    l_g_h1 = sim_linear_RF(p = 10, n = n, SNR = 2.5, k = 6, H_0 = F),
    fs_f_h0 = sim_flattened_sine_fish_RF(p = 10, n = n, k = 6,  SNR = SNR, H_0 = T),
    fs_f_h1 = sim_flattened_sine_fish_RF(p = 10, n = n, k = 6, SNR = SNR, H_0 = F),
    l_f_h0 = sim_linear_fish_RF(p = 10, n = n, SNR = SNR, k= 6, H_0 = T),
    l_f_h1 = sim_linear_fish_RF(p = 10, n = n, SNR = SNR, k = 6,  H_0 = F)
  )
  lapply(validated, FUN = function(x) x[["Agg_Pvalue"]])

}

cluster <- T
if((1-function_validate)*cluster){
save.image("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/MultiSplit_Pval_wkspace.RData")
}