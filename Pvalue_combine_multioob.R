### This is an adaptation of the P-value procedure from Meinshausen et al. (2008)
### "p-Values for High-Dimensional Regression"

source('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/R Stuff/MSE_Test_File.R')

N <- 500
Nvar <- 5
N_test <- 50
name_vec <- paste("X", 1:(2*Nvar), sep = "")

# training data:
X.design <- data.frame(replicate(Nvar, runif(N)),
                       replicate(Nvar, cut(runif(N), 3,
                                           labels = as.character(1:3))))

X <- X.design %>% mutate(Y = sqrt(X2) + rnorm(N, sd = .45))

names(X) <- c(name_vec, "Y")




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



stable_MSE_Test.default <- function(X, y, var, B = 1000, NTree = 500, p = 1/2, n_val = nrow(X), oob_size = 1,
                                    base.learner = "rpart", mtry = ncol(X), importance = T, 
                                    alpha = if(base.learner == "lm") 1,
                                    glm_cv = if(base.learner == "lm") "external" else "none", 
                                    lambda = if(glm_cv == "none" & base.learner == "lm") 1 else NULL, ranger = F){
  
  if(mtry > ncol(X)){
    warning("mtry is greater than number of columns in design matrix \n setting mtry = ncol(X)")
    mtry <- ncol(X)
  }
  y <- unlist(y)
  N <- nrow(X)
  X.train <- X
  X.train.pm <- X
  y.train <- y
  X.train.pm[,var] <- X.train.pm[sample(nrow(X.train)), var]
  
  
  #form.resp <- as.formula(paste(resp, "~.", sep = ""))
  form.resp <- y~.
  
  rf_og <- bag.s(X = X.train, y = y.train, base.learner = base.learner, return_ind = T,
                 ntree = NTree, k = ceiling(N^p), mtry = mtry, form = form.resp, 
                 alpha = alpha, lambda = lambda, ranger = ranger, glm_cv = glm_cv)
  rf_pm <- bag.s(X = X.train.pm, y = y.train, base.learner = base.learner, return_ind = T,
                 ntree = NTree, k = ceiling(N^p), mtry = mtry, form = form.resp,
                 alpha = alpha, lambda = lambda, ranger = ranger, glm_cv = glm_cv)
  
  cat("Finished model training ...\n")
  
  # Model to record out of bag predictions for each point
  oob_pred_f <- function(model, ind){
    model <- lapply(model, function(x){
      if(!any(ind %in% x[["ind"]])){
        if(class(model) %in% c("bagged_rtree", "bagged_ctree", "bagged_rpart")){
          P <- unlist(lapply(list(x[["model"]]), FUN = function(x) predict(x, newdata = X.train[ind,])))
        }
        else if(class(model) %in% c("bagged_rtree_ranger","bagged_ctree_ranger")){
          P <- unlist(lapply(list(x[["model"]]), FUN = function(x){ P <- predict(x, data = X.train[ind,]); P["predictions"] }))
        }
        else{
          form.resp <- as.formula("y~.")
          P <- predict(x[["model"]], newx = model.matrix(form.resp, data.frame(X.train[ind,], "y" = y.train[ind])))
        }
        return(P)
      }
      else{
        return(NULL)
      }
    })
    return(do.call(rbind, model[lengths(model) > 0]))
  }
  
  # Permutation test for each collection of oob points
  oob_perm_f <- function(ind){
    P <- oob_pred_f(rf_og, ind) 
    PR <- oob_pred_f(rf_pm, ind)
    
    num_tree_0 <- nrow(P)
    num_tree_R <- nrow(PR)
    forest_size <- min(num_tree_0, num_tree_R)
    #(c("B_0" = num_tree_0, "B_R" = num_tree_R))
    if(forest_size < 10){ # Really want at least 10 trees in each forest
      return(NULL)
    }
    # Scaling down the larger of the two forests, and predicting
    Pred_0 <- apply(P, FUN = mean, MARGIN = 2)
    Pred_R_0 <- apply(PR, FUN = mean, MARGIN = 2)
    
    
    MSE_0 <- mean((Pred_0 - y[ind])^2)
    MSE_R_0 <- mean((Pred_R_0 - y[ind])^2)
    diff.0 <- MSE_R_0 - MSE_0
    Pool <- rbind(P, PR)
    MSE <- data.frame("Full_MSE" = rep(B, 0), "Reduced_MSE" = rep(B, 0)) # Intialization
    for(i in 1:B){
      samps <- sample(1:(num_tree_0 + num_tree_R), num_tree_0, replace = F)
      Pred_t <- mean(Pool[samps,])
      Pred_R_t <- mean(Pool[-samps,])
      MSE[i,] <-c(mean((Pred_t - y[ind])^2), mean((Pred_R_t - y[ind])^2))
    }
    MSE.diffs <- MSE$Reduced_MSE - MSE$Full_MSE
    Pval <- c("Pval" = mean(c(diff.0 < MSE.diffs, 1)))
    if(importance){
      sdimp <- (diff.0 - mean(MSE.diffs))/sd(MSE.diffs)
      zimp <- pnorm(sdimp)
      importances <- c("Standard Deviation Importance" = sdimp,
                       "Standard Normal Importance" = zimp)
    }
    else{ 
      importances <- NULL
    }
    return(c(Pval, importances, forest_size = forest_size))
  }
  
  # Running the prediction stage
  # Sampling multiple sets without replacement
  # available_inds <- 1:nrow(X)
  # validation_indices <- sample(available_inds, oob_size, replace = F)
  # available_inds <- available_inds[-available_inds]
  # for(i in 1:n_val){
  #   temp <- sample(available_inds, oob_size, replace = F)
  #   validation_indices <- cbind(validation_indices, temp)
  #   available_inds
  # }
  # validation_indices <- replicate(n_val, sample(nrow(X), oob_size, replace = F), simplify = F)
  if(oob_size * n_val > nrow(X)){
    stop("oob_size * n_val <= nrow(X) is required")
  }
  ind_pool <- sample(nrow(X), oob_size*n_val, replace = F)
  validation_indices <- split(ind_pool, ceiling(seq_along(ind_pool)/oob_size), drop = T)
  out <- do.call(rbind, lapply(validation_indices, oob_perm_f))

  gamma_min <- 0.05
  Q_g <- function(g) min(c(1, quantile(out[,1]/g, probs = g)))
  p_out <- min(1, (1-log(gamma_min)) * min(sapply(seq(gamma_min, 1, 0.01), FUN = function(g) Q_g(g))))
  
  hist(out[,1], col = 'gray', breaks = 9)
  #print(mean(out[,1] < 0.05))
  #print(p_out)
  cat("Mean Forest Size:", mean(data.frame(out)[["forest_size"]]), "\n")
  
    return(list(Agg_Pvalue = p_out, Pvalues = out[,1], forest_size = mean(data.frame(out)$forest_size)))
}

N <- 1500
Nvar <- 5
N_test <- 50
sd <- 0.25
name_vec <- paste("X", 1:(2*Nvar), sep = "")

# training data:
X.design <- data.frame(replicate(Nvar, runif(N)),
                       replicate(Nvar, cut(runif(N), 3,
                                           labels = as.character(1:3))))

X <- X.design %>% mutate(Y = sqrt(X2) + rnorm(N, sd = sd))

SNR <- var(sqrt(X$X2))/(sd^2)
print(SNR)
names(X) <- c(name_vec, "Y")


H0_stable_10 <- stable_MSE_Test.default(X = X.design, y = X$Y, var = c("X1"), oob_size = 2,
                                                n_val = 75, base.learner = "rtree", NTree = 100, p = .5, mtry = 5, B = 100, ranger = T)
print(H0_stable_10[["Agg_Pvalue"]])
H1_stable_10 <- stable_MSE_Test.default(X = X.design, y = X$Y, var = c("X2"), oob_size = 2,
                                        n_val = 75, base.learner = "rtree", NTree = 100, p = .5, mtry = 5, B = 100, ranger = F)
print(H1_stable_10[["Agg_Pvalue"]])





#### Validating against a multisplit aggregation
marginal_combine <- function(n_val, nt, p){
  
  N <- 1500
  Nvar <- 5
  N_test <- 50
  sd <- 0.25
  name_vec <- paste("X", 1:(2*Nvar), sep = "")
  
  # training data:
  X.design <- data.frame(replicate(Nvar, runif(N)),
                         replicate(Nvar, cut(runif(N), 3,
                                             labels = as.character(1:3))))
  
  X <- X.design %>% mutate(Y = sqrt(X2) + rnorm(N, sd = sd))
  
  SNR <- var(sqrt(X$X2))/(sd^2)
  print(SNR)
  names(X) <- c(name_vec, "Y")
  
  single_MSE_run <- function(sim){
    H0<- MSE_Test.default(X = X.design, y = X$Y, var = c("X1"), NTest = nt,
                          base.learner = "rtree", NTree = 150, p = p, mtry = 5, B = 100, ranger = F)
    H1<- MSE_Test.default(X = X.design, y = X$Y, var = c("X2"), NTest = nt,
                          base.learner = "rtree", NTree = 150, p = p, mtry = 5, B = 100, ranger = F)
    if(n_val %% sim == 10){
    print(out <- c("Null" = H0[["Pvalue"]], "Alternative" = H1[["Pvalue"]]))
    }
    else{
      out <- c("Null" = H0[["Pvalue"]], "Alternative" = H1[["Pvalue"]])
    }
    return(out)
  }
  out <- sapply(1:n_val, single_MSE_run)
  par(mfrow = c(1,2))
  hist(out[1,], main = "Null P-value Distribution", col = 'gray')
  hist(out[2,], main = "Alternative P-value Distribution", col = 'gray')
  par(mfrow = c(1,1))
  aggregation <- function(x){
    gamma_min <- 0.05
    Q_g <- function(g) min(c(1, quantile(x/g, probs = g)))
    p_out <- min(1, (1-log(gamma_min)) * min(sapply(seq(gamma_min, 1, 0.01), FUN = function(g) Q_g(g))))
  }
  agg_H0 <- aggregation(out[1,])
  agg_H1 <- aggregation(out[2,])
  cat("Aggregated P-values:\n")
  print(aggregated <-  c("Null" = agg_H0, "Alternative" = agg_H1))
  return(list(out, aggregated))
}

marginal_combine(n_val = 100, nt = 50, p  = 0.5)
