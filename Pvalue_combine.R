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
  
X <- X.design %>% mutate(Y = 5*(X3) + 4*X2^2 + rnorm(N, sd = .45))

names(X) <- c(name_vec, "Y")

multi_split_combine <- function(nsplit = 10, gamma_min = 0.075){
  Pvals <- sapply(1:nsplit, function(x) {
    out <- MSE_Test(X = X.design, y = X$Y, NTree = 100, B = 500, NTest = N_test, var = c("X3"))
                  return(out$Pvalue)
    })
  hist(Pvals, col = 'gray', breaks = 30)
  Q_g <- function(g) min(c(1, quantile(Pvals/g, probs = g)))
  p_out <- min(1, (1-log(gamma_min)) * min(sapply(seq(gamma_min, 1, 0.01), FUN = function(g) Q_g(g))))
  print(p_out)
  return(Pvals)
}

# #mspc <- multi_split_combine(50)
# 
# gamma_min <- 0.075
# Pvals <- rbeta(100, shape1 = 1, shape2 = 25)
# min(1, (1-log(gamma_min)) * min(sapply(seq(gamma_min, 1, 0.01), FUN = function(g) Q_g(g))))
# plot(seq(0, 1, 0.001), sapply(seq(0, 1, .001), FUN = Q_g), type = 'l')
# 



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


### Testing this out
b.rpart <- bag.s(X = data.frame(replicate(5, rnorm(250))), y = rnorm(250),
                 base.learner = "rpart", ntree = 10, k = 50, mtry = 5, form = Y~.)
b.rf.inds <- bag.s(X = data.frame(replicate(5, rnorm(250))), y = rnorm(250), return_ind = T,
                 base.learner = "rtree", ntree = 10, k = 50, mtry = 5, form = Y~.)
b.rf <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
              base.learner = "rtree", ntree = 10, k = N^.95, mtry = 2, Y~., ranger = F)
b.glmnet <-  bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                   base.learner = "lm", ntree = 10, k = N^.95, mtry = 2)


stable_MSE_Test.default <- function(X, y, var, B = 1000, NTree = 500, p = 1/2, n_val = nrow(X), base.learner = "rpart",
                             mtry = ncol(X), importance = T, alpha = if(base.learner == "lm") 1,
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
      if(!(ind %in% x[["ind"]])){
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
    return(unlist(model[lengths(model) > 0]))
  }
  
  # One by One holdout of each point
  oob_perm_f <- function(ind){
    P <- oob_pred_f(rf_og, ind) 
    PR <- oob_pred_f(rf_pm, ind)
    
    num_tree_0 <- length(P)
    num_tree_R <- length(PR)
    
    Pred_0 <- mean(P)
    Pred_R_0 <- mean(PR)
    #print(c("B_0" = num_tree_0, "B_R" = num_tree_R))
    
    MSE_0 <- mean((Pred_0 - y[ind])^2)
    MSE_R_0 <- mean((Pred_R_0 - y[ind])^2)
    diff.0 <- MSE_R_0 - MSE_0
    Pool <- c(P, PR)
    MSE <- data.frame("Full_MSE" = rep(B, 0), "Reduced_MSE" = rep(B, 0)) # Intialization
    for(i in 1:B){
      samps <- sample(1:(num_tree_0 + num_tree_R), num_tree_0, replace = F)
      P_t <- Pool[samps]
      PR_t <- Pool[-samps]
      Pred_t <- mean(P_t)
      Pred_R_t <- mean(PR_t)
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
      importances = NULL
    }
    return(c(Pval, importances))
  }
  
  # Running the prediction stage
  out <- do.call(rbind, lapply(1:n_val, oob_perm_f))
  print(head(out))
  
  gamma_min <- 0.05
  Q_g <- function(g) min(c(1, quantile(out[,1]/g, probs = g)))
  p_out <- min(1, (1-log(gamma_min)) * min(sapply(seq(gamma_min, 1, 0.01), FUN = function(g) Q_g(g))))
  
  hist(out[,1], col = 'gray', breaks = 7)
  #print(mean(out[,1] < 0.05))
  #print(p_out)

  
  # result <- list("variables" = var,
  #                "originalStat" = c("Original MSE" = MSE_0, 
  #                                   "Permuted MSE" = MSE_R_0),
  #                "PermDiffs" = MSE.diffs,
  #                "Importance" = importances,
  #                "Pvalue" = Pvals,
  #                "test_pts" = X.test,
  #                "weak_learner" = base.learner,
  #                "model_original" = rf_og,
  #                "model_permuted" = rf_pm,
  #                "test_stat" = "MSE",
  #                "call" = match.call())
  # class(result) <- "MSE_Test"
  # result
  return(list(Agg_Pvalue = p_out, Pvalues = out[,1]))
}

# rep_H0 <- sapply(1:100, function(x) stable_MSE_Test.default(X = data.frame(replicate(25, rnorm(500))), y = rnorm(500), var = c("X2"),
#                         n_val = 100, base.learner = "rtree", NTree = 500, p = .6, mtry = 5, B = 250, ranger = F))
# 
# 
# rep_H1 <-  sapply(1:100, function(x) {
#   X <- data.frame(replicate(25, rnorm(500)))
#   y <- .1*x*X[,1] + rnorm(500, sd = 0.5)
#   stable_MSE_Test.default(X, y = y, var = c("X1"),  n_val = 100, base.learner = "rtree", NTree = 250, p = .6, mtry = 5, B = 250, ranger = F)})
