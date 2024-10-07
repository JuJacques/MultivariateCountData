poissonmix <- function(X, G = 3, trace = FALSE)
{
  if (is.vector(X))
  {
    X <- matrix(X, length(X), 1)
  }
  N <- nrow(X)
  P <- ncol(X)

  # Cluster the data using k-means
  # Initial parameter "guesses" based on k-means results
  
  fit0 <- kmeans(X, centers = G, nstart = 50, iter.max = 50)
  tau <- table(fit0$cl) / N
  lambda <- fit0$centers

  #fit0 <- Mclust(X, G = G, prior = priorControl(), verbose = FALSE)
  #tau <- fit0$parameters$pro
  #lambda <- fit0$parameters$mean
  
  # Controls for EM Algorithm
  crit <- TRUE
  logL0 <- -Inf
  eps <- 10 ^ (-3)
  
  while (crit)
  {
    # E-step
    W <- matrix(NA, N, G)
    lW <- matrix(NA, N, G)
    lWa <- lW
    for (g in 1:G)
    {
      if (P > 1)
      {
        Lambda <- matrix(lambda[g, ], N, P, byrow = TRUE)
        lW[, g] <- log(tau[g]) + apply(dpois(X, Lambda, log = TRUE), 1, sum)
        #lWa[, g] <- log(tau[g]) + apply(apply(X, 1, dpois, lambda[g, ], log = TRUE), 2, sum)
      } else {
        #W[ ,g] <- tau[g] * dpois(X, lambda[g, ])
        lW[ ,g] <- log(tau[g]) + dpois(X, lambda[g, ], log = TRUE)
      }
    }
    lWmax <- apply(lW, 1, max)
    W <- exp(lW-lWmax)
    Wsum <- apply(W, 1, sum)
    Z <- W / Wsum
    
    # Compute log-likelihood
    #logL <- sum(log(Wsum))
    logL <- sum(lWmax) + sum(log(Wsum))
    if (trace) 
    {
      print(data.frame(G, logL))
    }
    # M-step
    Ng <- apply(Z, 2, sum)
    tau <- Ng / N
    lambda <- (t(Z) %*% X) / Ng
    
    # Check convergence
    crit <- (logL - logL0 > eps)
    logL0 <- logL	
  }
  
  BIC <- 2 * logL - log(N) * ((G - 1) + G * P)
  
  l <- map(Z)
  
  res <- list(Z = Z, classification = l, tau = tau, lambda = lambda, G = G, BIC = BIC, logL = logL)
  return(res)
}

poissonmix_all <- function(X, Gvec = 1:10, trace = FALSE)
{
  bestfit <- list()
  Gmax <- max(Gvec)
  BICvec <-rep(NA, length = Gmax)
  BICmax <- -Inf
  names(BICvec) <- 1:Gmax
  for (g in Gvec)
  {
    fit <- poissonmix(X, G = g, trace = trace)
    BICvec[g] <- fit$BIC
    if (fit$BIC > BICmax)
    {
      bestfit <- fit
      BICmax <- fit$BIC
    }
  }
  return(list(bestfit = bestfit, BIC = BICvec[Gvec]))
}

poissonmix_varsel <- function(X, jchosen = NULL, Gvec = 1:10)
{
  if (is.data.frame(X)) {X <- as.matrix(X)}
  if (is.null(jchosen)) {jchosen <- poissonmix_screen(X = X, Gvec = Gvec, threshold = 0)$jchosen}
  N <- nrow(X)
  P <- ncol(X)
  criterion <- TRUE
  iter <- 1
  
  list_jchosen=list(jchosen)
  
  print(paste("Initial Selected Variables: ", paste(jchosen, sep="," , collapse=","), sep=""))
  while (criterion)
  {
    print(paste("Iteration: ", iter, sep = ""))
    
    if (length(jchosen) < P)
    {
      # Forward
      fitpoissonNG <- poissonmix_all(X[ ,jchosen], Gvec = Gvec)
      BICdiff <- rep(-Inf, P)
      for (j in (1:P)[-jchosen])
      {
        fitpoissonG <- poissonmix_all(X[ ,c(jchosen, j)], Gvec = Gvec)
        
        y <- X[ ,j]
        Xregr <- X[ ,jchosen]
        datregr <- data.frame(y, Xregr)
        fitglm <- glm(y ~ ., family = poisson(link = "log"), data = datregr)
        fit1 <- glm(y ~ 1, family = poisson(link = "log"), data = datregr)
        modelscope <- list(lower = formula(fit1), upper = formula(fitglm))
        fitglmf <- step(fit1 , scope = modelscope, k = log(N), trace = 0, direction = "both")
        fitglmb <- step(fitglm, scope = modelscope, k = log(N), trace = 0, direction = "both")
        if (BIC(fitglmf) < BIC(fitglmb)) {fitglms <- fitglmf} else {fitglms <- fitglmb}
        
        predglms <- predict(fitglms, type = "response")
        BICglms <- 2 * sum(dpois(y, predglms, log = TRUE)) - log(N) * length(coef(fitglms))
        BICdiff[j] <- fitpoissonG$bestfit$BIC - (fitpoissonNG$bestfit$BIC + BICglms)
      }
      if (max(BICdiff) > 0)
      {
        j <- which.max(BICdiff)
        print(paste("Add Variable:", j, " BIC Difference:", round(max(BICdiff), 1)))
        jchosen <- c(jchosen, j)
        accept_forward <- TRUE
        
      } else {
        j <- which.max(BICdiff)
        print(paste("Add Variable: NONE", j, "BIC Difference:", round(max(BICdiff), 1)))
        accept_forward <- FALSE
      }
    } else {
      print("Add Variable: NONE  All variables in the model")
      accept_forward <- FALSE
    }
    
    # Backward
    if (length(jchosen) > 1)
    {
      if (length(jchosen) < P)
      {
        jremove <- setdiff(1:P,jchosen)
        fitpoissonG <- poissonmix_all(X[, -jremove], Gvec = Gvec)
        jrange <- (1:P)[-jremove]
      } else
      {
        jremove <- NULL
        fitpoissonG <- poissonmix_all(X, Gvec = Gvec)
        jrange <- (1:P)
      }
      BICdiff <- rep(-Inf, P)
      for (j in jrange)
      {
        fitpoissonNG <- poissonmix_all(X[ ,-c(jremove, j)], Gvec = Gvec)
        
        y <- X[ ,j]
        Xregr <- X[ ,-c(j, jremove)]
        datregr <- data.frame(y, Xregr)
        fitglm <- glm(y ~ ., family = poisson(link = "log"), data = datregr)
        fit1 <- glm(y ~ 1, family = poisson(link = "log"), data = datregr)
        modelscope <- list(lower = formula(fit1), upper = formula(fitglm))
        fitglmf <- step(fit1, scope = modelscope, k = log(N), trace = 0, direction = "both")
        fitglmb <- step(fitglm, scope = modelscope, k = log(N), trace = 0, direction = "both")
        if (BIC(fitglmf) < BIC(fitglmb)) {fitglms <- fitglmf} else {fitglms <- fitglmb}
        
        predglms <- predict(fitglms, type = "response")
        BICglms <- 2 * sum(dpois(y, predglms, log = TRUE)) - log(N) * length(coef(fitglms))
        BICdiff[j] <- (fitpoissonNG$bestfit$BIC + BICglms) - fitpoissonG$bestfit$BIC 
      }
      if (max(BICdiff) > 0)
      {
        j <- which.max(BICdiff)
        print(paste("Remove Variable:", j, " BIC Difference:", round(max(BICdiff), 1)))
        jchosen <- setdiff(jchosen, j)
        accept_backward <- TRUE
      } else {
        j <- which.max(BICdiff)
        print(paste("Remove Variable: NONE", j, "BIC Difference:", round(max(BICdiff), 1)))
        accept_backward <- FALSE
      }
    } else {
      print("Remove Variable: NONE. Only one variable in the model.")
      accept_backward <- FALSE
    }
    criterion <- (accept_forward) | (accept_backward)
    iter <- iter + 1
    print(paste("Current Selected Variables: ", paste(jchosen, sep="," , collapse=","), sep=""))
  
    # test for infinite loop
    if (criterion){
      isalready=FALSE
      for (i in 1:length(list_jchosen)) isalready=isalready+sum(identical(jchosen,list_jchosen[[i]]))
      if (isalready){
        criterion = FALSE;
        print("Stop due to infinite loop")
      } else {
        list_jchosen[[length(list_jchosen)+1]]=jchosen
      } 
    }
  }
  
  jchosen <- sort(jchosen)
  bestfit <- poissonmix_all(X[,jchosen], Gvec = Gvec)$bestfit
  return(list(jchosen = jchosen, bestfit = bestfit))
}

poissonmix_screen <- function(X, Gvec = 1:10, threshold = 0)
{
  P <- ncol(X)
  Gmax <- max(Gvec)
  BICmat <- matrix(NA, Gmax, P)
  for (j in 1:P)
  {
    Gtry <- min(length(table(X[, j])), Gmax)
    fit <- poissonmix_all(X[, j], Gvec = Gvec[Gvec <= Gtry])
    BICmat[Gvec[Gvec <= Gtry], j] <- fit$BIC
  }
  BICdiff <- apply(BICmat, 2, max, na.rm = TRUE) - BICmat[1, ]
  Gvec <- apply(BICmat, 2, which.max)
  jchosen <-  (1:P)[BICdiff > threshold]
  list(jchosen = jchosen, BICdiff = BICdiff, BICmat = BICmat)
}

