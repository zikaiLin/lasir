#' Model fit with SEM
#' @title Latent Subgroup Image on Scalar Regression using SEM
#' @author Zikai Lin
#' @description Use stochastic expectation maximization (SEM) algorithm to obtain maximum likelihood estimate of image-on-scalar regression
#' @param Y Multivariate outcome matrix
#' @param X Covariate matrix
#' @param Z Control variable matrix
#' @param sites Integer vector, sites indicator, treated as study site fixed effect
#' @param nK Integer, number of subgroups
#' @param max.iter Integer, Maximum number of iteration in stochastic 
#' @param loglik_tol Numeric, tolerance of log-likelihood change, SEM will stop if the change of log-likelihood falls below this value
#' @return A list object
#' \describe{
#'  \item{dat}{Input data}
#'   \item{convergence}{Whether SEM converge (0) or failed to converged (1)}
#'   \item{cluster}{Estimated subgroup indicator}
#'   \item{main_basis_coef}{Estimated subgroup specific coefficient of  main effect}
#'   \item{fixed_basis_coef}{Estimated coefficient of fixed effect}
#'   \item{sites_basis_coef}{Estimated site-specific coefficient}
#'   \item{logit_wts}{Estimated coefficients of second-level logit model}
#'   \item{sigma_vec}{Variance parameters of each dimension}
#'   \item{sscpe}{Estimated error covariance}
#'   \item{prob_mat.y}{a-Posterior probability}
#'   \item{loglik_list}{Track of log-likelihood at each SEM iteration.}
#'   \item{BIC}{BIC criterion of the model}
#' }
#' 
#' @examples 
#' Y <- matrix(rnorm(75),15,5)  # a outcome matrix with 3 samples (rows), each with dimension of 4
#' X <- matrix(rnorm(15),15,1) # Covariate matrix
#' Z <- matrix(rnorm(15),15,1) # Fixed effect covariate matrix
#' Sites <- sample(1:3, size = 15, replace = T) # Site indicators, assumes 15 subjects from 3 study sites
#' nK <- 2 # Assume 2 subgroups
#' 
#' LASIR.fit(Y,X,Z,Sites,nK,max.iter=20,loglik_tol=1e-2)
#' @importFrom stats dnorm
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom stats sd
#' @importFrom nnet multinom
#' @importFrom nnet predict.multinom
#' @importFrom nnet coef.multinom
#' @importFrom nnet logLik.multinom
#' @importFrom stats residuals
#' @export
#' 
LASIR.fit =  function(Y, 
                   X, 
                   Z,
                   sites,
                   nK = 3,
                   max.iter = 20, 
                   loglik_tol = 0.01) {

  ##################################################################
  #                              Initialize                        #
  ##################################################################
  
  x = X
  z = Z
  g = sites
  y = Y
  N = nrow(X)
  L = ncol(Y)
  J = ncol(X) + 1
  Q = ncol(Z) 
  nG = length(unique(sites))
  
  prob_mat.logit <- matrix(1.0/nK, N, nK) # weight matrix (from logit model)
  prob_mat.y <- matrix(1.0/nK, N, nK) # weight matrix
  sigma_vec <- apply(y, 2, sd, na.rm = T) # for each dimension, there's a variance term
  
  comp_indicator <- sample(c(1:nK), size = N, replace = TRUE)
  
  # Residuals for regressing out fix effect
  fit.fix.resids <- matrix(0, N, L)
  fit.resids <- matrix(0, N, L)
  
  # basis function coefficients
  main_basis <- array(0, dim = c(J, L, nK))
  fixed_basis <- matrix(0, nrow = Q, ncol = L)
  sites_basis <- matrix(0, nrow = nG, ncol = L)
  
  # multinomial weights
  logit_wts <- matrix(0, nrow = nK-1, ncol = L)
  
  # Log likelihoods
  loglik.logit <- -Inf
  loglik.y <- -Inf
  loglik <- -Inf
  loglik_k.mat <- matrix(-Inf, N, nK)
  sigma_k <- matrix(NA, nK, L)
  loglik_list = NULL

  
  # Check tolerance function, if the change loglikelihood is below a specific tolerance
  
  check.tol = function(fmax,fmin,ftol){
    delta = abs(fmax - fmin)
    accuracy = (abs(fmax) + abs(fmin))*ftol
    return(delta < (accuracy + .Machine$double.eps))
  }
  
  
  # Check tolerance function, if the change loglikelihood is below a specific tolerance
  getBIC = function(loglik, N, nK, J, nG, Q, L){
    numParams = L*(nK*(J+1) + nG + Q) + (nK-1)*(Q+1) + L
    
    numParams*log(N*L) - 2*loglik
  }

   
  if(nK == 1)  comp_indicator <- rep(1,N)
  
  convergence = 1
  for(iter in 1:max.iter){
    loglik0 <- loglik
    
    #---------- fit.fix() ----------#
    y.x = t(sapply(1:N, function(i){
      y[i,] - cbind(1,x[i,]) %*% main_basis[,,comp_indicator[i]]
    }))
    
    fit.fix.lm = lm(y.x ~ factor(g) + z)
    
    sites_basis <- coef(fit.fix.lm)[1:nG,]
    
    # Put restriction on sites basis
    sites_basis <- sites_basis - rowMeans(sites_basis)
    
    # Update fixed_basis
    fit.fix.lm.refit = lm(y.x - sites_basis[g, ] ~ z - 1)
    fixed_basis <- coef(fit.fix.lm.refit)
    fit.fix.resids <- y - sites_basis[g, ] - z%*%fixed_basis
    
    
    #---------- fit.main() ----------#
    fit.lm = list()
    
    for(k in 1:nK){
      if(sum(comp_indicator == k) == 0){
        main_basis[,,k] <- 0
        next()
      } 
      
      fit.lm[[k]] = lm(fit.fix.resids ~ x, subset = (comp_indicator == k))
      main_basis[,,k] <- coef(fit.lm[[k]])
      fit.resids[comp_indicator == k,] <- residuals(fit.lm[[k]])
      
    }
    
    
    #---------- fit.logit() ----------#
    
    # Second level multinomial logit model
    if(nK > 1){
      multinom.fit = nnet::multinom(factor(comp_indicator) ~ z, trace = FALSE)
      logit_wts <- coef(multinom.fit)
      loglik.logit <- as.numeric(logLik(multinom.fit))
      prob_mat.logit <- matrix(predict(multinom.fit, type = "probs"), N, nK)
    }

    
    
    
    #---------- update.sigma() ----------#
    
    sigma_vec <- sqrt(colSums(fit.resids^2) / (N - (J+1+Q+nG)))
    sigma_mat <- crossprod(fit.resids)/(N-(J+1+Q+nG))
    
    
    #---------- update.prob() ----------#
    if(nK > 1){
      loglik_k.mat <- sapply(1:nK, function(k){
        
        res_k = fit.fix.resids - cbind(1,x) %*% main_basis[,,k] 
        sigma_k[k,] <- sqrt(colSums(res_k[comp_indicator == k,]^2)/(sum(comp_indicator == k) - (J+1)))
        rowSums(dnorm(res_k, mean = 0, matrix(rep(sigma_k[k,], N), N, ncol = L, byrow = T), log = T), na.rm = T)
        # rowSums(dnorm(res_k, mean = 0, matrix(rep(sigma_k[k,], N), N, ncol = L, byrow = T), log = T))
      })
      
      
      prob_mat.y <- log(prob_mat.logit) + loglik_k.mat
      
      max_log_prob <- apply(prob_mat.y, 1, max)
      
      prob_mat.y <- exp(prob_mat.y - max_log_prob)
      sum_prob <- apply(prob_mat.y, 1, sum)
      prob_mat.y <- prob_mat.y / sum_prob
      
    }else{
      prob_mat.y <- exp(prob_mat.y)
    }


    
    #---------- update.comp.ind() ----------#
    prob_mat.y[is.na(prob_mat.y)] <- 0
    if(nK > 1)  comp_indicator <- apply(prob_mat.y, 1, function(x) sample(1:nK, size = 1, prob = x))
    
    #---------- update.lik.y()  ----------#
    if(nK == 1){
      loglik.y.mat <- loglik_k.mat
    }else{
      loglik.y.mat <- log(prob_mat.y) + loglik_k.mat
    }
    
    
    if(nK > 1){
      maxLoglik.y <- apply(loglik.y.mat, 1, max)
      explik.y.mat <- exp(loglik.y.mat - maxLoglik.y)
      loglik.y <- log(apply(explik.y.mat, 1, sum)) + maxLoglik.y
      
    }else{
      loglik.y <- sapply(1:1, function(k){
        res_k = fit.fix.resids - cbind(1,x) %*% main_basis[,,k] #- sites_basis[g,] - z %*% fixed_basis
        sigma_k = sqrt(colSums(res_k[comp_indicator == k,]^2)/(sum(comp_indicator == k) - (J+1)))
        rowSums(dnorm(res_k, mean = 0, matrix(rep(sigma_k, N), N, ncol = L, byrow = T), log = T))
      })
    }
    
    
    
    #---------- update.loglik()  ----------#
    loglik.logit.obs <- log(prob_mat.y) +  log(prob_mat.logit)
    max.loglik.logit.obs = apply(loglik.logit.obs, 1, max)
    explik.logit.obs <- exp(loglik.logit.obs - max.loglik.logit.obs)
    loglik.logit.obs <- log(apply(explik.logit.obs, 1, sum)) + max.loglik.logit.obs
    
    loglik <- ifelse(nK>1, sum(loglik.y + loglik.logit.obs), sum(loglik.y))
    
    
    # cat("iter: ", iter, "loglik = ", loglik, "\n")
    loglik_list = c(loglik_list,loglik)
    if(check.tol(loglik0,loglik,loglik_tol)){
      convergence = 0
      break
    }
    
    
    
  }
  
  dat = list(Y = Y,
             X = X,
             Z = Z,
             sites = sites)
  
  if(nK > 1){
    res = list(dat = dat,
               convergence=convergence,
               cluster = comp_indicator,
               main_basis_coef = main_basis,
               fixed_basis_coef = fixed_basis,
               sites_basis_coef = sites_basis,
               logit_wts = logit_wts,
               sigma_vec=sigma_vec,
               sscpe = sigma_mat,
               prob_mat.y = prob_mat.y,
               iter=iter,
               loglik_list=loglik_list,
               BIC = getBIC(loglik,N,nK,J,nG,Q,L))
  }else{
    res = list(dat = dat,
               convergence=convergence,
               cluster = comp_indicator,
               main_basis_coef = main_basis,
               fixed_basis_coef = fixed_basis,
               sites_basis_coef = sites_basis,
               sigma_vec=sigma_vec,
               sscpe = sigma_mat,
               iter=iter,
               loglik_list=loglik_list,
               BIC = getBIC(loglik,N,nK,J,nG,Q,L))
    
    
  }
  class(res) = "lasir"
  return(res)
}





#' Model selection on number of cluster with BIC
#' @title Select Optimal Number of Subgroup for Latent Subgroup Image on Scalar Regression using SEM
#' @author Zikai Lin
#' @description Use BIC to select the optimal number of subgroups in LAtent Subgroup Image-on-scalar Regression (LASIR).
#' @param Y Multivariate outcome matrix
#' @param X Covariate matrix
#' @param Z Control variable matrix
#' @param sites Integer vector, sites indicator, treated as study site fixed effect
#' @param nK_list Integer vector, candidates of numbers of subgroups \eqn{K}
#' @param nrep Integer, number of repetition runs for each candidates \eqn{k\in \{1,\ldots,K\}}.
#' @param max.iter Integer, Maximum number of iteration in stochastic 
#' @param loglik_tol Numeric, tolerance of log-likelihood change, SEM will stop if the change of log-likelihood falls below this value
#' @return A list object
#' \describe{
#'   \item{best_model}{Model with minimum BIC}
#'   \item{BIC_list}{BIC criterion of the model}
#' }
#' 
#' @examples 
#' # a outcome matrix with 3 samples (rows), each with dimension of 2
#' Y <- rbind(matrix(rnorm(50),10,5),
#'            matrix(rnorm(25, mean = 2, sd = 3),5,5))  
#' X <- matrix(rnorm(15),15,1) # Covariate matrix
#' Z <- matrix(rnorm(15),15,1) # Fixed effect covariate matrix
#' Sites <- sample(1:3, size = 15, replace = T) # Site indicators, assumes 15 subjects from 3 study sites
#' 
#' selectK(Y,X,Z,Sites,1:2,max.iter=20,loglik_tol=1e-2)
#' @export
#' 

selectK = function(Y, 
                   X, 
                   Z,
                   sites,
                   nK_list = 1:3,
                   nrep = 5,
                   max.iter = 20, 
                   loglik_tol = 0.01){
  
  
  if(class(nK_list) != "integer") stop("Non-integer value(s) detected in nK_list, please check.")
  
  # A matrix containing list of BIC for each candidate model
  BIC_list = matrix(Inf, nrow = length(nK_list), ncol = nrep)
  minBIC = Inf
  bestK = NA
  # Best Model candidate
  Best_Model = NULL
  for(k in 1:length(nK_list)){
    for(r in 1:nrep){
      fit_kr = SEM.ios.fit(Y, X, Z, sites, nK_list[k], max.iter, loglik_tol)
      BIC_list[k,r] = fit_kr$BIC
      if(fit_kr$BIC < minBIC){
        Best_Model = fit_kr
        bestK = nK_list[k]
        minBIC = fit_kr$BIC
      } 
    }
  }
  
  
  return(list(BIC_list = BIC_list,
              best_model = Best_Model))
  
}


