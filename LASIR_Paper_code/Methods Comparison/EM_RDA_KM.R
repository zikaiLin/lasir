#set.seed(2020)
library(flexmix)
library(gtools)
library(lasir)
#library(lasir)
rm(list = ls())
####### Input ##############
input_dat <-
  readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")



# Number of cluster (pre-defined)
G         <- input_dat$G
nG        <- input_dat$nG

# Number of observation
N         <- input_dat$N

J         <- input_dat$J
L         <- input_dat$L

# Covariates
X         <- input_dat$X[,-1]
Z         <- input_dat$Z[,-1]



V         <- nrow(input_dat$EigenFuns)
Y    <- input_dat$Y


#####################################################################
#      Get initial estimate of constant main effect                 #
#####################################################################

# 
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
#source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/updated_imgGPfit.R")
params_mbkm = list(batch_size = 10, init_fraction = 0.3, early_stop_iter = 10)
optKM = Optimal_Clusters_KMeans(input_dat$Y, max_clusters = 8, tol = 0.001, verbose =  T, mini_batch_params = params_mbkm, criterion = "distortion_fK")

kmeans_gp_fit = function(Y,
                         X, 
                         grids,
                         Z = NULL,
                         concomitant_type = c("none", "const", "mult"),
                         random_var=NULL,
                         ncomp = 3,
                         d = 3,
                         poly_degree = 5,
                         a = 0.0001,
                         smoothness = 0.05,
                         step.nrep = 10,
                         verbose = TRUE,
                         intercept = FALSE,
                         eigens_trans = NULL,
                         ...){
  
  concomitant_type = match.arg(concomitant_type)
  
  # dimension of input data
  if(intercept){
    N = nrow(X)
    V = ncol(Y)
    J = ncol(X) - 1
    Q = ncol(Z) - 1
    X_fit = matrix(X[,-1], N, J)
  }else{
    N = nrow(X)
    V = ncol(Y)
    J = ncol(X)
    Q = ncol(Z)
    X_fit = matrix(X, N, J)
  }
  
  Z_fit = matrix(Z, N, Q)
  
  
  
  
  # Dimension reduction 
  if(is.null(eigens_trans)){
    Y.eigen.trans = eigen.transform(
      Y = Y,
      grids = grids,
      smoothness = smoothness,
      poly_degree = poly_degree)
  }else{
    Y.eigen.trans = eigens_trans
  }
  
  Y_star = Y %*% Y.eigen.trans$svd$u
  L = ncol(Y_star)
  
  
  # Step 1: regress out demogrphic variable to obtain initial estimate of fixed main coef
  
  fixed_main_fit = .voxelwise_lm_fit_all(Y_star,  Z_fit) # The intercept should be included
  fixed_main_coef = fixed_main_fit$coef
  fixed_main_residuals = fixed_main_fit$residual
  
  
  # Model selection (# of components)
  km.fit = ClusterR::KMeans_rcpp(Y, ncomp, verbose = T)
  
  
  # The optimal number of components
  nK = ncomp
  
  # Best model with least BIC
  grp_ind   <- km.fit$cluster
  
  # Extract theta_alpha coefficient and covariance matrix from kmean model
  y_star_by_grp  <- list()
  theta_alpha_post = array(dim=c(J+1, L, nK))
  
  for (k in 1:nK) {
    y_star_by_grp[[k]] = fixed_main_residuals[km.fit$cluster == k, ]
    theta_alpha_post[,,k] = voxelwise_lm_fit(y_star_by_grp[[k]],  X_fit[km.fit$cluster == k, ])
  }
  
  
  
  if(!is.null(random_var)){
    nG = max(unique(random_var))
    
    Y_fit_gamma = fixed_main_residuals
    theta_alpha_post_gamma = theta_alpha_post
    cov_post = array(NA, dim = c(ncol(Y_star), ncol(Y_star), nK))
    
    Y_fit_gamma.random = Y_fit_gamma - t(sapply(1:N,  function(i)
      c(1, X_fit[i, ]) %*% theta_alpha_post_gamma[, , grp_ind[i]]))
    
    theta_gamma_post = t(sapply(1:nG, function(g) {
      if(sum(random_var == g) == 0){
        return(rep(0, ncol(Y_fit_gamma)))
      }else if(sum(random_var == g) == 1){
        return(Y_fit_gamma.random[random_var == g, ])
      }else{
        apply(Y_fit_gamma.random[random_var == g, ], 2, mean)
      }
    })) 
    
  }else{
    gamma_coef_post = NULL
  }
  
  
  # Iteratively update random intercept, fixed effect and main effect
  
  Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
  tol = 1e-4
  eps = Inf
  pre_fixed_main_coef = fixed_main_coef
  pre_fixed_coef      = theta_alpha_post
  pre_random_coef     = theta_gamma_post
  
  k = 0
  while ((eps >= tol) & (k<10)) {
    k = k+1
    cat("Iteration: ", k, "....\n")
    
    cat("Updating Alpha.\n")
    # 1. update alpha
    post_fixed_coef = update_alpha_coef(
      gamma_coef_post_all = pre_random_coef[random_var, ],
      fixed_coef_post = pre_fixed_main_coef,
      X = cbind(1, X_fit),
      Z =  Z_fit,
      Y = Y_star,
      grp_ind = grp_ind,
      nK = nK,
      sscp = T
    )
    
    
    cat("Updating main coef. \n")
    post_fixed_main_coef = update_fixed_coef(
      alpha_coef_post = post_fixed_coef$alpha_coef,
      gamma_post_all = pre_random_coef[random_var, ],
      X = cbind(1, X_fit),
      Z = Z_fit,
      Y = Y_star,
      grp_ind = grp_ind,
      nK = nK
    )
    
    cat("Updating random coef. \n")
    post_random_coef = update_gamma_coef(
      alpha_coef_post = post_fixed_coef$alpha_coef, 
      fixed_coef = post_fixed_main_coef$coef,
      X = cbind(1, X_fit),
      Z =  Z_fit,
      Y = Y_star,
      grp_ind = grp_ind,
      sites = random_var,
      nK = nK,
      nG = nG
    )
    
    eps = mean((post_fixed_coef$alpha_coef - pre_fixed_coef)^2)
    pre_fixed_coef = post_fixed_coef$alpha_coef
    pre_random_coef = post_random_coef
    pre_fixed_main_coef = post_fixed_main_coef$coef
    covariance_post = post_fixed_coef$SSCPE
    cat("eps = ", eps, "\n")
  }
  
  
  # Reconstruct fixed effect coefficient
  theta_alpha_post = post_fixed_coef$alpha_coef
  gamma_coef_post = post_random_coef %*% t(Y.eigen.trans$svd$u)
  main_coef_post = post_fixed_main_coef$coef %*% t(Y.eigen.trans$svd$u)
  
  alpha_coef_post = array(NA, dim = c(J + 1, V, nK))
  
  for (k in 1:nK) {
    alpha_coef_post[, , k] = (theta_alpha_post[, , k] %*% t(Y.eigen.trans$svd$u))
  }
  
  
  # # Obtain uncorrected p-value
  
  # alpha_coef_post_var = fixed_p_corr(
  #   theta_alpha_post,
  #   covariance_post,
  #   X = cbind(1, X_fit),
  #   grp_ind = grp_ind,
  #   svd_u = Y.eigen.trans$svd$u
  # )
  # 
  # alpha_coef_post_pval = alpha_coef_post_tval = array(dim=dim(alpha_coef_post_var))
  # alpha_coef_post_sd = sqrt(alpha_coef_post_var)
  # 
  # for(j in 1:(J+1)){
  #   alpha_coef_post_tval[j,,] = abs(alpha_coef_post[j,,])/alpha_coef_post_sd[j,,]
  #   alpha_coef_post_pval[j,,] = (1-stats::pnorm(abs(alpha_coef_post[j,,])/alpha_coef_post_sd[j,,],0,1))*2
  # }

  concomitant_mod = nnet::multinom(factor(grp_ind) ~ Z_fit)
  
  res = list(
    model_best = km.fit,
    clusters = grp_ind,
    main_coef = alpha_coef_post,
    # fixed_pval = alpha_coef_post_pval,
    # fixed_tval = alpha_coef_post_tval,
    # fixed_sd   = alpha_coef_post_sd,
    fixed_coef = main_coef_post,
    site_coef = gamma_coef_post,
    concomitant_mod = concomitant_mod
  )
  
  return(res)
}  

#####################################################################
#                         EM Classification                         #
#####################################################################
eigens_trans = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/transform_res_p14_rho05.rds")

m1 = kmeans_gp_fit(
  Y = input_dat$Y,
  X = input_dat$X_g,
  grids = input_dat$loc,
  d = 3,
  poly_degree = 14,
  random_var = input_dat$G,
  ncomp = 3,
  smoothness = 0.05,
  step.nrep = 10,
  Z = input_dat$X_demogrphic,
  concomitant_type = "mult",
  criteria = "BIC",
  verbose = T,
  intercept = F,
  eigens_trans = eigens_trans
)




list_of_vars_tosave = c("m1")

save(list = list_of_vars_tosave,
     file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/result_p14_rho05_KM.RData")


