
kmeans_gp_fit = function(Y,
                         X, 
                         grids,
                         Z = NULL,
                         random_var=NULL,
                         ncomp = 3,
                         d = 3,
                         poly_degree = 5,
                         a = 0.01,
                         smoothness = 0.05,
                         step.nrep = 10,
                         verbose = TRUE,
                         intercept = FALSE,
                         eigens_trans = NULL,
                         ...){
  
  #concomitant_type = match.arg(concomitant_type)
  
  # dimension of input data
  if(intercept){
    N = nrow(X)
    V = ncol(Y)
    J = ncol(X) - 1
    Q = ncol(Z) - 1
    X_fit = matrix(X[,-1], N, J)
    Z_fit = matrix(Z[,-1], N, Q)
  }else{
    N = nrow(X)
    V = ncol(Y)
    J = ncol(X)
    Q = ncol(Z)
    X_fit = matrix(X, N, J)
    Z_fit = matrix(Z, N, Q)
  }
  
  
  
  
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
  
  fixed_main_fit = voxelwise_lm_fit_all(Y_star, cbind(1, Z_fit)) # The intercept should be included
  fixed_main_coef = fixed_main_fit$coef
  fixed_main_residuals = fixed_main_fit$residual
  
  
  # Model selection (# of components)
  opt_km = ClusterR::Optimal_Clusters_KMeans(fixed_main_residuals, max_clusters = max(ncomp), criterion = "distortion_fK")
  km.fit = ClusterR::KMeans_rcpp(fixed_main_residuals, which.min(opt_km))
  
  
  # The optimal number of components
  nK = length(unique(km.fit$clusters))
  
  # Best model with least BIC
  grp_ind   <- km.fit$cluster
  
  # Extract theta_alpha coefficient and covariance matrix from kmean model
  y_star_by_grp  <- list()
  theta_alpha_post = array(dim=c(J+1, L, nK))
  
  for (k in 1:nK) {
    if(sum(km.fit$cluster == k) == 1){
      # Only one observation observed
      warning("Cluster size equal to 1.")
      y_star_by_grp[[k]] = fixed_main_residuals[km.fit$cluster == k, ] # Should be a 1xL element
      theta_alpha_post[1,,k] = y_star_by_grp[[k]]
      theta_alpha_post[2,,k] = rep(0, L)
      next()
    }
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
      main_coef_post = pre_fixed_main_coef,
      X = cbind(1, X_fit),
      Z =  cbind(1,Z_fit),
      Y = Y_star,
      grp_ind = grp_ind,
      nK = nK,
      sscp = T
    )
    
    
    cat("Updating main coef. \n")
    post_fixed_main_coef = update_main_coef(
      alpha_coef_post = post_fixed_coef$alpha_coef,
      gamma_post_all = pre_random_coef[random_var, ],
      X = cbind(1, X_fit),
      Z = cbind(1,Z_fit),
      Y = Y_star,
      grp_ind = grp_ind,
      nK = nK
    )
    
    cat("Updating random coef. \n")
    post_random_coef = update_gamma_coef(
      alpha_coef_post = post_fixed_coef$alpha_coef,
      main_coef = post_fixed_main_coef$coef,
      X = cbind(1, X_fit),
      Z =  cbind(1,Z_fit),
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
    fixed_coef = main_coef_post,
    # fixed_pval = alpha_coef_post_pval,
    # fixed_tval = alpha_coef_post_tval,
    # fixed_sd   = alpha_coef_post_sd,
    main_coef = alpha_coef_post,
    random_coef = gamma_coef_post,
    concomitant_mod = concomitant_mod
  )
  
  return(res)
}  
