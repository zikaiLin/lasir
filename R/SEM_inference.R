#' Model Inference of LASIR
#' @title Model Inference of Main Effect Coefficients obtained from LASIR
#' @author Zikai Lin
#' @description Obtain standard error estimates and spatially varying coefficient (SVCs) from basis coefficients estimates
#' @param lasir_obj Object obtained from
#' @param eigen_transformation Eigen transformation from original image space to orthonormal space, can be obtained from eigen.transformation function in the package
#' @return A list object
#' \describe{
#'   \item{main_effect_svcs}{Main effect spatially varying coefficients}
#'   \item{sd}{Estimated standard error at each location of main effect coefficients}
#'   \item{pvals}{Uncorrected p-values of main effect at each image location.}
#' }
#'
#' @examples
#' Y <- matrix(rnorm(75),15,5)  # a outcome matrix with 3 samples (rows), each with dimension of 4
#' X <- matrix(rnorm(15),15,1) # Covariate matrix
#' Z <- matrix(rnorm(15),15,1) # Fixed effect covariate matrix
#“ Y_fit <- matrix(rnorm(60),15,4)  # a outcome matrix with 3 samples (rows), each with 4 voxels (columns)
#' grids <- matrix(c(0.1,0.1,0.5,0.2, 0.3,0.3, 0.1,0.1),4,2) # 4 voxel locations
#' Sites <- sample(1:3, size = 15, replace = T) # Site indicators, assumes 15 subjects from 3 study sites
#' nK <- 2 # Assume 2 subgroups
#'
#' eigen_trans <- eigen.transform(Y, grids, poly_degree = 1)
#' lasir_fit <- LASIR.fit(Y,X,Z,Sites,nK,max.iter=20,loglik_tol=1e-2)
#' lasir.infer <- LASIR.inference(lasir_fit, eigen_trans)
#'
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @export
#'

LASIR.inference = function(lasir_obj, eigen_transformation){

  if(class(lasir_obj) != "lasir") stop("Input is not a lasir object.")

  main_basis_coef = lasir_obj$main_basis_coef

  dim_main_basis = dim(main_basis_coef)
  J = dim_main_basis[1]
  L = dim_main_basis[2]
  nK = dim_main_basis[3]
  V = nrow(eigen_transformation$svd$u)

  X = lasir_obj$dat$X

  n_pars_main_basis = J*L


  clusters = lasir_obj$cluster

  # Observed information matrix on basis coefficient
  observed_info_mat = inv_observed_info_mat = array(0, dim = c(n_pars_main_basis, n_pars_main_basis, nK))
  for(k in 1:nK){
    kronecker_cross= kronecker(crossprod(cbind(1, X[clusters == k,])), lasir_obj$sscpe)
    diag(kronecker_cross) = diag(kronecker_cross+0.001)
    inv_observed_info_mat[,,k] = solve(kronecker_cross)
  }


  # Variance estimation
  var_est_main_effect = array(NA, dim = c(J, V, nK))

  for(k in 1:nK){
    for(j in 1:J){
      subblock_idx = ((j-1)*L+1):(j*L)
      inv_obs_info_mat_subblock = inv_observed_info_mat[subblock_idx, subblock_idx,k]
      var_est_main_effect[j, ,k] = rowSums(eigen_transformation$svd$u * t(inv_obs_info_mat_subblock %*% t(eigen_transformation$svd$u)) )
    }
  }

  sd_est_main_effect = sqrt(var_est_main_effect)


  # Main effect SVCs estimation
  main_effect_coef =  array(NA, dim = c(J, V, nK))

  for(k in 1:nK){
    main_effect_coef[,,k] = lasir_obj$main_basis_coef[,,k] %*% t(eigen_transformation$svd$u)
  }

  # Obtain p-value estimation
  pval_est_main_effect = (1-pnorm(abs((main_effect_coef/sd_est_main_effect)),0,1) )*2 # 2-tail z-test

  res = list(main_effect_svcs = main_effect_coef,
             sd = sd_est_main_effect,
             pvals = pval_est_main_effect)
  return(res)
}



