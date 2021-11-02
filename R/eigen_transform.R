#' Basis function dimension reduction
#' @title Dimension reduction using eigen-function
#' @author Zikai Lin
#' @param Y A numeric matrix
#' @param grids A matrix where rows represent points and columns are coordinates.
#' @param smoothness A positive real number specifies the smoothness parameter in the modified exponential squared kernel. The larger value the smoother Gaussian process is, the default is 0.5
#' @param poly_degree A integer number specifies the highest degree of Hermite polynomials. The default value is 5.
#' @param concentration A positive real number specifies the concentration parameter in the modified exponetial squared kernel. The larger value the more the GP concentrates around the center. The default value is 0.01.
#' @return A list object with
#' \itemize{
#'   \item Y_star a smaller matrix with mutually independent columns.
#'   \item svd single value decomposition of \code{grids}.
#' }
#' @examples Y <- matrix(rnorm(60),15,4)  # a outcome matrix with 3 samples (rows), each with 4 voxels (columns)
#' grids <- matrix(c(0.1,0.1,0.5,0.2, 0.3,0.3, 0.1,0.1),4,2) # 4 voxel locations
#' eigen.transform(Y, grids, poly_degree = 1)
#' 
#' @importFrom BayesGPfit GP.eigen.funcs.fast
#' 

eigen.transform <- function(Y, grids, smoothness = 0.5, poly_degree = 5, concentration = 0.01){
  b = 1/(2*smoothness^2) # applied to GP.eigen.funcs.fast()
  
  # Compute eigen functions for the standard modified exponential squared correlation kernel
  eigenfuns_fit = BayesGPfit::GP.eigen.funcs.fast(grids, poly_degree = poly_degree, a = concentration, b = b)
  
  # Compute SVD on the eigenfuns
  eigenfuns_fit.svd = svd(eigenfuns_fit)
  
  # transform outcome matrix Y
  Y.u  = Y  %*% eigenfuns_fit.svd$u
  
  return(list(Y_star = Y.u, svd = eigenfuns_fit.svd))
}