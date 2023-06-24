# Cross validation wtih residual test
library(lasir)
library(flexmix)
library(BayesGPfit)




source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mvEMfix_r copy.R")

# Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
# Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/scripts/fit_predict.cpp")
source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/section2-EM-run/kmeans_gp_fit.R")
input_dat = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")

load("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/folds.RData")

BIC_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/BIC_list_res.rds")
model_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/model_list_res.rds")
best_mod_idx = which(BIC_list_res == min(BIC_list_res), arr.ind = T)
mod = model_list_res[[best_mod_idx[1]]][[best_mod_idx[2]]]



# 
n_train = floor(input_dat$N*0.9)
V = ncol(input_dat$Y)
n_fold = 50

# --------------------- step 1: generate partitions -------------------------- #



mse_y_lasir_shuffle = rep(NA, n_fold)

eigens_trans =  readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/transform_res_p14_rho05.rds")

# The latent class is known
clusters = mod$cluster
clusters_shuffle = clusters[permute::shuffle(length(clusters))]
nK = length(unique(clusters_shuffle))

for(f in 1:50){
  try({
    cat("Iteration: ", f, "\n")
    
    ##### Training ##########
    fold_ind = folds_50_train[[f]]
    Y_train = input_dat$Y[fold_ind,]
    eigens_trans_f = eigens_trans
    eigens_trans_f$Y_star = eigens_trans_f$Y_star[fold_ind,]
    
    X_gp_train = matrix(input_dat$X_g[fold_ind,], ncol = 1)
    X_demo_train = input_dat$X_demogrphic[fold_ind,]
    G_train  = input_dat$G[fold_ind]
    nG = length(unique(G_train))
    
    clusters_f = clusters_shuffle[fold_ind]
    
    # ---------- LASIR ---------- #
    
    N = nrow(X_gp_train)
    V = ncol(Y_train)
    J = ncol(X_gp_train)
    Q = ncol(X_demo_train)
    X_fit = matrix(X_gp_train, N, J)
    Z_fit = matrix(X_demo_train, N, Q)
    
    
    
    main_coef = array(NA, dim = c(J+1, V, nK))
    fixed_coef = matrix(NA, Q , V)
    site_coef = matrix(NA, nG, V)
    
    Y_star = eigens_trans_f$Y_star
    L = ncol(Y_star)
    
    
    EM <- mvnMixSEMFix_fix$new(outcome= Y_star,
                               exposed_var = X_fit,
                               control_var = Z_fit,
                               sites = G_train,
                               num_components = nK,
                               preSpec_cluster = clusters_f)
    
    res_EM <- EM$run.EM(loglik_tol=1e-4)
    
    for(k in 1:nK){
      main_coef[,,k] = res_EM$main_basis_coef[,,k] %*% t(eigens_trans$svd$u) 
    }
    
    site_coef = res_EM$sites_basis_coef %*% t(eigens_trans$svd$u) 
    fixed_coef= res_EM$fixed_basis_coef %*% t(eigens_trans$svd$u) 
    
    

    
    # ------------------------- Testing -------------------------- #
    fold_ind_test = folds_50_test[[f]]
    Y_test = input_dat$Y[fold_ind_test,]
    X_gp_test = matrix(input_dat$X_g[fold_ind_test,], ncol = 1)
    X_demo_test = input_dat$X_demogrphic[fold_ind_test,]
    G_test  = input_dat$G[fold_ind_test]
    
    
    # -------- Testing of LASIR --------- #
    Y_pred = predict_gp_known_latent(x = cbind(1, X_gp_test), 
                                     z = X_demo_test,
                                     G = G_test,
                                     fixed_coef  = fixed_coef,
                                     main_coef  = main_coef,
                                     gamma_coef = site_coef,
                                     grp_ind = clusters_shuffle[fold_ind_test]
    )
    
    mse_y_lasir_shuffle[f] = mean((Y_test - Y_pred$y_pred)^2)
    
    
    cat("MSE (lasir shuffle): ", mse_y_lasir_shuffle[f], "\n")
    
    
  })
  
  
  
}

save(list = "mse_y_lasir_shuffle",
     file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_lasir_cross_validation_shuffle.RData")

# save(list = "mse_y_svcm",
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_svcm_cross_validation.RData")
# 
# 
# save(list = "mse_y_kmlr",
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_kmlr_cross_validation.RData")

# save(list = c("mse_y_lasir","mse_y_svcm", "mse_y_kmlr", "folds_50_train", "folds_50_test"),
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_cross_validation.RData")

# cv_mse_bind = cbind(mse_y_lasir, mse_y_lasir_shuffle, mse_y_svcm)
# colnames(cv_mse_bind) = c("LASIR", "KMLR", "SVCM")
# 
# boxplot(cv_mse_bind, ylab ="Square of residuals",  cex.lab = 1.5, cex.axis = 1.5)
