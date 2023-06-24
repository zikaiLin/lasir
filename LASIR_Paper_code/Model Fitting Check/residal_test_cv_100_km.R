# Cross validation wtih residual test
library(lasir)
library(flexmix)
library(BayesGPfit)

Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/scripts/fit_predict.cpp")
source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/section2-EM-run/kmeans_gp_fit.R")
input_dat = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")
load("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/folds.RData")

# load("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/section4-cross-validation/mse_cross_validation.RData")

# 
# n_train = floor(input_dat$N*0.9)
# V = ncol(input_dat$Y)
# n_fold = 50
# 
# # step 1: generate partitions
# folds_50_train = list()
# folds_50_test = list()
# 
# for(i in 1:n_fold){
#   folds_50_train[[i]] = sample(input_dat$N, size = n_train)
#   folds_50_test[[i]] = setdiff(1:input_dat$N, folds_50_train[[i]])
# }
#  
# mse_y_lasir = mse_y_svcm = mse_y_kmlr = rep(NA, n_fold)

eigens_trans = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/transform_res_p14_rho05.rds")

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
    

    # # ---------- KMLR ---------- #
    mod_kmlr =  kmeans_gp_fit(
      Y = Y_train,
      X = X_gp_train,
      grids = input_dat$loc,
      d = 3,
      poly_degree = 12,
      random_var = G_train,
      ncomp = 1:6,
      smoothness = 0.05,
      step.nrep = 5,
      Z = X_demo_train,
      concomitant_type = "mult",
      criteria = "ICL",
      eigens_trans = eigens_trans_f,
      verbose = T,
      intercept = F
    )
    
    
    
    
    # ------------------------- Testing -------------------------- #
    # fold_ind_test = folds_50_test[[f]]
    # Y_test = input_dat$Y[fold_ind_test,]
    # X_gp_test = matrix(input_dat$X_g[fold_ind_test,], ncol = 1)
    # X_demo_test = input_dat$X_demogrphic[fold_ind_test,]
    # G_test  = input_dat$G[fold_ind_test]
    
    
    # -------- Testing of KMLR --------- #
    # nK    = length(unique(mod_kmlr$clusters))
    # N     = length(fold_ind_test)
    # # Probability allocation
    # prob_mat = matrix(NA, N, nK)
    # 
    # for(i in 1:N){
    #   prob_mat[i,] = c(1, X_demo_test[i,]) %*% cbind(0,t(coef(mod_kmlr$concomitant_mod)))
    #   prob_mat[i,] = exp(prob_mat[i,] - max(prob_mat[i,]))
    #   prob_mat[i,] = prob_mat[i,]/sum(prob_mat[i,], na.rm = T)
    # }
    
    Y_pred_kmlr = predict_gp_known_latent(x = cbind(1, X_gp_train), 
                                          z = cbind(1, X_demo_train),
                                          G = G_train,
                                          fixed_coef  = mod_kmlr$fixed_coef,
                                          main_coef  = mod_kmlr$main_coef,
                                          gamma_coef = mod_kmlr$random_coef,
                                          # prob_mat = prob_mat,
                                          # y_train = Y_train,
                                          grp_ind = mod_kmlr$clusters
                                          # train_class = mod_kmlr$clusters,
                                          # svcm = F
    )
    
    mse_y_kmlr[f] = mean((Y_train - Y_pred_kmlr$y_pred)^2, na.rm = T)
    
    # cat("MSE (lasir): ", mse_y_lasir[f], "\n")
    # cat("MSE (SVCM): ", mse_y_svcm[f], "\n")
    cat("MSE (KMLR): ", mse_y_kmlr[f], "\n")
    
    
  })
  
  
  
}

save(list = "mse_y_kmlr",
     file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_kmlr_cross_validation.RData")
# 
# # 
# cv_mse_bind = cbind(mse_y_lasir, mse_y_kmlr, mse_y_svcm)
# colnames(cv_mse_bind) = c("LASIR", "KMLR", "SVCM")
# boxplot(cv_mse_bind, ylab ="Mean square error of predicted outcome on test-set",  cex.lab = 1.5, cex.axis = 1.5)
