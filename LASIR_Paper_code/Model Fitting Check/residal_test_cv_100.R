# Cross validation wtih residual test
library(lasir)
library(flexmix)
library(BayesGPfit)

Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/scripts/fit_predict.cpp")
source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/section2-EM-run/kmeans_gp_fit.R")
input_dat = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")

load("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/section4-cross-validation/mse_cross_validation.RData")

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
#mse_y_lasir = rep(NA, n_fold)

eigens_trans = eigen.transform(input_dat$Y, input_dat$loc, smoothness = 0.05, poly_degree = 12)

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

    
    # ---------- LASIR ---------- #
    mod =  image_GPfit(
      Y = Y_train,
      X = X_gp_train,
      grids = input_dat$loc,
      d = 3,
      poly_degree = 12,
      random_var = G_train,
      ncomp = 2:6,
      smoothness = 0.05,
      step.nrep = 5,
      Z = X_demo_train,
      concomitant_type = "mult",
      criteria = "ICL",
      eigens_trans = eigens_trans_f,
      verbose = T,
      intercept = F,
      control = list(
        minprior = 0,
        tol = 1e-4,
        iter.max = 100,
        classify = "SEM"
      )
    )
    # -------------------------- #



    
    
    # ------------------------- Testing -------------------------- #
    fold_ind_test = folds_50_test[[f]]
    Y_test = input_dat$Y[fold_ind_test,]
    X_gp_test = matrix(input_dat$X_g[fold_ind_test,], ncol = 1)
    X_demo_test = input_dat$X_demogrphic[fold_ind_test,]
    G_test  = input_dat$G[fold_ind_test]
    
    
    # -------- Testing of LASIR --------- #

    mod_lasir = mod$model_best
    nK    = length(mod$model_best@size)
    N     = length(fold_ind_test)
    # Probability allocation
    prob_mat = matrix(NA, N, nK)

    for(i in 1:N){
      prob_mat[i,] = exp(c(1, X_demo_test[i,]) %*% mod_lasir@concomitant@coef)
      prob_mat[i,] = prob_mat[i,]/sum(prob_mat[i,])
    }
    grp_ind = apply(prob_mat, 1, function(x) sample(nK, 1, prob = x))
    Y_pred = predict_gp_known_latent(x = cbind(1, X_gp_train), 
                        z = cbind(1, X_demo_train),
                        G = G_train,
                        fixed_coef  = mod$fixed_coef,
                        main_coef  = mod$main_coef,
                        gamma_coef = mod$random_coef,
                        # prob_mat = prob_mat,
                        # y_train = Y_train,
                        grp_ind = mod$clusters
                        # train_class = mod$clusters,
                        # svcm = F
                        )

    mse_y_lasir[f] = mean((Y_train - Y_pred$y_pred)^2)


    cat("MSE (lasir): ", mse_y_lasir[f], "\n")
    
    
  })

  
  
}

save(list = "mse_y_lasir",
     file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_lasir_cross_validation.RData")

# save(list = "mse_y_svcm",
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_svcm_cross_validation.RData")
# 
# 
# save(list = "mse_y_kmlr",
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_kmlr_cross_validation.RData")

# save(list = c("mse_y_lasir","mse_y_svcm", "mse_y_kmlr", "folds_50_train", "folds_50_test"),
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section4-cross-validation/mse_cross_validation.RData")
cv_mse_bind = data.frame(PredictiveMSE = c(mse_y_lasir,mse_y_svcm, mse_y_lasir_shuffle),
                         Method = factor(rep(c("Within Subgroup",  "W/O Subgroup", "Randomized Subgroup"), 
                                        each = 50)))

g <- ggplot(cv_mse_bind, aes(Method, PredictiveMSE))
g + geom_boxplot(aes(fill=Method)) + theme_bw()+
  theme(axis.text.x = element_text(angle=0, vjust=1, size = 14),
        axis.title = element_text(size = 16,vjust=1),
        axis.text.y = element_text(angle=0, vjust=1, size = 14),
        title = element_text(size = 18), 
        legend.text = element_text(size = 12)) + 
  labs(title="Predictive MSE Given by Each Approach", 
       y="Predictive MSE of Testing Set",
       x="")



boxplot(cv_mse_bind, ylab ="Predictive Mean Square Error (MSE)",  cex.lab = 1.5, cex.axis = 1.5)
