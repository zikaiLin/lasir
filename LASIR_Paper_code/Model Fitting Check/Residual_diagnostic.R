#  --------------------- Residual Diagnostic and comparison with KMLR --------- #
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD//EM_RDA/scripts/fit_predict.cpp")

input_dat = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")

# --------------------------------------------------------------------------
# Step 1: Create a standard normal distribution sample for comparison of KL divergence
# --------------------------------------------------------------------------


# dnorm_sd = dnorm(sort(input_dat$N))
dnorm_sd = dnorm(sort(rnorm(2021)))
rnorm_sd = rnorm(2021)
# --------------------------------------------------------------------------
# Step 2(a): Use existing fit of LASIR to calculate the residual diagnostic.
# --------------------------------------------------------------------------

transform_res_p14_rho05 <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/transform_res_p14_rho05.rds")
BIC_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/BIC_list_res.rds")
model_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/model_list_res.rds")
best_mod_idx = which(BIC_list_res == min(BIC_list_res), arr.ind = T)
# mod = model_list_res[[best_mod_idx[1]]][[best_mod_idx[2]]]
mod = model_list_res[[4]][[3]]

# Transform estimation to volumetric space (coefficient)
nK    = length(table(mod$cluster))
V = ncol(input_dat$Y)

main_effect_coef =  array(NA, dim = c(input_dat$J+1, V, nK))

for(k in 1:nK){
  main_effect_coef[,,k] = mod$main_basis_coef[,,k] %*% t(transform_res_p14_rho05$svd$u)
}

site_coef = mod$sites_basis_coef %*% t(transform_res_p14_rho05$svd$u)
fixed_coef = mod$fixed_basis_coef %*% t(transform_res_p14_rho05$svd$u)




Y_pred = predict_gp_known_latent(x = cbind(1, input_dat$X_g),
                                 z = input_dat$X_demogrphic, 
                                 G = input_dat$G,
                                 fixed_coef = fixed_coef,
                                 main_coef = main_effect_coef,
                                 gamma_coef = site_coef,
                                 grp_ind = mod$cluster)

residuals_y_lasir = input_dat$Y- Y_pred$y_pred
residuals_y_sd_lasir = apply(residuals_y_lasir, 2, function(x) x/sd(x))
qqnorm(residuals_y_sd_lasir[,20000])

KLdiv_lasir = rep(NA, nrow(input_dat$loc))

for(v in 1:nrow(input_dat$loc)){
  KLdiv_lasir[v] = flexmix::KLdiv(cbind(dnorm(sort(residuals_y_sd_lasir[,v])), dnorm_sd))[1,2]
}

saveRDS(KLdiv_lasir, file = "./EM_RDA/section4-cross-validation/KLdiv_lasir.rds")


# --------------------------------------------------------------------------
# Step 2(b): fit a SVCM with only 1 group
# --------------------------------------------------------------------------
transform_res_p14_rho05 <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/transform_res_p14_rho05.rds")
BIC_list_res_svcm <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/BIC_list_svcm_res.rds")
model_list_res_svcm <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/model_list_svcm_res.rds")
best_mod_idx_svcm = which(BIC_list_res_svcm == min(BIC_list_res_svcm), arr.ind = T)
mod_svcm = model_list_res_svcm[[best_mod_idx_svcm[1]]][[best_mod_idx_svcm[2]]]


# Transform estimation to volumetric space (coefficient)
nK    = 1
V = ncol(input_dat$Y)

main_effect_coef_svcm =  array(NA, dim = c(input_dat$J+1, V, nK))

for(k in 1:nK){
  main_effect_coef_svcm[,,k] = mod_svcm$main_basis_coef[,,k] %*% t(transform_res_p14_rho05$svd$u)
}

site_coef_svcm = mod_svcm$sites_basis_coef %*% t(transform_res_p14_rho05$svd$u)
fixed_coef_svcm = mod_svcm$fixed_basis_coef %*% t(transform_res_p14_rho05$svd$u)



grp_ind = rep(1, input_dat$N)

Y_pred = predict_gp_known_latent(x = cbind(1, input_dat$X_g),
                                 z = input_dat$X_demogrphic, 
                                 G = input_dat$G,
                                 fixed_coef = fixed_coef_svcm,
                                 main_coef = main_effect_coef_svcm,
                                 gamma_coef = site_coef_svcm,
                                 grp_ind = grp_ind)


residuals_y_svcm = input_dat$Y- Y_pred$y_pred
residuals_y_sd_svcm = apply(residuals_y_svcm, 2, function(x) x/sd(x))



KLdiv_svcm = rep(NA, nrow(input_dat$loc))

for(v in 1:nrow(input_dat$loc)){
  KLdiv_svcm[v] = flexmix::KLdiv(cbind(dnorm(sort(residuals_y_sd_svcm[,v])), dnorm_sd))[1,2]
}
saveRDS(KLdiv_svcm, file = "./EM_RDA/section4-cross-validation/KLdiv_svcm.rds")


# --------------------------------------------------------------------------
# Step 2(c): fit a KMLR with multiple group
# --------------------------------------------------------------------------
load("./EM_RDA/section2-EM-run/result_p14_rho05_KM.RData")

nK    = length(unique(m1$clusters))
grp_ind = m1$clusters

Y_pred = predict_gp_known_latent(x = cbind(1, input_dat$X_g),
                    z = input_dat$X_demogrphic,
                    G = input_dat$G,
                    fixed_coef = m1$fixed_coef,
                    main_coef = m1$main_coef,
                    gamma_coef = m1$site_coef,
                    grp_ind = grp_ind)

residuals_y_km = input_dat$Y- Y_pred$y_pred
residuals_y_sd_km = apply(residuals_y_km, 2, function(x) x/sd(x))
KLdiv_kmlr = rep(NA, nrow(input_dat$loc))

for(v in 1:nrow(input_dat$loc)){
  KLdiv_kmlr[v] =  flexmix::KLdiv(cbind(dnorm(sort(residuals_y_sd_km[,v])), dnorm_sd))[1,2]
}
saveRDS(KLdiv_kmlr, file = "./EM_RDA/section4-cross-validation/KLdiv_kmlr.rds")

KLdiv_bind = cbind(KLdiv_lasir, KLdiv_kmlr, KLdiv_svcm)
colnames(KLdiv_bind) = c("LASIR", "KMLR", "SVCM")
boxplot(KLdiv_bind, ylab ="Kullback–Leibler divergence on voxels",  cex.lab = 1.5, cex.axis = 1.5)


######### Use ggplot for qqplot ###############

  sum((KLdiv_lasir < KLdiv_svcm) & (KLdiv_lasir < KLdiv_kmlr))

# library(ggplot2)
# df <- data.frame(y = residuals_y_lasir[,v])
# p <- ggplot(df, aes(sample = y))
# p + stat_qq() + stat_qq_line() + theme_bw()
# 
# 
# library(ggplot2)
# df <- data.frame(y = residuals_y_svcm[,v])
# p <- ggplot(df, aes(sample = y))
# p + stat_qq() + stat_qq_line() + theme_bw()
# 

v = 6083

v = 20000

library(ggplot2)
df <- data.frame(Residuals = c(rnorm(2021),residuals_y_sd_lasir[,v], residuals_y_sd_svcm[,v], residuals_y_sd_km[,v]), Method = rep(c( "Standard Normal", "LASIR", "SVCM", "KMLR"), each = 2021))
p <- ggplot(df, aes(sample = Residuals, colour = Method))
p + stat_qq() + stat_qq_line() + theme_bw()

p2 <- ggplot(df, aes(x = Residuals, fill = Method, alpha = I(0.4)))
p2 + geom_histogram(bins = 100, position = "identity")
