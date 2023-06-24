# p-value correction inference

input_dat = readRDS("./EM_RDA/Archive/RDA_Apr4/section1-preprocess/input_2_0_thres08.rds")
BIC_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/BIC_list_res_no_intercept.rds")
model_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/model_list_res_no_intercept.rds")
transform_res_p14_rho05 <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/transform_res_p14_rho05.rds")

best_mod_idx = which(BIC_list_res == min(BIC_list_res), arr.ind = T)
mod = model_list_res[[best_mod_idx[1]]][[best_mod_idx[2]]]
nK = length(unique(mod$cluster))
L = ncol(mod$fixed_basis_coef)
J = input_dat$J
nG = input_dat$nG
Q = input_dat$Q
X = input_dat$X_g
Z = input_dat$X_demogrphic
G = input_dat$G

mod_cluster = mod$cluster
n_pars_alpha = (J+1) *L 

observed_info_mat = array(0, dim = c(n_pars_alpha, n_pars_alpha, nK))

inv_observed_info_mat = array(0, dim = c(n_pars_alpha, n_pars_alpha, nK))
for(k in 1:nK){
  inv_observed_info_mat[,,k] = solve(kronecker( crossprod(cbind(1, X[mod$cluster == k,])), mod$sigma_mat))
}

# Transform to volumetric space
V = ncol(input_dat$Y)


var_est_main_effect = array(NA, dim = c(2, V, nK))


for(k in 1:nK){
  for(j in 1:(J+1)){
    subblock_idx = ((j-1)*L+1):(j*L)
    inv_obs_info_mat_subblock = inv_observed_info_mat[subblock_idx, subblock_idx,k]
    var_est_main_effect[j, ,k] = rowSums(transform_res_p14_rho05$svd$u * t(inv_obs_info_mat_subblock %*% t(transform_res_p14_rho05$svd$u)) )
  }
}

sd_est_main_effect = sqrt(var_est_main_effect)



# Transform estimation to volumetric space (coefficient)

main_effect_coef =  array(NA, dim = c(J+1, V, nK))

for(k in 1:nK){
  main_effect_coef[,,k] = mod$main_basis_coef[,,k] %*% t(transform_res_p14_rho05$svd$u)
}


# Obtain p-value estimation

pval_est_main_effect = (1-pnorm(abs((main_effect_coef/sd_est_main_effect)),0,1) )*2 # tail z-test
pval_est_main_effect.adjust = pval_est_main_effect
# P.value adjustment

for(j in 1:(J+1)){
  for(k in 1:nK){
    pval_est_main_effect.adjust[j,,k] <- p.adjust(pval_est_main_effect[j,,k], method = "BH")
  }
}



# Transform estimation to volumetric space (coefficient)

main_effect_coef =  array(NA, dim = c(J+1, V, nK))

for(k in 1:nK){
  main_effect_coef[,,k] = mod$main_basis_coef[,,k] %*% t(transform_res_p14_rho05$svd$u)
}


# Obtain p-value estimation

pval_est_main_effect = (1-pnorm(abs((main_effect_coef/sd_est_main_effect)),0,1) )*2 # tail z-test


saveRDS(pval_est_main_effect, file = "./EM_RDA/section3-post-process/pval_est_main_effect_no_intercept.rds")
saveRDS(sd_est_main_effect, "./EM_RDA/section3-post-process/est_sd_main_effect_no_intercept.rds")
saveRDS(main_effect_coef, "./EM_RDA/section3-post-process/est_main_effect_no_intercept.rds")
