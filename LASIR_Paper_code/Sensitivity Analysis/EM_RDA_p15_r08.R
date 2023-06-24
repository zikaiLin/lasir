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


library(lasir)
rho = 8
poly_degree = 15
if(file.exists(sprintf("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section5-sensitivity/Y_eigen_trans_p%d_r0%d.rds", poly_degree, rho))){
  Y.eigen.trans = readRDS(sprintf("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section5-sensitivity/Y_eigen_trans_p%d_r0%d.rds", poly_degree, rho))
}else{
  Y.eigen.trans = eigen.transform(
    Y = input_dat$Y,
    grids = input_dat$loc,
    concentration = 0.0001,
    smoothness = 0.08,
    poly_degree = 15)
  
  saveRDS(Y.eigen.trans, "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section5-sensitivity/Y_eigen_trans_p15_r08.rds")
  
  
}

#####################################################################
#                         EM Classification                         #
#####################################################################


source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/mvEMfix_r.R")
model_list_p15_r08 = list()
BIC_list_p15_r08 = matrix(Inf, 6, 5)
for(k in 1:6){
  cat("k = ", k, "\n")
  model_list_p15_r08[[k]] = list()
  for(iter in 1:5){
    cat("* ")
    res_EM <- lasir::LASIR.fit(Y.eigen.trans$Y_star,
                               input_dat$X_g,
                               input_dat$X_demogrphic,
                               sites = input_dat$G,
                               nK = as.integer(k))
    # EM <- mvnMixSEMFix$new(outcome= Y.eigen.trans$Y_star,
    #                        exposed_var = input_dat$X_g,
    #                        control_var = input_dat$X_demogrphic,
    #                        sites = input_dat$G,
    #                        num_components = as.integer(k))
    # 
    # res_EM <- EM$run.EM(loglik_tol=1e-4)
    BIC_list_p15_r08[k,iter] = res_EM$BIC
    model_list_p15_r08[[k]][[iter]] = res_EM
    
  }
  cat("\n")
}


list_of_vars_tosave = c("model_list_p15_r08", "BIC_list_p15_r08")

save(list = list_of_vars_tosave,
     file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section5-sensitivity/result_p15_rho08.RData")

