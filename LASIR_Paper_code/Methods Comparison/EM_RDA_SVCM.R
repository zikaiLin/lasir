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
# Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
# Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")
# #source("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/updated_imgGPfit.R")



#####################################################################
#                         EM Classification                         #
#####################################################################
library(lasir)
Y.eigen.trans = eigen.transform(
  Y = input_dat$Y,
  grids = input_dat$loc,
  a = 0.0001,
  smoothness = 0.05,
  poly_degree = 14)


source("./mvEMfix_r.R")
model_list = list()
BIC_list = matrix(Inf, 1, 20)
model_list[[k]] = list()
  for(iter in 1:20){
    cat("* ")
    EM <- mvnMixSEMFix$new(outcome= Y.eigen.trans$Y_star,
                           exposed_var = input_dat$X_g,
                           control_var = input_dat$X_demogrphic,
                           sites = input_dat$G,
                           num_components = as.integer(1))
    
    res_EM <- EM$run.EM(loglik_tol=1e-4)
    BIC_list[k,iter] = res_EM$BIC
    model_list[[k]][[iter]] = res_EM
    
  }
cat("\n")

saveRDS(model_list, "./EM_RDA/section2-EM-run/model_list_res.rds")
saveRDS(BIC_list, "./EM_RDA/section2-EM-run/BIC_list_res.rds")
# list_of_vars_tosave = c("m1")
# 
# save(list = list_of_vars_tosave,
#      file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/result_p13_rho05_m1.RData")


