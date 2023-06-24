#set.seed(2020)
library(flexmix)
library(gtools)
library(lasir)
rm(list = ls())
####### Input ##############
input_dat <-
  readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")



#####################################################################
#      Get initial estimate of constant main effect                 #
#####################################################################

# 
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_utils.cpp")
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_em_descent_update.cpp")




#####################################################################
#                         EM Classification                         #
#####################################################################


m1 = image_GPfit(
  Y = input_dat$Y,
  X = input_dat$X_g,
  grids = input_dat$loc,
  d = 3,
  poly_degree = 12,
  random_var = input_dat$G,
  ncomp = 1:10,
  smoothness = 0.08,
  step.nrep = 5,
  Z = input_dat$X_demogrphic,
  concomitant_type = "mult",
  criteria = "ICL",
  verbose = T,
  intercept = F,
  control = list(
    minprior = 0,
    tol = 1e-4,
    iter.max = 50,
    classify = "SEM"
  )
)


list_of_vars_tosave = c("m1")

save(list = list_of_vars_tosave,
     file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section5-sensitivity/result_p12_rho08.RData")





# sens_size_tab = rbind(c(1493,267,261),
#                       c(1480,296,245),
#                       c(1464,315,242),
#                       c(1509,262,250),
#                       c(1474,307,240),
#                       c(1484,298,239),
#                       c(1486,287,248),
#                       c(1479,290,252),
#                       c(1468,307,246))
