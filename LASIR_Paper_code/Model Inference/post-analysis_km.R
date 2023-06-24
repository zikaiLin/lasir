rm(list = ls())
load("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/result_p12_rho05_varDF_km.RData")

library(dplyr)

############################################################################
#             Save image                                                   #
############################################################################

#--------------- Read in a new image --------------#
names_X = c("Intercept", "g")
nK = length(unique(m1$clusters))
J = 2
# names_X = c( "Family_Size_2" ,  "Age_10" ,                                      
#              "Female" , "MaritalStatus" , "Parental_Edu_HS_GED"  ,
#                "Parental_Edu_Post_Graduate_Degree", "Parental_Edu_College", "Parental_Edu_Bachelor","Income_100K." ,                    
#              "Income_50K_100K.", "Race_Hispanic" ,"Race_White", "Race_Black")

pval_dir = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section3-post-process/alpha_pval_km"
coef_dir = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section3-post-process/alpha_coef_km"

if(!dir.exists(pval_dir)) dir.create(pval_dir)
if(!dir.exists(coef_dir)) dir.create(coef_dir)

# -------------Write pvalue image ----------------#
library(oro.nifti)
std_img = oro.nifti::readNIfTI("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/ABCD_3mm/ImagingData/img_2bk-baseline_con/sub-NDARINV00CY2MDM_2bk-baseline_con_3mm.nii.gz")
slot(std_img,".Data")[!is.nan(std_img)] = NaN 

yeo_img = oro.nifti::readNIfTI("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/Resources/Parcellations/yeo7_MNI152_3mm.nii.gz")
yeo_img_mask = yeo_img != 0


alpha_img_list_pval = list()

for(j in 1:length(names_X)){
  alpha_img_list_pval[[j]] = list()
  filename = names_X[j]
  
  for(k in  1:nK){
    # For each subgroup first write out the p-value image
    a1_img_sig = std_img
    
    slot(a1_img_sig,".Data")[input_dat$vox_eff] = m1$fixed_pval[j,,k] 
    
    alpha_img_list_pval[[j]][[k]] = a1_img_sig
    filename1_sig = file.path(pval_dir, sprintf("alpha_%s_pval_k%d", filename, k))
    
    oro.nifti::writeNIfTI(alpha_img_list_pval[[j]][[k]], filename=filename1_sig)
  }
}



# Mask with FDR-corrected p-value mark with cortical parcellation image as well
alpha_img_list_sig = list()

for(j in 1:J){
  alpha_img_list_sig[[j]] = list()
  for(k in 1:nK){
    
    filename = filename = names_X[j]
    a1_img = std_img
    
    slot(a1_img,".Data")[input_dat$vox_eff] = m1$fixed_coef[j,, k] 
    
    # correct p-value
    filename1_sig = file.path(pval_dir, sprintf("alpha_%s_pval_k%d", filename, k))
    p_thres = as.numeric(system(sprintf("/usr/local/fsl/bin/fdr -i %s -q 0.0001", filename1_sig), intern = T)[2])
    slot(a1_img,".Data")[alpha_img_list_pval[[j]][[k]]>=p_thres] = NaN
    
    
    
    # mask with cortical mask
    #slot(a1_img,".Data")[!yeo_img_mask] = NaN
    
    
    alpha_img_list_sig[[j]][[k]] = a1_img
    
    filename1_sig_coef = file.path(coef_dir, sprintf("alpha_%s_coef_k%d", filename, k))
    oro.nifti::writeNIfTI(alpha_img_list_sig[[j]][[k]], filename=filename1_sig_coef)
  }
}




# ----------- A script to output significant regions in Power-264 --------- #

library(readxl)
power.264 = read_xlsx("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/Resources/Parcellations/Power264_Consensus264_simplified.xlsx",sheet = 1, col_names = T)
power.264.img = oro.nifti::readNIfTI("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/Resources/Parcellations/power264-master/power264-master/power264MNI_3mm.nii.gz")

alpha_sig_regions = list()

for(j in 1:J){
  alpha_sig_regions[[names_X[j]]] = list()
  for(k in 1:nK){
    alpha_sig_regions[[names_X[j]]][[k]] = list()
    sum_table = power.264[unique(power.264.img[power.264.img >0 & !is.nan(alpha_img_list_sig[[j]][[k]])]),]
    
    if(nrow(sum_table) == 0) next()
    
    sum_table$roi_sign = rep(NA, nrow(sum_table))
    for(roi in 1:length(sum_table$roi_sign)){
      ROI = sum_table$ROI[roi]
      sum_table$roi_sign[roi] = sign(mean(alpha_img_list_sig[[j]][[k]][power.264.img == ROI], na.rm = T))
    }
    
    alpha_sig_regions[[names_X[j]]][[k]]$all_roi = sum_table[,c("ROI", "X","Y","Z", "roi_sign", "Suggested System")]
    
    alpha_sig_regions[[names_X[j]]][[k]]$networks = unique(sum_table[,"Suggested System"])
    
  }
}

save(list = "alpha_sig_regions", file = "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section3-post-process/alpha_sig_regions_rho005_p12_km.RData")




############################################################################
#             Report concomitant coefficients                              #
############################################################################
grp_ind_fct = relevel(factor(m1$clusters), ref = "1")
library(nnet)
fit_mlt = nnet::multinom(grp_ind_fct ~ input_dat$X_demogrphic)
fit_mlt_coef = coef(fit_mlt)
colnames(fit_mlt_coef) = c("Intercept", colnames(input_dat$X_demogrphic))
fit_mlt_sd   = summary(fit_mlt)$standard.errors
fit_mlt_tval = abs(fit_mlt_coef)/fit_mlt_sd
fit_mlt_pval = (1-stats::pnorm(fit_mlt_tval,0,1))*2
fit_mlt_pval
# xtable::xtable(t(fit_mlt_pval), digits = 3)





############################################################################
#             Random effect clustering and report                          #
############################################################################

nK = length(unique(m1$clusters))
sites_table = matrix(NA, 21, nK)

for(g in 1:21){
  sites_count =  table(m1$clusters[input_dat$G == g])
  for(k in 1:nK){
    if(!(as.character(k) %in% names(sites_count))){
      # Set as 0
      sites_table[g,k] = 0
    }else{
      sites_table[g,k] = sites_count[as.character(k)]
    }
  }
}

# Perform chisquare test for pairwise combination
chis_sq_table = matrix(NA, 21, 21)
sites_comb = gtools::combinations(21, 2)

for(s_comb in 1:nrow(sites_comb)){
  try({chisq.test_sites = fisher.test(rbind(sites_table[sites_comb[s_comb,1],], sites_table[sites_comb[s_comb,2],]), simulate.p.value = T)
  if(chisq.test_sites$p.value < 0.05){
    chis_sq_table[sites_comb[s_comb,1], sites_comb[s_comb,2]]  = chis_sq_table[sites_comb[s_comb,2], sites_comb[s_comb,1]] = 1
  }
  })
}

library(corrplot)
library(ape)
corrplot(chis_sq_table, method = "square", type = "upper", na.label = " ")



chis_sq_table_1 = chis_sq_table
rownames(chis_sq_table_1) = colnames(chis_sq_table_1) = as.character(1:21)
chis_sq_table_1[is.na(chis_sq_table_1)] = 0
rm_site = which(apply(sites_table, 1, sum) < 15)
rownames(chis_sq_table_1) = colnames(chis_sq_table_1) = paste("site", (1:21))
hc = hclust(as.dist(chis_sq_table_1[-rm_site, -rm_site]))

colors = c("red", "blue", "black")
clus4 = cutree(hc, 3)
plot(as.phylo(hc), type = "fan", tip.color = colors[clus4], edge.width = 5,
     use.edge.length = T, cex = 2.5, no.margin = F)



site_cluster1 = c(18,13,4,8)
site_cluster2 = c(21,20,14,12,11,9,2,5)
site_cluster3 = c(16,3,15)
site_cluster4 = c(10,6)

sites_table[site_cluster1,] # skew to the largest subgroup
site_cluster1_random_efect = apply(m1$random_coef[site_cluster1, ], 2, mean)
site_cluster2_random_efect = apply(m1$random_coef[site_cluster2, ], 2, mean)
site_cluster3_random_efect = apply(m1$random_coef[site_cluster3, ], 2, mean)
site_cluster4_random_efect = apply(m1$random_coef[site_cluster4, ], 2, mean)


# Write nifti image
library(oro.nifti)
std_img = oro.nifti::readNIfTI("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/ABCD_3mm/ImagingData/img_2bk-baseline_con/sub-NDARINV00CY2MDM_2bk-baseline_con_3mm.nii.gz")
slot(std_img,".Data")[!is.nan(std_img)] = NaN 


site_cluster1_random_efect_img = site_cluster2_random_efect_img= site_cluster3_random_efect_img = site_cluster4_random_efect_img = std_img
slot(site_cluster1_random_efect_img, ".Data")[input_dat$vox_eff] = site_cluster1_random_efect
slot(site_cluster2_random_efect_img, ".Data")[input_dat$vox_eff] = site_cluster2_random_efect
slot(site_cluster3_random_efect_img, ".Data")[input_dat$vox_eff] = site_cluster3_random_efect
slot(site_cluster4_random_efect_img, ".Data")[input_dat$vox_eff] = site_cluster4_random_efect

writeNIfTI(site_cluster1_random_efect_img, "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/result/site_cluster1_random_efect_img")
writeNIfTI(site_cluster2_random_efect_img, "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/result/site_cluster2_random_efect_img")
writeNIfTI(site_cluster3_random_efect_img, "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/result/site_cluster3_random_efect_img")
writeNIfTI(site_cluster4_random_efect_img, "/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/result/site_cluster4_random_efect_img")






############################################################################
#                    residual diagnostic                                   #
############################################################################
Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/scripts/fit_predict.cpp")
mod_residuals = residual.gp.em(
  input_dat$Y,
  cbind(1,input_dat$X_g),
  cbind(1,input_dat$X_demogrphic),
  grp_ind = m1$clusters,
  alpha_coef_post = m1$fixed_coef,
  gamma_post_all = m1$random_coef[input_dat$G,],
  main_coef = m1$fixed_main_coef
)


nK    = length(m1$model_best@size)

Y_pred = predict_gp(x = cbind(1, input_dat$X_g), z = cbind(1, input_dat$X_demogrphic), G = input_dat$G,
                    main_coef = m1$fixed_main_coef,
                    alpha_coef = m1$fixed_coef,
                    gamma_coef = m1$random_coef,
                    grp_ind = m1$clusters)

residuals_em = input_dat$Y - Y_pred$y_pred
std_resid_em = residuals_em
V = ncol(input_dat$Y)
dnorm_test = dnorm(rnorm(input_dat$N))
KL_div_resid_vox_em = rep(NA, V)
for(v in 1:V){
  # The off-diagnonal will be the KL divergence between standardized residual and std-normal
  std_resid_em[,v] = residuals_em[,v]/sd(residuals_em[,v])
  KL_div_resid_vox_em[v] = LaplacesDemon::KLD(dnorm_test, dnorm(std_resid_em[,v]))$mean.sum.KLD
}

library(oro.nifti)
std_img = oro.nifti::readNIfTI("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/ABCD_3mm/ImagingData/img_2bk-baseline_con/sub-NDARINV00CY2MDM_2bk-baseline_con_3mm.nii.gz")
slot(std_img,".Data")[!is.nan(std_img)] = NaN 
cutoff_KL = quantile(KL_div_resid_vox_em, 0.95, na.rm = T)
KL_div_resid_vox_em[KL_div_resid_vox_em<cutoff_KL] = NaN
std_img[input_dat$vox_eff] = KL_div_resid_vox_em
oro.nifti::writeNIfTI(std_img, sprintf("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section3-post-process/KL_div_RDA"))