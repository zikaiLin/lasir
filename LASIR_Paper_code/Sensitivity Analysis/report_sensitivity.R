############################################################################
#             section5-sensitivity                                        #
############################################################################
rho = c("01", "05", "08")
poly = c("11", "12", "13")
sens_cluster_tab = matrix(NA, 9, 3)
i = 1
for(r in rho){
  for (p in poly) {
    fileName.res = sprintf("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section5-sensitivity/result_p%s_rho%s.RData", p, r)
    load(fileName.res)
    
    sens_cluster_tab[i, ] = sort(m1$model_best@size)
    i = i+1
    
    
    #--------------- Read in a new image --------------#
    names_X = c("Intercept", "g")
    nK = length(unique(m1$clusters))
    J = length(names_X) 
    # names_X = c( "Family_Size_2" ,  "Age_10" ,                                      
    #              "Female" , "MaritalStatus" , "Parental_Edu_HS_GED"  ,
    #                "Parental_Edu_Post_Graduate_Degree", "Parental_Edu_College", "Parental_Edu_Bachelor","Income_100K." ,                    
    #              "Income_50K_100K.", "Race_Hispanic" ,"Race_White", "Race_Black")
    
    
    # -------------Write pvalue image ----------------#
    home_dir = sprintf("EM_RDA/section5-sensitivity/rho%s_p%s/", r, p)
    pval_dir = sprintf("EM_RDA/section5-sensitivity/rho%s_p%s/alpha_pval", r, p)
    coef_dir = sprintf("EM_RDA/section5-sensitivity/rho%s_p%s/alpha_coef", r, p)
    dir.create(home_dir)
    dir.create(pval_dir)
    dir.create(coef_dir)
    
    library(oro.nifti)
    std_img = oro.nifti::readNIfTI("./ABCD_3mm/ImagingData/img_2bk-baseline_con/sub-NDARINV00CY2MDM_2bk-baseline_con_3mm.nii.gz")
    slot(std_img,".Data")[!is.nan(std_img)] = NaN 
    
    yeo_img = oro.nifti::readNIfTI("./Resources/Parcellations/yeo7_MNI152_3mm.nii.gz")
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
        p_thres = as.numeric(system(sprintf("/usr/local/fsl/bin/fdr -i %s -q 0.000001", filename1_sig), intern = T)[2])
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
    power.264 = read_xlsx("./Resources/Parcellations/Power264_Consensus264_simplified.xlsx",sheet = 1, col_names = T)
    power.264.img = oro.nifti::readNIfTI("./Resources/Parcellations/power264-master/power264-master/power264MNI_3mm.nii.gz")
    
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
    
    save(list = "alpha_sig_regions", file = file.path(home_dir, "alpha_sig_regions_rho005.RData"))
    
  }
}


# 
# 1    2    3 
# 282 1519  220 
# 1    2    3                                                                                                                               
# 284 1527  210 
# 1    2    3                                                                                                                               
# 1511  226  284 
# 1    2    3                                                                                                                               
# 280  223 1518 
# 1    2    3                                                                                                                               
# 276  223 1522 
# 1    2    3                                                                                                                               
# 267  238 1516 
# 1    2    3                                                                                                                               
# 289  225 1507 
# 1    2    3                                                                                                                               
# 1517  217  287 
# 1    2    3                                                                                                                               
# 276  238 1507 
