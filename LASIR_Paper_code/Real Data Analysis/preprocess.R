##########################################################
#                    Read in image data                  #
##########################################################

load(file = "./ABCD_3mm/subjects_with_task/imageDat_mat_4D_2bk_0bk_realigned.RData")
imageDat_mat_2bk_0bk <- imageDat_mat_2bk_0bk_realigned

# First 3 dimension: x, y, z axis
# Last dimension: number of individual N

N = dim(imageDat_mat_2bk_0bk)[4]
Vx = dim(imageDat_mat_2bk_0bk)[1]
Vy = dim(imageDat_mat_2bk_0bk)[2]
Vz = dim(imageDat_mat_2bk_0bk)[3]


# first assume that all voxels can be removed due to NaN value
vox_to_rm = array(TRUE, dim = c(Vx, Vy, Vz))


for(i in 1:N){
  # Once an image contains a non-Nan value, the flag will be turned to FALSE
  # Meaning that it can not be removed due to the information it contains
  vox_to_rm = vox_to_rm & is.nan(imageDat_mat_2bk_0bk[,,,i])
  
}
vox_eff = !vox_to_rm


################## Checking for missingness ######################
eff_voxel_nan_count = array(dim = dim(vox_eff))
eff_voxel_var = array(dim = dim(vox_eff))
for(i in 1:Vx){
  for(j in 1:Vy){
    for(k in 1:Vz){
      if(vox_to_rm[i,j,k]) next() # if this specific voxel is missing for all individuals, we don't need to count it 
      eff_voxel_nan_count[i,j,k] = mean(is.nan(imageDat_mat_2bk_0bk_realigned[i,j,k,]))
      if(eff_voxel_nan_count[i,j,k]<=0.70){
        eff_voxel_var[i,j,k] = var(imageDat_mat_2bk_0bk_realigned[i,j,k,], na.rm = T)
      }
    }
  }
}

vox_eff[eff_voxel_nan_count>0.30] = FALSE
std_img = readNIfTI(fname = "./ABCD_3mm/subjects_with_task/0bk-baseline_nii/sub-NDARINV00X2TBWJ_0bk-baseline_con_3mm.nii")
slot(std_img, ".Data") = eff_voxel_var
oro.nifti::writeNIfTI(std_img, "./EM_RDA/section1-preprocess/voxel_var")
slot(std_img, ".Data") = eff_voxel_nan_count
oro.nifti::writeNIfTI(std_img, "./EM_RDA/section1-preprocess/voxel_nan")




############################################
imageDat_mat_eff = matrix(NA, nrow = N, ncol = sum(vox_eff))
for(i in 1:N){
  imageDat_mat_eff[i,] = imageDat_mat_2bk_0bk[,,,i][vox_eff]
}


imageDat_mat_eff[is.na(imageDat_mat_eff)]  = 0


## Voxel location matrix, should be a V*3
register_vox = function(vec, center, length_scale){
  (vec - center)/length_scale 
}

loc = t(apply(which(vox_eff, arr.ind = T), 1, function(x){
  register_vox(x, center = c(61, 73, 61)/2, length_scale=Vy)
}))



## Predictors - covariate of interest
library(dplyr)
data_full  = read.csv("./Biostats_nPack_varsInterest.csv")
data_full  = fastDummies::dummy_cols(data_full,
                                     select_columns = "RaceEthnicity",
                                     remove_selected_columns = T)

## Transform age into year
data_full  = data_full %>% mutate(Age = floor(Age)/12)
data_full = fastDummies::dummy_cols(data_full, select_columns = "Age", remove_selected_columns = T, remove_first_dummy = T)
names(data_full)[which(names(data_full) == "Age_0.833333333333333")] = "Age_10"

## Create dummy variable for fam_size
data_full = fastDummies::dummy_cols(data_full, select_columns = "fam_size", remove_selected_columns = T, remove_first_dummy = T)
data_full = data_full %>% 
  mutate("fam_size_2+" = fam_size_2 + fam_size_3 + fam_size_5) %>%
  dplyr::select(-c("fam_size_2", "fam_size_3", "fam_size_5")) %>%
  mutate(`HighestParentalEducation_..HS.Diploma` = `HighestParentalEducation_..HS.Diploma`+`HighestParentalEducation_Other`)
  


fixed_pred = c("fam_size_2+", "Age_10", "Female", "HouseholdMaritalStatus",
               "HighestParentalEducation_HS.Diploma.GED", 
               "HighestParentalEducation_Post.Graduate.Degree",
               "HighestParentalEducation_Some.College",
               "HighestParentalEducation_Bachelor",
               "HouseholdIncome_...100K.", 
               "HouseholdIncome_...50K....100K.",
               "RaceEthnicity_Hispanic",
               "RaceEthnicity_White",
               "RaceEthnicity_Black")

# Sites will be the site-indicator vector for every individual, this will be a size-N vector
Sites     = data_full %>% dplyr::select(sort(c("Site_site14", "Site_site02", "Site_site06", "Site_site13", "Site_site08", 
                                               "Site_site15", "Site_site20", "Site_site04", "Site_site11", "Site_site18",
                                               "Site_site05", "Site_site16", "Site_site09", "Site_site21", "Site_site03", 
                                               "Site_site12", "Site_site10", "Site_site07", "Site_site01", "Site_site17",
                                               "Site_site19")))
Sites     = apply(Sites, 1, function(x) which(x == 1))





x_pred   = c("fam_size_2+", "Age_10", "Female", "HouseholdMaritalStatus",
             "HighestParentalEducation_HS.Diploma.GED", 
             "HighestParentalEducation_Post.Graduate.Degree",
             "HighestParentalEducation_Some.College", 
             "HighestParentalEducation_Bachelor",
             "HouseholdIncome_...100K.", 
             "HouseholdIncome_...50K....100K.",
             "RaceEthnicity_Hispanic",
             "RaceEthnicity_White",
             "RaceEthnicity_Black")

x_g   = c("g")

###########################################################

library(MASS)
input_path_dir = "./EM_RDA/input_dat/"
X_demogrphic   = data_full[,x_pred]
X_g           = data_full[, x_g]

if(is.null(dim(X_g))) X_g = matrix(X_g, N, 1)

Z              = cbind(1, data_full[,x_pred])
Q              = ncol(Z) - 1
J              = ncol(X_g)


# clean up missing value
sapply(Z, function(x) sum(is.na(x)))
apply(X_g,2, function(x) sum(is.na(x)))

## all missing value located in marital status, imputed as 0
Z[is.na(Z)] = 0

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

X_demogrphic = apply(X_demogrphic,2, function(x){
  x[is.na(x)] = getmode(x)
  return(x)
})

X_g = apply(X_g,2, function(x){
  x[is.na(x)] = mean(x,na.rm=T)
  return(x)
})





dat = list(N = N,
           nG = length(unique(Sites)),
           Q = ncol(Z) - 1,
           J = ncol(X_g),
           Y = imageDat_mat_eff,
           X_demogrphic = X_demogrphic,
           X_g = X_g,
           Z = Z,
           G = Sites,
           loc = loc,
           vox_eff = vox_eff)



saveRDS(dat, file = "./EM_RDA/input_dat/input_2_0_thres08.rds")
