setwd("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD")

model_list_res_true <- readRDS("./EM_RDA/section2-EM-run/model_list_res.rds")
BIC_list_res_true <- readRDS("./EM_RDA/section2-EM-run/BIC_list_res.rds")

best_idx_true = which(BIC_list_res_true == min(BIC_list_res_true), arr.ind = T)
model_best_true = model_list_res_true[[best_idx_true[1]]][[best_idx_true[2]]]



rho = c("02","05","08")
polydegree = 11:15

NMI_with_true = matrix(NA, nrow = length(rho), length(polydegree))
for(i in 1:length(rho)){
  r = rho[i]
  for(j in 1:length(polydegree)){
    
    p = polydegree[j]
    
    if(p == 14 & r == "05"){
      NMI_with_true[i,j] = 1
      next()
    } 
    filename = file.path("./EM_RDA/section5-sensitivity/", sprintf("result_p%d_rho%s.RData",p, r))
    load(filename)
    BIC_list_res = get(sprintf("BIC_list_p%d_r%s", p, r))
    model_list_res = get(sprintf("model_list_p%d_r%s", p, r))
    # Obtain best model fitting with minimum BIC
    best_idx = which(BIC_list_res == min(BIC_list_res), arr.ind = T)
    model_best = model_list_res[[best_idx[1]]][[best_idx[2]]]    
   
    # Compute NMI with the results presented in draft
    NMI_with_true[i,j] = aricode::NMI(model_best_true$cluster, model_best$cluster) 
  }
}
