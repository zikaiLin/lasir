input_dat = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")

createSubSamples_cluster = function(N, clusters, n_samps = 50, size_prop = 0.9){
  # Combine sites with small number of 
  
  folds_g = list()
  
  cluster_cat = unique(clusters)
  cluster_size = table(clusters)
  nK = length(cluster_size)
  
  n_each_fold = floor(N*size_prop)
  n_each_fold_by_clust = cluster_size * size_prop
  
  for(iter in 1:n_samps){
    fold_iter = unlist(lapply(1:nK, function(x){
      sample(which(clusters == x), n_each_fold_by_clust[x])
    } ))
    names(fold_iter) = NULL
    folds_g[[iter]] = fold_iter
  }
  return(folds_g)
}

folds_50_train = createSubSamples_cluster(input_dat$N, mod$cluster, n_samps = 50, size_prop = 0.9)
folds_50_test = lapply(folds_50_train, function(x){
  setdiff(1:input_dat$N, x)
})

save(list = c("folds_50_test","folds_50_train"), file = "./EM_RDA/section4-cross-validation/folds.RData")