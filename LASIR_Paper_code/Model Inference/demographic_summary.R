
library(dplyr)
model_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/model_list_res.rds")
BIC_list_res <- readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section2-EM-run/BIC_list_res.rds")

best_mod_idx = which(BIC_list_res == min(BIC_list_res), arr.ind = T)
mod  = model_list_res[[best_mod_idx[1]]][[best_mod_idx[[2]]]]
input_dat = readRDS("/Volumes/easystore/Dropbox/U-Mich/Research/ABCD/EM_RDA/section1-preprocess/input_2_0_thres08.rds")

data_full_post = as.data.frame(input_dat$X_demogrphic)
data_full_post$grp_ind = factor(mod$cluster)

race_post = list()
for(k in 1:length(unique(mod$cluster))){
  race_post[[k]] = data_full_post %>% 
    select(grp_ind, RaceEthnicity_Hispanic, RaceEthnicity_White, RaceEthnicity_Black) %>% 
    filter(grp_ind == k) %>%
    mutate(Race = 1*RaceEthnicity_Black + 2*RaceEthnicity_White + 3*RaceEthnicity_Hispanic) %>%
    mutate(Race = factor(Race))%>%
    mutate(Race = plyr::revalue(Race,replace = c("0"="Other+Asian", "1"="Black","2"="White", "3"="Hispanic")))%>%
    group_by(Race) %>% 
    summarise(Count = n()) %>% mutate(pct = Count / sum(Count) * 100)
}




#############     Education      ####################
edu_post = list()
for(k in 1:length(unique(mod$cluster))){
  edu_post[[k]] = data_full_post %>% 
    select(grp_ind, HighestParentalEducation_HS.Diploma.GED, HighestParentalEducation_Post.Graduate.Degree, HighestParentalEducation_Some.College, HighestParentalEducation_Bachelor) %>% 
    filter(grp_ind == k) %>%
    mutate(Edu = 1*`HighestParentalEducation_HS.Diploma.GED` + 4*`HighestParentalEducation_Post.Graduate.Degree` + 2*`HighestParentalEducation_Some.College`+3*`HighestParentalEducation_Bachelor`) %>%
    mutate(Edu = factor(Edu))%>%
    mutate(Edu = plyr::revalue(Edu,replace = c("0"="High School",
                                               "1"="HS.GED","2"="College", "3"="Bachelors", "4"="Post Graduate")))%>%
    group_by(Edu) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(Edu)) %>%
    mutate(prop = Count / sum(Count) *100) 
}


#############     Income      ####################
income_post = list()
for(k in 1:length(unique(mod$cluster))){
income_post[[k]] = data_full_post %>% 
  select(grp_ind, `HouseholdIncome_...50K....100K.`, `HouseholdIncome_...100K.`) %>% 
  filter(grp_ind == k) %>%
  mutate(income = 1*`HouseholdIncome_...50K....100K.` + 2*`HouseholdIncome_...100K.`) %>%
  mutate(income = factor(income))%>%
  mutate(income = plyr::revalue(income,replace = c("0"="<50K",
                                                   "1"="50K-100K",
                                                   "2"="100K")))%>%
  group_by(income) %>% 
  summarise(Count = n()) %>% 
  arrange(desc(income)) %>%
  mutate(prop = Count / sum(Count) *100) 
}


#############     Income      ####################
income_post = list()
for(k in 1:length(unique(mod$cluster))){
  income_post[[k]] = data_full_post %>% 
    select(grp_ind, `HouseholdIncome_...50K....100K.`, `HouseholdIncome_...100K.`) %>% 
    filter(grp_ind == k) %>%
    mutate(income = 1*`HouseholdIncome_...50K....100K.` + 2*`HouseholdIncome_...100K.`) %>%
    mutate(income = factor(income))%>%
    mutate(income = plyr::revalue(income,replace = c("0"="<50K",
                                                     "1"="50K-100K",
                                                     "2"="100K")))%>%
    group_by(income) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(income)) %>%
    mutate(prop = Count / sum(Count) *100) 
}



#############     Marrital      ####################
married_post = list()
for(k in 1:length(unique(mod$cluster))){
married_post[[k]] = data_full_post %>% 
  select(grp_ind, `HouseholdMaritalStatus`) %>% 
  filter(grp_ind == k) %>%
  mutate(married = factor(HouseholdMaritalStatus))%>%
  mutate(married = plyr::revalue(married,replace = c("0"="Single",
                                                     "1"="Married")))%>%
  group_by(grp_ind, married) %>% 
  summarise(Count = n()) %>% 
  arrange(desc(married)) %>%
  mutate(prop = Count / sum(Count) *100) 
}




#############     Famsize      ####################

famsize_post = list()
for(k in 1:length(unique(mod$cluster))){
  famsize_post[[k]] = data_full_post %>% 
    select(grp_ind, `fam_size_2+`) %>% 
    filter(grp_ind == k) %>%
    mutate(Famsize = factor(`fam_size_2+`))%>%
    mutate(Famsize = plyr::revalue(Famsize,replace = c("0"="<2",
                                                       "1"=">=2")))%>%
    group_by(grp_ind, Famsize) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(Famsize)) %>%
    mutate(prop = Count / sum(Count) *100) 
}


#############     Age      ####################

age_post = list()
for(k in 1:length(unique(mod$cluster))){
  age_post[[k]] = data_full_post %>% 
    select(grp_ind, `Age_10`) %>% 
    filter(grp_ind == k) %>%
    mutate(age = factor(`Age_10`))%>%
    mutate(age = plyr::revalue(age,replace = c("0"="<10",
                                               "1"=">=10")))%>%
    group_by(grp_ind, age) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(age)) %>%
    mutate(prop = Count / sum(Count) *100)
  
}



#############     Sex      ####################

female_post = list()
for(k in 1:length(unique(mod$cluster))){
  female_post[[k]] = data_full_post %>% 
    select(grp_ind, `Female`) %>% 
    filter(grp_ind == k) %>%
    mutate(female = factor(`Female`))%>%
    mutate(female = plyr::revalue(female,replace = c("0"="Male",
                                                     "1"="Female")))%>%
    group_by(grp_ind, female) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(female)) %>%
    mutate(prop = Count / sum(Count) *100) 
  
}


########### Count ###############
race_count = abind::abind(lapply(race_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(race_count) = race_post[[1]]$Race

edu_count = abind::abind(lapply(edu_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(edu_count) = edu_post[[1]]$Edu


income_count = abind::abind(lapply(income_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(income_count) = income_post[[1]]$income

married_count = abind::abind(lapply(married_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(married_count) = married_post[[1]]$married

famsize_count = abind::abind(lapply(famsize_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(famsize_count) = famsize_post[[1]]$Famsize

age_count = abind::abind(lapply(age_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(age_count) = age_post[[1]]$age

female_count = abind::abind(lapply(female_post, FUN = function(x) x[,"Count"]), along = 2)
rownames(female_count) = female_post[[1]]$female


overall_demographic_count = rbind(race_count, edu_count, income_count, married_count,
                                  famsize_count, age_count, female_count)



########### Prop ###############
race_prop = abind::abind(lapply(race_post, FUN = function(x) x[,"pct"]), along = 2)
rownames(race_prop) = race_post[[1]]$Race

edu_prop = abind::abind(lapply(edu_post, FUN = function(x) x[,"prop"]), along = 2)
rownames(edu_prop) = edu_post[[1]]$Edu


income_prop = abind::abind(lapply(income_post, FUN = function(x) x[,"prop"]), along = 2)
rownames(income_prop) = income_post[[1]]$income

married_prop = abind::abind(lapply(married_post, FUN = function(x) x[,"prop"]), along = 2)
rownames(married_prop) = married_post[[1]]$married

famsize_prop = abind::abind(lapply(famsize_post, FUN = function(x) x[,"prop"]), along = 2)
rownames(famsize_prop) = famsize_post[[1]]$Famsize

age_prop = abind::abind(lapply(age_post, FUN = function(x) x[,"prop"]), along = 2)
rownames(age_prop) = age_post[[1]]$age

female_prop = abind::abind(lapply(female_post, FUN = function(x) x[,"prop"]), along = 2)
rownames(female_prop) = female_post[[1]]$female


overall_demographic_prop = rbind(race_prop, edu_prop, income_prop, married_prop,
                                 famsize_prop, age_prop, female_prop)


########### Make summary table ###############
overall_demographic_tab = as.data.frame(matrix(NA, nrow = nrow(overall_demographic_count), ncol = ncol(overall_demographic_count)))
rownames(overall_demographic_tab) = rownames(overall_demographic_count)

for(i in 1:nrow(overall_demographic_count)){
  for(j in 1:ncol(overall_demographic_prop)){
    overall_demographic_tab[i,j] = sprintf("%d (%.1f%s)", overall_demographic_count[i,j], overall_demographic_prop[i,j],"%")
  }
}

xtable::xtable(overall_demographic_tab)




########### Demographic Testing ###############


fisher.test(race_count, simulate.p.value = T)
fisher.test(edu_count, simulate.p.value = T)
fisher.test(income_count, simulate.p.value = T)
fisher.test(married_count, simulate.p.value = T)
fisher.test(famsize_count, simulate.p.value = T)
fisher.test(age_count, simulate.p.value = T)

i_comb = gtools::combinations(nrow(race_count), 2)
j_comb = gtools::combinations(ncol(race_count), 2)

# for(j in 1:nrow(j_comb)){
# for(i in 1:nrow(i_comb)){
#     ft = fisher.test(race_count[i_comb[i,],j_comb[j,]])
#     if(ft$p.value < 0.05){
#       cat("i=",i_comb[i,],", j=",j_comb[j,], ", pval=",ft$p.value, "\n")
#     }
#   }
# }
# 

for(j in 1:nrow(j_comb)){
  ft = fisher.test(race_count[,j_comb[j,]])
  if(ft$p.value < 0.01){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}

i_comb = gtools::combinations(nrow(edu_count), 2)
j_comb = gtools::combinations(ncol(edu_count), 2)



# for(i in 1:nrow(i_comb)){
#   for(j in 1:nrow(j_comb)){
#     ft = fisher.test(edu_count[i_comb[i,],j_comb[j,]])
#     if(ft$p.value < 0.01){
#       cat("i=",i_comb[i,],", j=",j_comb[j,], ", pval=",ft$p.value, "\n")
#     }
#   }
# }

for(j in 1:nrow(j_comb)){
  ft = fisher.test(edu_count[,j_comb[j,]], simulate.p.value = T)
  if(ft$p.value < 0.01){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}


i_comb = gtools::combinations(nrow(income_count), 2)
j_comb = gtools::combinations(ncol(income_count), 2)

for(j in 1:nrow(j_comb)){
  ft = fisher.test(income_count[,j_comb[j,]])
  if(ft$p.value < 0.01){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}




i_comb = gtools::combinations(nrow(married_count), 2)
j_comb = gtools::combinations(ncol(married_count), 2)

for(j in 1:nrow(j_comb)){
  ft = fisher.test(married_count[,j_comb[j,]])
  if(ft$p.value < 0.01){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}


i_comb = gtools::combinations(nrow(age_count), 2)
j_comb = gtools::combinations(ncol(age_count), 2)

for(j in 1:nrow(j_comb)){
  ft = fisher.test(age_count[,j_comb[j,]])
  if(ft$p.value < 0.01){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}


i_comb = gtools::combinations(nrow(female_count), 2)
j_comb = gtools::combinations(ncol(female_count), 2)

for(j in 1:nrow(j_comb)){
  ft = fisher.test(female_count[,j_comb[j,]])
  if(ft$p.value < 0.05){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}


i_comb = gtools::combinations(nrow(famsize_count), 2)
j_comb = gtools::combinations(ncol(famsize_count), 2)

for(j in 1:nrow(j_comb)){
  ft = fisher.test(famsize_count[,j_comb[j,]])
  if(ft$p.value < 0.05){
    cat("j=",j_comb[j,], ", pval=",ft$p.value, "\n")
  }
}

