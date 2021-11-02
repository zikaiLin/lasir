# 3lasir

LAtent Subgroup Image-on-scalar Regression (LASIR)



# Installation

```R
library(devtools)
install_github("https://github.com/zikaiLin/lasir")
```



# Example

To demonstrate how LASIR works, here's an toy example containing 500 individuals from 5 study sites.

```R
N = 500
n_sites = 5
```



Here we consider 2D images, with 50 voxels at each dimensions:

```R

# Simulate a set of 2-d images
grids = BayesGPfit::GP.generate.grids(d = 2, num_grids = 50)
V = 50*50 # Total number of voxels
```



We consider two main effect spatial varying coefficient structure (sparse, but with different sign): 

```R
main_effect_coef_1 = matrix(0, 50, 50)
main_effect_coef_1[15:35,15:35] = 2
main_effect_coef_1 = fields::image.smooth(main_effect_coef_1)$z


main_effect_coef_2 = matrix(0, 50, 50)
main_effect_coef_2[20:30,20:30] = -2
main_effect_coef_2 = fields::image.smooth(main_effect_coef_2)$z

main_effect_coef_intercept = matrix(0.2, 50, 50)


```

```R
# Visualize SVCs of main effect
library(ggplot2)

# Dummy dat
data <- expand.grid(X = 1:50, Y = 1:50)
data$Z1 <- c(main_effect_coef_1)
data$Z2 <- c(main_effect_coef_2)

# Heatmap 
ggplot(data, aes(X, Y, fill= Z1)) + 
  geom_tile()
ggplot(data, aes(X, Y, fill= Z2)) + 
  geom_tile()


```

![image-20211102000636728](/Users/ziklin/Library/Application Support/typora-user-images/image-20211102000636728.png)



Generate covariates:
$$
\begin{align}
X_i&\sim N(2,1)\\
Z_i&\sim N(1,1)\\
U_i &\sim Cat(3, 1/3)\\
\gamma_s(v) &\sim N(0,1); \quad \forall v\\

\end{align}
$$

```R
# Simulate fixed effect and sites effect
X = rnorm(N, 2, sd = 1) 
Z = rnorm(N, mean = 1, sd = 1)
sites = sample(1:n_sites, size = N, replace = T)

fixed_coef = rnorm(V, mean = 0, sd = 1)
sites_coef = matrix(0, n_sites, V)
for(s in 1:n_sites){
  sites_coef[s,] = rnorm(V, mean = 0, sd = 1)  
}
```

Finally, simulate the outcome based on the following true model:
$$
y_i(v) = \beta_{}
$$


```R
# Simulate group indicator for each individual
nK = 2 # 2-subgroups
grp_ind = sample(1:2, size = N, replace = T)

grp1_idx = (grp_ind==1)
grp2_idx = (grp_ind==2)

# Simulate outcome image
Y = matrix(NA, nrow = N, ncol = V)
Y[grp1_idx,] = cbind(1, X[grp1_idx]) %*% rbind(c(main_effect_coef_intercept), c(main_effect_coef_1)) + sites_coef[sites[grp1_idx],] + Z[grp1_idx]%*% t(fixed_coef)
Y[grp2_idx,] = cbind(1, X[grp2_idx]) %*% rbind(c(main_effect_coef_intercept), c(main_effect_coef_2)) + sites_coef[sites[grp2_idx],] + Z[grp2_idx]%*% t(fixed_coef)
```



```R

# Eigen-transformation applied on 
eigen_trans = eigen.transform(Y, grids,poly_degree = 6)

# SEM fit
lasir_fit = LASIR.fit(eigen_trans$Y_star, as.matrix(X), as.matrix(Z), sites, nK = 2, max.iter = 20)

# Model Inference
lasir_inference = LASIR.inference(lasir_fit, eigen_trans)

# Dummy dat
data_fit <- expand.grid(X = 1:50, Y = 1:50)
data_fit$Z1 <- c(lasir_inference$main_effect_svcs[2,,1])
data_fit$Z2 <- c(lasir_inference$main_effect_svcs[2,,2])

# Heatmap 
ggplot(data_fit, aes(X, Y, fill= Z1)) + 
  geom_tile()
ggplot(data_fit, aes(X, Y, fill= Z2)) + 
  geom_tile()
# Plot
aricode::NMI(lasir_fit$cluster, grp_ind)


```

