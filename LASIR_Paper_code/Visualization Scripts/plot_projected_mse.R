load("./EM_RDA/section4-cross-validation/mse_lasir_cross_validation_shuffle.RData")
load("./EM_RDA/section4-cross-validation/mse_lasir_cross_validation.RData")
load("./EM_RDA/section4-cross-validation/mse_svcm_cross_validation.RData")
# load("./EM_RDA/section4-cross-validation/mse_socioClust_cross_validation.RData")
# 

# projected_MSE = data.frame(mse = c(mse_y_lasir, mse_y_lasir_shuffle, mse_y_svcm, mse_y_socioClust),
           # Approach = factor(rep(c("Within Subgroup", "W/O Subgroup", "Randomized Subgroup", "SocioClust"), each = 50)))

projected_MSE = data.frame(mse = c(mse_y_lasir, mse_y_lasir_shuffle, mse_y_svcm),
                           Approach = factor(rep(c("Within Subgroup", "W/O Subgroup", "Randomized Subgroup"), each = 50)))

library(ggthemes)
library(ggplot2)
g <- ggplot(projected_MSE, aes(Approach, mse))+ geom_boxplot() + theme_bw()+
  theme(axis.text.x = element_text(angle=0, vjust=0.6, size = 24, family = "Times New Roman",face = "bold"),
        axis.text.y = element_text(angle=0, vjust=0.6, size = 18, family = "Times New Roman",face = "bold"),
        axis.title =  element_text(angle=0, vjust = 0.1, size = 24, face = "bold", family = "Times New Roman", 
                                   margin = margin(t = 30, r = 30, b = 30, l = 30)),
        # legend.text = element_text(angle=0, vjust=0.6, size = 18, family = "Times New Roman"),
        # legend.title = element_text(angle=0, vjust=0.6, size = 24, face = "bold", family = "Times New Roman"),
        title =  element_text(angle=0,size = 24, family = "Times New Roman")) + 
  labs(
    # title="Projected Mean Square Error Given by each Approach", 
       x="",
       y="Projectd Mean Square Error")
g


ggsave("./EM_RDA/section4-cross-validation/mse_cv.eps",
       plot = g, device = "eps")

ggsave("./EM_RDA/section4-cross-validation/mse_cv.png",
       plot = g, device = "png")
