# Display the neighborhood correlation for different set of hyperparameter
library(BayesGPfit)
setwd("/Volumes/easystore/Dropbox/U-Mich/Research/WIR_ABCD")
ref_input = readRDS("./Simulation/Inputs/input_2_0_thres08.rds")
grids = ref_input$loc
V = nrow(grids)
head(grids)
grids_point_1 = grids[1, ]
grids_point_2 = grids[2, ]

# Compute the correlation between the two points
a = 0.01
b = 80
corr_se = function (grid1, grid2, a, b) {
    dist = sqrt(sum((grid1 - grid2)^2))
    norm1 = sqrt(sum(grid1^2))
    norm2 = sqrt(sum(grid2^2))
    return (exp(-a * (norm1^2 + norm2^2) - b*dist^2))
}

corr_se(grids_point_1, grids_point_2, 0.01, 80)
corr_se(grids_point_1, grids_point_2, 0.01, 120)
corr_se(grids_point_1, grids_point_2, 0.01, 200)
corr_se(grids_point_1, grids_point_2, 0.01, 300)
corr_se(grids_point_1, grids_point_2, 0.01, 1250)

