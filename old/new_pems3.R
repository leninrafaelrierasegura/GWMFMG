EXPNUM <- 1


# # Set seed for reproducibility
set.seed(1938)
# # Set global options for all code chunks
# knitr::opts_chunk$set(
#   # Disable messages printed by R code chunks
#   message = FALSE,
#   # Disable warnings printed by R code chunks
#   warning = FALSE,
#   # Show R code within code chunks in output
#   echo = TRUE,
#   # Include both R code and its results in output
#   include = TRUE,
#   # Evaluate R code chunks
#   eval = FALSE,
#   # Enable caching of R code chunks for faster rendering
#   cache = FALSE,
#   # Align figures in the center of the output
#   fig.align = "center",
#   # Enable retina display for high-resolution figures
#   retina = 2,
#   # Show errors in the output instead of stopping rendering
#   error = TRUE,
#   # Do not collapse code and output into a single block
#   collapse = FALSE
# )
# # Start the figure counter
# fig_count <- 0
# # Define the captioner function
# captioner <- function(caption) {
#   fig_count <<- fig_count + 1
#   paste0("Figure ", fig_count, ": ", caption)
# }
# 
# 
# # Define the function to truncate a number to two decimal places
truncate_to_two <- function(x) {
  truncated <- floor(x * 100) / 100
  sprintf("%.2f", truncated)
}


## ----eval = TRUE------------

# inla.upgrade(testing = TRUE)
# remotes::install_github("inlabru-org/inlabru", ref = "devel")
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(INLA)
library(inlabru)
library(rSPDE)
library(MetricGraph)

library(dplyr)
library(plotly)
library(scales)
library(patchwork)

library(ggplot2)
library(cowplot)
library(ggpubr) #annotate_figure()
library(grid) #textGrob()
library(ggmap)

library(viridis)
library(OpenStreetMap)


library(tidyr)
library(sf)

library(here)
library(rmarkdown)
library(grateful) # Cite all loaded packages


library(slackr)
source("keys.R")
slackr_setup(token = token) # token comes from keys.R


## ---------------------------
# Load the data
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_pems_repl1_data.RData"))
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/tau_from_graphlme.RData"))

load(here::here("data_files/new_pems_repl1_data.RData"))
load(here::here("data_files/tau_from_graphlme.RData"))

# ---------------------------
# Non-stationary parameters
B.tau = cbind(0, 1, 0, cov, 0)
B.kappa = cbind(0, 0, 1, 0, cov)



log_tau_from_graphlme <- log(tau_from_graphlme)
log_kappa_from_graphlme <- log(kappa_from_graphlme)


#####################################
#############nu=0.5##################
#####################################

# ---------------------------
# Build the model
GRAPH_LME_statnu0.5 <-  graph_lme(y ~ mean_value,
                                  graph = graph,
                                  which_repl = 1:13,
                                  parallel = TRUE,
                                  #improve_hessian = TRUE,
                                  # model_options = list(start_tau = tau_from_graphlme,
                                  #                      start_kappa = kappa_from_graphlme),
                                  model = list(type = "WhittleMatern",
                                               fem = TRUE,
                                               alpha = 1)
                                  )
summary(GRAPH_LME_statnu0.5)

POST <- posterior_crossvalidation(object = GRAPH_LME_statnu0.5, mode = "loo", true_CV = FALSE)
MSE_GRAPH_LME_statnu0.5 <- POST$scores$rmse^2

slackr_msg(text = paste0("MSE_GRAPH_LME_statnu0.5 = ", 
                         MSE_GRAPH_LME_statnu0.5), 
           channel = "#research")

slackr_msg(
  text = paste(
    capture.output(print(GRAPH_LME_statnu0.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

sigma_e_statnu0.5 <- GRAPH_LME_statnu0.5$coeff$measurement_error |> as.numeric()
tau_statnu0.5 <- GRAPH_LME_statnu0.5$coeff$random_effects["tau"] |> as.numeric()
kappa_statnu0.5 <- GRAPH_LME_statnu0.5$coeff$random_effects["kappa"] |> as.numeric()

save(sigma_e_statnu0.5, tau_statnu0.5, kappa_statnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_statnu0.5_parameters.RData"))
save(GRAPH_LME_statnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_statnu0.5.RData"))


# ---------------------------
# Build the model
GRAPH_LME_nonstatnu0.5 <-  graph_lme(y ~ mean_value,
                                  graph = graph,
                                  which_repl = 1:13,
                                  parallel = TRUE,
                                  #improve_hessian = TRUE,
                                  model = list(type = "WhittleMatern",
                                               fem = TRUE,
                                               alpha = 1,
                                               B.tau = B.tau,
                                               B.kappa =  B.kappa))
summary(GRAPH_LME_nonstatnu0.5)

POST <- posterior_crossvalidation(object = GRAPH_LME_nonstatnu0.5, mode = "loo", true_CV = FALSE)
MSE_GRAPH_LME_nonstatnu0.5 <- POST$scores$rmse^2


slackr_msg(text = paste0("MSE_GRAPH_LME_nonstatnu0.5 = ", 
                         MSE_GRAPH_LME_nonstatnu0.5), 
           channel = "#research")
slackr_msg(
  text = paste(
    capture.output(print(GRAPH_LME_nonstatnu0.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

sigma_e_nonstatnu0.5 <- GRAPH_LME_nonstatnu0.5$coeff$measurement_error |> as.numeric()
theta_nonstatnu0.5 <- c(0,GRAPH_LME_nonstatnu0.5$coeff$random_effects[2:5]) |> as.vector()

save(sigma_e_nonstatnu0.5, theta_nonstatnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_nonstatnu0.5_parameters.RData"))
save(GRAPH_LME_nonstatnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_nonstatnu0.5.RData"))


#####################################
#############nu=1.5##################
#####################################

# ---------------------------
# Build the model
GRAPH_LME_statnu1.5 <-  graph_lme(y ~ mean_value,
                                  graph = graph,
                                  which_repl = 1:13,
                                  parallel = TRUE,
                                  #improve_hessian = TRUE,
                                  # model_options = list(start_tau = tau_from_graphlme,
                                  #                      start_kappa = kappa_from_graphlme),
                                  model = list(type = "WhittleMatern",
                                               fem = TRUE,
                                               alpha = 2))
summary(GRAPH_LME_statnu1.5)

POST <- posterior_crossvalidation(object = GRAPH_LME_statnu1.5, mode = "loo", true_CV = FALSE)
MSE_GRAPH_LME_statnu1.5 <- POST$scores$rmse^2

slackr_msg(text = paste0("MSE_GRAPH_LME_statnu1.5 = ", 
                         MSE_GRAPH_LME_statnu1.5), 
           channel = "#research")

slackr_msg(
  text = paste(
    capture.output(print(GRAPH_LME_statnu1.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

sigma_e_statnu1.5 <- GRAPH_LME_statnu1.5$coeff$measurement_error |> as.numeric()
tau_statnu1.5 <- GRAPH_LME_statnu1.5$coeff$random_effects["tau"] |> as.numeric()
kappa_statnu1.5 <- GRAPH_LME_statnu1.5$coeff$random_effects["kappa"] |> as.numeric()

save(sigma_e_statnu1.5, tau_statnu1.5, kappa_statnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_statnu1.5_parameters.RData"))
save(GRAPH_LME_statnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_statnu1.5.RData"))


# ---------------------------
# Build the model
GRAPH_LME_nonstatnu1.5 <-  graph_lme(y ~ mean_value,
                                  graph = graph,
                                  which_repl = 1:13,
                                  parallel = TRUE,
                                  #improve_hessian = TRUE,
                                  model = list(type = "WhittleMatern",
                                               fem = TRUE,
                                               alpha = 2,
                                               B.tau = B.tau,
                                               B.kappa =  B.kappa))
summary(GRAPH_LME_nonstatnu1.5)

POST <- posterior_crossvalidation(object = GRAPH_LME_nonstatnu1.5, mode = "loo", true_CV = FALSE)
MSE_GRAPH_LME_nonstatnu1.5 <- POST$scores$rmse^2


slackr_msg(text = paste0("MSE_GRAPH_LME_nonstatnu1.5 = ", 
                         MSE_GRAPH_LME_nonstatnu1.5), 
           channel = "#research")

slackr_msg(
  text = paste(
    capture.output(print(GRAPH_LME_nonstatnu1.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

sigma_e_nonstatnu1.5 <- GRAPH_LME_nonstatnu1.5$coeff$measurement_error |> as.numeric()
theta_nonstatnu1.5 <- c(0,GRAPH_LME_nonstatnu1.5$coeff$random_effects[2:5]) |> as.vector()

save(sigma_e_nonstatnu1.5, theta_nonstatnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_nonstatnu1.5_parameters.RData"))
save(GRAPH_LME_nonstatnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_nonstatnu1.5.RData"))


#####################################
#############nu=est##################
#####################################

# ---------------------------
# Build the model
GRAPH_LME_statnuest <-  graph_lme(y ~ mean_value,
                                  graph = graph,
                                  which_repl = 1:13,
                                  parallel = TRUE,
                                  improve_hessian = TRUE,
                                  # model_options = list(start_tau = tau_from_graphlme,
                                  #                      start_kappa = kappa_from_graphlme),
                                  model = list(type = "WhittleMatern",
                                               fem = TRUE))
summary(GRAPH_LME_statnuest)

POST <- posterior_crossvalidation(object = GRAPH_LME_statnuest, mode = "loo", true_CV = FALSE)
MSE_GRAPH_LME_statnuest <- POST$scores$rmse^2


slackr_msg(text = paste0("MSE_GRAPH_LME_statnuest = ", 
                         MSE_GRAPH_LME_statnuest), 
           channel = "#research")

slackr_msg(
  text = paste(
    capture.output(print(GRAPH_LME_statnuest)),
    collapse = "\n"
  ),
  channel = "#research"
)

sigma_e_statnuest <- GRAPH_LME_statnuest$coeff$measurement_error |> as.numeric()
alpha_statnuest <- GRAPH_LME_statnuest$coeff$random_effects["alpha"] |> as.numeric()
tau_statnuest <- GRAPH_LME_statnuest$coeff$random_effects["tau"] |> as.numeric()
kappa_statnuest <- GRAPH_LME_statnuest$coeff$random_effects["kappa"] |> as.numeric()

save(sigma_e_statnuest, tau_statnuest, kappa_statnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_statnuest_parameters.RData"))
save(GRAPH_LME_statnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_statnuest.RData"))


# ---------------------------
# Build the model
GRAPH_LME_nonstatnuest <-  graph_lme(y ~ mean_value,
                                  graph = graph,
                                  which_repl = 1:13,
                                  parallel = TRUE,
                                  #improve_hessian = TRUE,
                                  previous_fit = GRAPH_LME_nonstatnu1.5,
                                  model = list(type = "WhittleMatern",
                                               fem = TRUE,
                                               B.tau = B.tau,
                                               B.kappa =  B.kappa))
summary(GRAPH_LME_nonstatnuest)

POST <- posterior_crossvalidation(object = GRAPH_LME_nonstatnuest, mode = "loo", true_CV = FALSE)
MSE_GRAPH_LME_nonstatnuest <- POST$scores$rmse^2


slackr_msg(text = paste0("MSE_GRAPH_LME_nonstatnuest = ", 
                         MSE_GRAPH_LME_nonstatnuest), 
           channel = "#research")

slackr_msg(
  text = paste(
    capture.output(print(GRAPH_LME_nonstatnuest)),
    collapse = "\n"
  ),
  channel = "#research"
)

sigma_e_nonstatnuest <- GRAPH_LME_nonstatnuest$coeff$measurement_error |> as.numeric()
alpha_nonstatnuest <- GRAPH_LME_nonstatnuest$coeff$random_effects["alpha"] |> as.numeric()
theta_nonstatnuest <- c(0,GRAPH_LME_nonstatnuest$coeff$random_effects[2:5]) |> as.vector()

save(sigma_e_nonstatnuest, theta_nonstatnuest, alpha_nonstatnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_nonstatnuest_parameters.RData"))
save(GRAPH_LME_nonstatnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/GRAPH_LME_nonstatnuest.RData"))


#####################################
#############isoCov##################
#####################################

# ---------------------------
res_exp <- graph_lme(y ~ mean_value,
                     graph = graph,
                     parallel = TRUE,
                     which_repl = 1:13,
                     model = list(type = "isoCov"))
POST2 <- posterior_crossvalidation_loo(object = res_exp)
MSE_ISOCOV <- POST2$scores$rmse^2

slackr_msg(text = paste0("MSE_ISOCOV = ", 
                         MSE_ISOCOV), 
           channel = "#research")

slackr_msg(
  text = paste(
    capture.output(print(res_exp)),
    collapse = "\n"
  ),
  channel = "#research"
)

# ---------------------------
list_MSE_GRAPH_LME <- list(
  MSE_GRAPH_LME_statnu0.5 = MSE_GRAPH_LME_statnu0.5,
  MSE_GRAPH_LME_nonstatnu0.5 = MSE_GRAPH_LME_nonstatnu0.5,
  MSE_GRAPH_LME_statnu1.5 = MSE_GRAPH_LME_statnu1.5,
  MSE_GRAPH_LME_nonstatnu1.5 = MSE_GRAPH_LME_nonstatnu1.5,
  MSE_GRAPH_LME_statnuest = MSE_GRAPH_LME_statnuest,
  MSE_GRAPH_LME_nonstatnuest = MSE_GRAPH_LME_nonstatnuest,
  MSE_ISOCOV = MSE_ISOCOV
)
save(list_MSE_GRAPH_LME, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/list_MSE_GRAPH_LME.RData"))


# ---------------------------
# Load the data
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_pems_repl1_data.RData"))
load(here::here("data_files/new_pems_repl1_data.RData"))
# graph <- update_graph(graph)
# Extract the data from the graph
data <- graph$get_data()
data_simple <- data |> as.data.frame() |>
  select(y, mean_value, repl)

data_simple$repl <- factor(data_simple$repl)

n <- nrow(data_simple)

pred_loocv <- numeric(n)

for(i in 1:n){
  train_data <- data_simple[-i, ]
  test_data  <- data_simple[i, , drop = FALSE]

  model <- lm(y ~ mean_value, data = train_data)
  pred_loocv[i] <- predict(model, newdata = test_data)
  #print(paste("Processed observation", i, "out of", n))
}

mse_loocv_lm <- mean((data_simple$y - pred_loocv)^2)


# ---------------------------
# Load the data
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_pems_repl1_data.RData"))
load(here::here("data_files/new_pems_repl1_data.RData"))
# Extract the data from the graph

initial_data <- graph$get_data()
data <- initial_data |> as.data.frame() |> select(y, mean_value, repl)

n <- length(data |> filter(repl == 1) |> pull(y))
data$repl <- factor(rep(1:13, each = n))

library(FNN)

aux <- initial_data |> filter(repl == 1) |>
  rename(distance_on_edge = .distance_on_edge,
         edge_number = .edge_number) |>
  as.data.frame() |>
  dplyr::select(edge_number,
                distance_on_edge)

D <- graph$compute_geodist_PtE(
  PtE = aux,
  normalized = TRUE,
  include_vertices = FALSE)

# -----------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------
n_loc  <- 314
n_repl <- 13
n_obs  <- nrow(data)  # 4082

# Location index for each observation in data
# (assumes data is ordered: all 314 locations for repl 1, then repl 2, etc.)
loc_idx <- rep(1:n_loc, times = n_repl)

# -----------------------------------------------------------
# 1. Build normalized distance matrices
# -----------------------------------------------------------

# Spatial: expand 314x314 -> 4082x4082 using location indices
D_space_full <- D[loc_idx, loc_idx]
D_space_norm <- D_space_full / max(D_space_full)

# Covariate: pairwise distances on mean_value across all 4082 observations
D_cov_full <- as.matrix(dist(scale(data$mean_value)))
D_cov_norm <- D_cov_full / max(D_cov_full)

# -----------------------------------------------------------
# 2. Combined distance (precomputed, outside all loops)
# -----------------------------------------------------------
alpha <- 0  # 0 = pure spatial, 1 = pure covariate
D_combined <- alpha * D_cov_norm + (1 - alpha) * D_space_norm

# -----------------------------------------------------------
# 3. LOO cross-validation over k
# -----------------------------------------------------------
k_values <- 1:30

loo_mse <- sapply(k_values, function(k) {
  pred <- numeric(n_obs)

  for (i in 1:n_obs) {
    neighbors <- order(D_combined[i, -i])[1:k]
    pred[i]   <- mean(data$y[-i][neighbors])
  }

  mean((data$y - pred)^2)
})

# -----------------------------------------------------------
# 4. Results
# -----------------------------------------------------------

best_k    <- k_values[which.min(loo_mse)]
best_KNN_mse <- loo_mse[which.min(loo_mse)]


# ---------------------------
save(mse_loocv_lm, best_KNN_mse, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_simple_linear_regression_results.RData"))


## ----eval = TRUE------------
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_simple_linear_regression_results.RData"))
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/list_MSE_GRAPH_LME.RData"))

# format as data frame
list_MSE_GRAPH_LME_df <- data.frame(Model = names(list_MSE_GRAPH_LME), MSE = unlist(list_MSE_GRAPH_LME), row.names = NULL) |> 
  bind_rows(data.frame(Model = "Simple Linear Regression", MSE = mse_loocv_lm)) |>
  bind_rows(data.frame(Model = "kNN Regression", MSE = best_KNN_mse)) 

list_MSE_GRAPH_LME_df
