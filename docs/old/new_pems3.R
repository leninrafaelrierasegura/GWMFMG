## ---------------------------
# # Create a clipboard button on the rendered HTML page
# source(here::here("clipboard.R")); clipboard
# # Set seed for reproducibility
# set.seed(1938)
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
# truncate_to_two <- function(x) {
#   truncated <- floor(x * 100) / 100
#   sprintf("%.2f", truncated)
# }


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


## ---------------------------
# source_if_updated <- function(name) {
#   rmd_file  <- here::here(paste0(name, ".Rmd"))
#   hash_file <- here::here(paste0(name, ".hash"))
#   r_file    <- here::here(paste0(name, ".R"))
#   purl_file <- here::here("old/purl_log.txt")
# 
#   new_hash <- digest::digest(file = rmd_file, algo = "md5")
# 
#   old_hash <- if (file.exists(hash_file)) {
#     readLines(hash_file, warn = FALSE)
#   } else {
#     ""
#   }
# 
#   if (!identical(new_hash, old_hash)) {
#     capture.output(
#       knitr::purl(rmd_file, output = r_file),
#       file = purl_file
#     )
#     writeLines(new_hash, hash_file)
#   }
# 
#   source(r_file)
# }


## ---------------------------
# # capture.output(
# #   knitr::purl(here::here("functionality1.Rmd"), output = here::here("functionality1.R")),
# #   file = here::here("old/purl_log.txt")
# # )
# # source(here::here("functionality1.R"))
# source_if_updated("functionality1")


## ---------------------------
# process_model_results <- function(fit, model) {
#   fit_spde <- rspde.result(fit, "field", model, parameterization = "spde")
#   fit_matern <- rspde.result(fit, "field", model, parameterization = "matern")
#   df_for_plot_spde <- gg_df(fit_spde)
#   df_for_plot_matern <- gg_df(fit_matern)
#   param_spde <- summary(fit_spde)
#   param_matern <- summary(fit_matern)
#   param_fixed <- fit$summary.fixed[,1:6]
#   marginal.posterior.sigma_e = inla.tmarginal(
#     fun = function(x) exp(-x/2),
#     marginal = fit[["internal.marginals.hyperpar"]][["Log precision for the Gaussian observations"]])
#   quant.sigma_e <- capture.output({result_tmp <- inla.zmarginal(marginal.posterior.sigma_e)}, file = "/dev/null")
#   quant.sigma_e <- result_tmp
#   statistics.sigma_e <- unlist(quant.sigma_e)[c(1,2,3,5,7)]
#   mode.sigma_e <- inla.mmarginal(marginal.posterior.sigma_e)
#   allparams <- rbind(param_fixed, param_spde, param_matern, c(statistics.sigma_e, mode.sigma_e))
#   rownames(allparams)[nrow(allparams)] <- "sigma_e"
#   return(list(allparams = allparams, df_for_plot_spde = df_for_plot_spde, df_for_plot_matern = df_for_plot_matern))
# }


## ---------------------------
# # Load the data
# load(here("data_files/new_pems_repl1_data.RData"))
# # graph <- update_graph(graph)
# # Extract the data from the graph
# data <- graph$get_data()


## ---------------------------
# # Non-stationary parameters
# B.tau = cbind(0, 1, 0, cov, 0)
# B.kappa = cbind(0, 0, 1, 0, cov)


## ---------------------------
# load(here::here("data_files/tau_from_graphlme.RData"))
# log_tau_from_graphlme <- log(tau_from_graphlme)
# log_kappa_from_graphlme <- log(kappa_from_graphlme)


## ---------------------------
# distance = seq(from = 0, to = 10, by = 0.1)


## ---------------------------
# # Define aux data frame to compute the distance matrix
# aux <- data |> filter(repl == 1) |>
#   rename(distance_on_edge = .distance_on_edge, edge_number = .edge_number) |> # Rename the variables (because graph$compute_geodist_PtE() requires so)
#   as.data.frame() |> # Transform to a data frame (i.e., remove the metric_graph class)
#   dplyr::select(edge_number, distance_on_edge)
# 
# # Compute the distance matrix
# distmatrix <- graph$compute_geodist_PtE(PtE = aux,
#                                              normalized = TRUE,
#                                              include_vertices = FALSE)
# 
# # Compute the groups for one replicate
# GROUPS <- list()
# GROUPS2 <- list()
# for (j in 1:length(distance)) {
#   GROUPS[[j]] = list()
#   GROUPS2[[j]] = list()
#   for (i in 1:nrow(aux)) {
#     GROUPS[[j]][[i]] <- which(as.vector(distmatrix[i, ]) <= distance[j])
#     GROUPS2[[j]][[i]] <- which(as.vector(distmatrix[i, ]) > distance[j])
#   }
# }
# # Compute the groups for all replicates, based on the groups of the first replicate
# nrowY <- length(unique(data$repl))
# ncolY <- nrow(filter(data, repl == 1))
# 
# NEW_GROUPS <- list()
# for (j in 1:length(distance)) {
#   my_list <- GROUPS[[j]]
#   aux_list <- list()
#   for (i in 0:(nrowY - 1)) {
#   added_vectors <- lapply(my_list, function(vec) vec + i*ncolY)
#   aux_list <- c(aux_list, added_vectors)
#   }
#   NEW_GROUPS[[j]] <- aux_list
# }
# GROUPS <- NEW_GROUPS
# 
# NEW_GROUPS2 <- list()
# for (j in 1:length(distance)) {
#   my_list <- GROUPS2[[j]]
#   aux_list <- list()
#   for (i in 0:(nrowY - 1)) {
#     added_vectors <- lapply(my_list, function(vec) vec + i*ncolY)
#     aux_list <- c(aux_list, added_vectors)
#   }
#   NEW_GROUPS2[[j]] <- aux_list
# }
# GROUPS2 <- NEW_GROUPS2
# 
# save(GROUPS, GROUPS2, file = here::here("data_files/new_groups_for_cv_new_pems_3.RData"))


## ---------------------------
# load(here::here("data_files/new_groups_for_cv_new_pems_3.RData"))
# AUXX0 <- GROUPS2
# train_test_indices_from_GROUP2 <- list()
# for (i in seq_along(AUXX0)) {
#   AUXX1 <- AUXX0[[i]]
#   train_test_indices_from_GROUP2[[i]] <- list()
#   for (j in seq_along(AUXX1)) {
#     train_test_indices_from_GROUP2[[i]][[j]] <- list(train = list(AUXX1[[j]]), test = list(j))
#   }
# }


## ---------------------------
# # Build the model
# GRAPH_LME_statnu0.5 <-  graph_lme(y ~ mean_value,
#                                   graph = graph,
#                                   which_repl = 1:13,
#                                   model = list(type = "WhittleMatern",
#                                                # start.ltau = log_tau_from_graphlme,
#                                                # start.lkappa = log_kappa_from_graphlme,
#                                                fem = TRUE,
#                                                alpha = 1))
# summary(GRAPH_LME_statnu0.5)
# 
# POST <- posterior_crossvalidation(object = GRAPH_LME_statnu0.5, mode = "loo", true_CV = FALSE)
# MSE_GRAPH_LME_statnu0.5 <- POST$scores$rmse^2
# 
# 
# save(GRAPH_LME_statnu0.5, file = "~/Desktop/GRAPH_LME_statnu0.5.RData")


## ---------------------------
# # Build the model
# GRAPH_LME_nonstatnu0.5 <-  graph_lme(y ~ mean_value,
#                                   graph = graph,
#                                   which_repl = 1:13,
#                                   model = list(type = "WhittleMatern",
#                                                fem = TRUE,
#                                                alpha = 1,
#                                                B.tau = B.tau,
#                                                B.kappa =  B.kappa,
#                                                improve_hessian = TRUE))
# summary(GRAPH_LME_nonstatnu0.5)
# 
# POST <- posterior_crossvalidation(object = GRAPH_LME_nonstatnu0.5, mode = "loo", true_CV = FALSE)
# MSE_GRAPH_LME_nonstatnu0.5 <- POST$scores$rmse^2
# 
# 
# save(GRAPH_LME_nonstatnu0.5, file = "~/Desktop/GRAPH_LME_nonstatnu0.5.RData")


## ---------------------------
# # Build the model
# GRAPH_LME_statnu1.5 <-  graph_lme(y ~ mean_value,
#                                   graph = graph,
#                                   which_repl = 1:13,
#                                   model = list(type = "WhittleMatern",
#                                                # start.ltau = log_tau_from_graphlme,
#                                                # start.lkappa = log_kappa_from_graphlme,
#                                                fem = TRUE,
#                                                alpha = 2))
# summary(GRAPH_LME_statnu1.5)
# 
# POST <- posterior_crossvalidation(object = GRAPH_LME_statnu1.5, mode = "loo", true_CV = FALSE)
# MSE_GRAPH_LME_statnu1.5 <- POST$scores$rmse^2
# 
# 
# save(GRAPH_LME_statnu1.5, file = "~/Desktop/GRAPH_LME_statnu1.5.RData")


## ---------------------------
# # Build the model
# GRAPH_LME_nonstatnu1.5 <-  graph_lme(y ~ mean_value,
#                                   graph = graph,
#                                   which_repl = 1:13,
#                                   model = list(type = "WhittleMatern",
#                                                fem = TRUE,
#                                                alpha = 2,
#                                                B.tau = B.tau,
#                                                B.kappa =  B.kappa))
# summary(GRAPH_LME_nonstatnu1.5)
# 
# POST <- posterior_crossvalidation(object = GRAPH_LME_nonstatnu1.5, mode = "loo", true_CV = FALSE)
# MSE_GRAPH_LME_nonstatnu1.5 <- POST$scores$rmse^2
# 
# save(GRAPH_LME_nonstatnu1.5, file = "~/Desktop/GRAPH_LME_nonstatnu1.5.RData")


## ---------------------------
# # Build the model
# GRAPH_LME_statnuest <-  graph_lme(y ~ mean_value,
#                                   graph = graph,
#                                   which_repl = 1:13,
#                                   improve_hessian = TRUE,
#                                   model = list(type = "WhittleMatern",
#                                                # start.ltau = log_tau_from_graphlme,
#                                                # start.lkappa = log_kappa_from_graphlme,
#                                                fem = TRUE))
# summary(GRAPH_LME_statnuest)
# 
# POST <- posterior_crossvalidation(object = GRAPH_LME_statnuest, mode = "loo", true_CV = FALSE)
# MSE_GRAPH_LME_statnuest <- POST$scores$rmse^2
# 
# save(GRAPH_LME_statnuest, file = "~/Desktop/GRAPH_LME_statnuest.RData")


## ---------------------------
# # Build the model
# GRAPH_LME_nonstatnuest <-  graph_lme(y ~ mean_value,
#                                   graph = graph,
#                                   which_repl = 1:13,
#                                   previous_fit = GRAPH_LME_nonstatnu1.5,
#                                   model = list(type = "WhittleMatern",
#                                                fem = TRUE,
#                                                B.tau = B.tau,
#                                                B.kappa =  B.kappa))
# summary(GRAPH_LME_nonstatnuest)
# 
# POST <- posterior_crossvalidation(object = GRAPH_LME_nonstatnuest, mode = "loo", true_CV = FALSE)
# MSE_GRAPH_LME_nonstatnuest <- POST$scores$rmse^2
# 
# save(GRAPH_LME_nonstatnuest, file = "~/Desktop/GRAPH_LME_nonstatnuest.RData")


## ---------------------------
# res_exp <- graph_lme(y ~ mean_value,
#                      graph = graph,
#                      parallel = TRUE,
#                      which_repl = 1:13,
#                      model = list(type = "isoCov"))
# POST2 <- posterior_crossvalidation_loo(object = res_exp)
# MSE_ISOCOV <- POST2$scores$rmse^2


## ---------------------------
# list_MSE_GRAPH_LME <- list(
#   MSE_GRAPH_LME_statnu0.5 = MSE_GRAPH_LME_statnu0.5,
#   MSE_GRAPH_LME_nonstatnu0.5 = MSE_GRAPH_LME_nonstatnu0.5,
#   MSE_GRAPH_LME_statnu1.5 = MSE_GRAPH_LME_statnu1.5,
#   MSE_GRAPH_LME_nonstatnu1.5 = MSE_GRAPH_LME_nonstatnu1.5,
#   MSE_GRAPH_LME_statnuest = MSE_GRAPH_LME_statnuest,
#   MSE_GRAPH_LME_nonstatnuest = MSE_GRAPH_LME_nonstatnuest,
#   MSE_ISOCOV = MSE_ISOCOV
# )
# save(list_MSE_GRAPH_LME, file = here::here("data_files/list_MSE_GRAPH_LME.RData"))


## ---------------------------
# # Load the data
# load(here("data_files/new_pems_repl1_data.RData"))
# # graph <- update_graph(graph)
# # Extract the data from the graph
# data <- graph$get_data()
# data_simple <- data |> as.data.frame() |>
#   select(y, mean_value, repl)
# 
# data_simple$repl <- factor(data_simple$repl)
# 
# n <- nrow(data_simple)
# 
# pred_loocv <- numeric(n)
# 
# for(i in 1:n){
#   train_data <- data_simple[-i, ]
#   test_data  <- data_simple[i, , drop = FALSE]
# 
#   model <- lm(y ~ mean_value, data = train_data)
#   pred_loocv[i] <- predict(model, newdata = test_data)
#   print(paste("Processed observation", i, "out of", n))
# }
# 
# mse_loocv_lm <- mean((data_simple$y - pred_loocv)^2)


## ---------------------------
# # Load the data
# load(here::here("data_files/new_pems_repl1_data.RData"))
# # Extract the data from the graph
# 
# initial_data <- graph$get_data()
# data <- initial_data |> as.data.frame() |> select(y, mean_value, repl)
# 
# n <- length(data |> filter(repl == 1) |> pull(y))
# data$repl <- factor(rep(1:13, each = n))
# 
# library(FNN)
# 
# aux <- initial_data |> filter(repl == 1) |>
#   rename(distance_on_edge = .distance_on_edge,
#          edge_number = .edge_number) |>
#   as.data.frame() |>
#   dplyr::select(edge_number,
#                 distance_on_edge)
# 
# D <- graph$compute_geodist_PtE(
#   PtE = aux,
#   normalized = TRUE,
#   include_vertices = FALSE)
# 
# # -----------------------------------------------------------
# # 0. Setup
# # -----------------------------------------------------------
# n_loc  <- 314
# n_repl <- 13
# n_obs  <- nrow(data)  # 4082
# 
# # Location index for each observation in data
# # (assumes data is ordered: all 314 locations for repl 1, then repl 2, etc.)
# loc_idx <- rep(1:n_loc, times = n_repl)
# 
# # -----------------------------------------------------------
# # 1. Build normalized distance matrices
# # -----------------------------------------------------------
# 
# # Spatial: expand 314x314 -> 4082x4082 using location indices
# D_space_full <- D[loc_idx, loc_idx]
# D_space_norm <- D_space_full / max(D_space_full)
# 
# # Covariate: pairwise distances on mean_value across all 4082 observations
# D_cov_full <- as.matrix(dist(scale(data$mean_value)))
# D_cov_norm <- D_cov_full / max(D_cov_full)
# 
# # -----------------------------------------------------------
# # 2. Combined distance (precomputed, outside all loops)
# # -----------------------------------------------------------
# alpha <- 0  # 0 = pure spatial, 1 = pure covariate
# D_combined <- alpha * D_cov_norm + (1 - alpha) * D_space_norm
# 
# # -----------------------------------------------------------
# # 3. LOO cross-validation over k
# # -----------------------------------------------------------
# k_values <- 1:30
# 
# loo_mse <- sapply(k_values, function(k) {
#   pred <- numeric(n_obs)
# 
#   for (i in 1:n_obs) {
#     neighbors <- order(D_combined[i, -i])[1:k]
#     pred[i]   <- mean(data$y[-i][neighbors])
#   }
# 
#   mean((data$y - pred)^2)
# })
# 
# # -----------------------------------------------------------
# # 4. Results
# # -----------------------------------------------------------
# 
# best_k    <- k_values[which.min(loo_mse)]
# best_KNN_mse <- loo_mse[which.min(loo_mse)]


## ---------------------------
# save(mse_loocv_lm, best_KNN_mse, file = here::here("data_files/new_simple_linear_regression_results.RData"))


## ----eval = TRUE------------
load(here::here("data_files/new_simple_linear_regression_results.RData"))
load(here::here("data_files/list_MSE_GRAPH_LME.RData"))

# format as data frame
list_MSE_GRAPH_LME_df <- data.frame(Model = names(list_MSE_GRAPH_LME), MSE = unlist(list_MSE_GRAPH_LME), row.names = NULL) |> 
  bind_rows(data.frame(Model = "Simple Linear Regression", MSE = mse_loocv_lm)) |>
  bind_rows(data.frame(Model = "kNN Regression", MSE = best_KNN_mse)) 

list_MSE_GRAPH_LME_df


## ----eval = TRUE------------
grateful::cite_packages(output = "paragraph", out.dir = ".")

