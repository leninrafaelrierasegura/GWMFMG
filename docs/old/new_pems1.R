EXPNUM <- 1

## ---------------------------
# # Create a clipboard button on the rendered HTML page
# source(here::here("clipboard.R")); clipboard
# Set seed for reproducibility
set.seed(1982)
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
# # Define the function to truncate a number to two decimal places
truncate_to_two <- function(x) {
  truncated <- floor(x * 100) / 100
  sprintf("%.2f", truncated)
}


## ----eval = TRUE------------
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



# ---------------------------
gets_summary_parameters <- function(fit, model) {
  param_spde <- summary(rspde.result(fit, "field", model, parameterization = "spde"))
  param_matern <- summary(rspde.result(fit, "field", model, parameterization = "matern"))
  param_fixed <- fit$summary.fixed[,1:6]
  marginal.posterior.sigma_e = inla.tmarginal(
    fun = function(x) exp(-x/2),
    marginal = fit[["internal.marginals.hyperpar"]][["Log precision for the Gaussian observations"]])
  quant.sigma_e <- capture.output({result_tmp <- inla.zmarginal(marginal.posterior.sigma_e)}, file = "/dev/null")
  quant.sigma_e <- result_tmp
  statistics.sigma_e <- unlist(quant.sigma_e)[c(1,2,3,5,7)]
  mode.sigma_e <- inla.mmarginal(marginal.posterior.sigma_e)
  allparams <- rbind(param_fixed, param_spde, param_matern, c(statistics.sigma_e, mode.sigma_e))
  rownames(allparams)[nrow(allparams)] <- "sigma_e"
  return(allparams)
}


# ---------------------------
Lines <- read_sf(here("data_files/data.pems/lines.shp"))
lines <- as_Spatial(Lines)

EtV <- read.csv(here("data_files/data.pems/E.csv"), header = T, row.names = NULL)
PtE <- read.csv(here("data_files/data.pems/PtE.csv"), header = T, row.names = NULL)
PtE[,1] <- PtE[,1] + 1
Y <- read.csv(here("data_files/data.pems/Y.csv"), header = T, row.names = NULL)
Y <- as.matrix(Y[,-1])
edge_length_m <- EtV[,4]
PtE[,2] = PtE[,2]/edge_length_m[PtE[,1]]


# ---------------------------
all_same <- function(col) {
  length(unique(col[1:13])) == 1
}
cols_to_keep <- apply(Y, 2, function(col) !all_same(col))

Y <- Y[, cols_to_keep]
Y_raw <- Y
PtE <- PtE[cols_to_keep,]
PtE_raw <- PtE


# ---------------------------
sampled_numbers <- 1:13
not_sampled_numbers <- 14:26
Y_for_summary <- Y[sampled_numbers,]
Y_logstd <- apply(Y_for_summary, 2, sd) |> as.vector() |> log()
Y_mean <- apply(Y_for_summary, 2, mean) |> as.vector()
Y <- Y[not_sampled_numbers,]
save(Y_mean,Y_raw, PtE_raw, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_Y_mean.RData"))
save(Y_logstd, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/Y_logstd.RData"))


# ---------------------------
# Build the graph
graph <- metric_graph$new(edges = pems_repl$edges, longlat = TRUE)
# Add the observations
graph$add_observations(data = data.frame(y = Y_logstd,
                                         edge_number = PtE[,1],
                                         distance_on_edge = PtE[,2]),
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE,
                       clear_obs = TRUE)
# Build the mesh
graph$build_mesh(h = 0.05)


# ---------------------------
# Build the model
rspde_fit_stat_logstd <-  graph_lme(y ~ 1,
                                    graph = graph,
                                    #improve_hessian = TRUE,
                                    model = list(type = "WhittleMatern",
                                                 fem = FALSE,
                                                 alpha = 2),
                                    #model_options = list(fix_sigma_e = 10^-2),
                                    parallel = TRUE)

save(rspde_fit_stat_logstd, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_rspde_fit_stat_and_model_logstd.RData"))


## ----eval = TRUE------------
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_rspde_fit_stat_and_model_logstd.RData"))
# Summarize the results
summary(rspde_fit_stat_logstd)

# ---------------------------
# Prediction locations
data_prd_list_mesh <- graph$mesh$VtE
# Compute the kriging predictor
y_pred <- predict(rspde_fit_stat_logstd, data.frame(edge_number = data_prd_list_mesh[,1], distance_on_edge = data_prd_list_mesh[,2]), normalized = TRUE)
cov_almost <- y_pred$mean
# Standardize the kriging predictor of the log standard deviation
cov <- (cov_almost - mean(cov_almost))/sd(cov_almost)

# this is just to check

SSS <- rename(PtE, .edge_number = E, .distance_on_edge = length)
auxxx <- predict(rspde_fit_stat_logstd, data.frame(edge_number = SSS$.edge_number, distance_on_edge = SSS$.distance_on_edge), normalized = TRUE)
tocompatewithlogstd <- auxxx$mean
save(tocompatewithlogstd, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/tocompatewithlogstd.RData"))


## ----eval =TRUE, fig.width = 24, fig.height = 4, fig.cap = captioner("Predicted log std and observed log std")----
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/Y_logstd.RData"))
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/tocompatewithlogstd.RData"))
plot(Y_logstd, type = "l", col = "darkblue")
lines(tocompatewithlogstd, col = "darkred")



# ---------------------------
graph$add_observations(data = data.frame(y = Y_mean,
                                         edge_number = PtE[,1],
                                         distance_on_edge = PtE[,2]),
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE,
                       clear_obs = TRUE)


# ---------------------------
# Build the model
rspde_fit_stat_mean <-  graph_lme(y ~ 1,
                                    graph = graph,
                                    #improve_hessian = TRUE,
                                    model = list(type = "WhittleMatern",
                                                 fem = FALSE,
                                                 alpha = 2),
                                    #model_options = list(fix_sigma_e = 10^-2),
                                    parallel = TRUE)
save(rspde_fit_stat_mean, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_rspde_fit_stat_and_model_mean.RData"))


## ----eval = TRUE------------
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_rspde_fit_stat_and_model_mean.RData"))
# Summarize the results
summary(rspde_fit_stat_mean)


# ---------------------------
mean_on_mesh_pred <- predict(rspde_fit_stat_mean, data.frame(edge_number = data_prd_list_mesh[,1], distance_on_edge = data_prd_list_mesh[,2]), normalized = TRUE)
cov_for_mean_to_plot <- mean_on_mesh_pred$mean

mean_on_loc_pred <- predict(rspde_fit_stat_mean, data.frame(edge_number = SSS$.edge_number, distance_on_edge = SSS$.distance_on_edge), normalized = TRUE)
cov_for_mean <- mean_on_loc_pred$mean
save(cov_for_mean_to_plot, cov_for_mean, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_cov_for_mean_to_plot_and_cov_for_mean.RData"))


## ----eval =TRUE, fig.width = 24, fig.height = 4, fig.cap = captioner("Predicted mean and observed mean.")----
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_Y_mean.RData"))
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_cov_for_mean_to_plot_and_cov_for_mean.RData"))
plot(Y_mean, type = "l", col = "darkblue")
lines(cov_for_mean, col = "darkred")



# ---------------------------
df_rep <- lapply(1:nrow(Y), function(i){data.frame(y = Y[i,],
                                                   mean_value = cov_for_mean,
                                                   edge_number = PtE[,1],
                                                   distance_on_edge = PtE[,2],
                                                   repl = i)})
df_rep <- do.call(rbind, df_rep)

graph$add_observations(data = df_rep,
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE,
                       clear_obs = TRUE,
                       group = "repl")

save(graph, cov, cov_for_mean_to_plot, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_pems_repl1_data.RData"))







Y_mu <- apply(Y_raw[1:13,], 2, mean) |> as.vector()


df_isocov <- data.frame(y = Y_mu, 
                        edge_number = PtE_raw[,1], 
                        distance_on_edge = PtE_raw[,2])

graph$add_observations(data = df_isocov,
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE, 
                       clear_obs = TRUE)

res_exp <- graph_lme(y ~ 1, 
                     graph = graph, 
                     model = list(type = "isoCov"), 
                     parallel = TRUE)
summary(res_exp)
u_est_exp_mean <- predict(res_exp, df_isocov[,c("edge_number", "distance_on_edge")], normalized = TRUE)$mean

plot(Y_mu, type = "l", col = "darkblue")
lines(u_est_exp_mean, col = "darkred")


Y2part <- Y_raw[14:26,]
DF_ISOCOV <- lapply(1:nrow(Y2part), function(i){data.frame(y = Y2part[i,],
                                                           mean_value = u_est_exp_mean,
                                                           edge_number = PtE_raw[,1],
                                                           distance_on_edge = PtE_raw[,2],
                                                           repl = i)})
DF_ISOCOV <- do.call(rbind, DF_ISOCOV)

graph$add_observations(data = DF_ISOCOV, 
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE, 
                       clear_obs = TRUE, 
                       group = "repl")

RES_EXP <- graph_lme(y ~ mean_value, 
                     graph = graph, 
                     which_repl = 1:13, 
                     model = list(type = "isoCov"), 
                     parallel = TRUE)
summary(RES_EXP)

tau_from_graphlme <- RES_EXP$coeff$random_effects[1] |> as.numeric()
kappa_from_graphlme <- RES_EXP$coeff$random_effects[2] |> as.numeric()
save(tau_from_graphlme, kappa_from_graphlme, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/tau_from_graphlme.RData"))




