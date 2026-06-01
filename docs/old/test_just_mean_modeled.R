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
graph$add_observations(data = data.frame(y = Y_mean, 
                                         edge_number = PtE[,1], 
                                         distance_on_edge = PtE[,2]),
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE, 
                       clear_obs = TRUE)
# Build the mesh
graph$build_mesh(h = 0.05)

##############################
###########Nu 0.5##############
###############################

# Build the model
rspde_model_stat_mean <- rspde.metric_graph(graph,
                                            nu = 0.5,
                                            parameterization = "spde")
# Prepare the data for fitting
data_rspde_bru_stat <- graph_data_rspde(rspde_model_stat_mean,
                                        bru = TRUE)
# Define the component
cmp_stat <- y ~ -1 +
  Intercept(1) +
  field(cbind(.edge_number, .distance_on_edge), model = rspde_model_stat_mean)
# Fit the model
rspde_fit_stat_mean <-
  bru(cmp_stat,
      data = data_rspde_bru_stat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

rspde_fit_stat_mean_nu0.5 <- rspde_fit_stat_mean
save(rspde_fit_stat_mean_nu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_stat_mean_nu0.5.RData"))


##############################
###########Nu 1.5##############
###############################

# Build the model
rspde_model_stat_mean <- rspde.metric_graph(graph,
                                            nu = 1.5,
                                            parameterization = "spde")
# Prepare the data for fitting
data_rspde_bru_stat <- graph_data_rspde(rspde_model_stat_mean,
                                        bru = TRUE)
# Define the component
cmp_stat <- y ~ -1 +
  Intercept(1) +
  field(cbind(.edge_number, .distance_on_edge), model = rspde_model_stat_mean)
# Fit the model
rspde_fit_stat_mean <-
  bru(cmp_stat,
      data = data_rspde_bru_stat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

rspde_fit_stat_mean_nu1.5 <- rspde_fit_stat_mean
save(rspde_fit_stat_mean_nu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_stat_mean_nu1.5.RData"))


##############################
###########Nu est##############
###############################

# Build the model
rspde_model_stat_mean <- rspde.metric_graph(graph,
                                            parameterization = "spde")
# Prepare the data for fitting
data_rspde_bru_stat <- graph_data_rspde(rspde_model_stat_mean,
                                        bru = TRUE)
# Define the component
cmp_stat <- y ~ -1 +
  Intercept(1) +
  field(cbind(.edge_number, .distance_on_edge), model = rspde_model_stat_mean)
# Fit the model
rspde_fit_stat_mean <-
  bru(cmp_stat,
      data = data_rspde_bru_stat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

rspde_fit_stat_mean_nuest <- rspde_fit_stat_mean
save(rspde_fit_stat_mean_nuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_stat_mean_nuest.RData"))

param_spde <- summary(rspde.result(rspde_fit_stat_mean_nuest, "field", rspde_model_stat_mean, parameterization = "spde"))
param_matern <- summary(rspde.result(rspde_fit_stat_mean_nuest, "field", rspde_model_stat_mean, parameterization = "matern"))


data <- graph$get_data()

distance = seq(from = 0, to = 10, by = 0.1)



# Define aux data frame to compute the distance matrix
aux <- data |> 
  rename(distance_on_edge = .distance_on_edge, edge_number = .edge_number) |> # Rename the variables (because graph$compute_geodist_PtE() requires so)
  as.data.frame() |> # Transform to a data frame (i.e., remove the metric_graph class)
  dplyr::select(edge_number, distance_on_edge)

# Compute the distance matrix
distmatrix <- graph$compute_geodist_PtE(PtE = aux,
                                        normalized = TRUE,
                                        include_vertices = FALSE)

# Compute the groups for one replicate
GROUPS <- list()
for (j in 1:length(distance)) {
  GROUPS[[j]] = list()
  for (i in 1:nrow(aux)) {
    GROUPS[[j]][[i]] <- which(as.vector(distmatrix[i, ]) <= distance[j])
  }
}



mse.statnu0.5 <- ls.statnu0.5 <- rep(0,length(distance))
mse.statnu1.5 <- ls.statnu1.5 <- rep(0,length(distance))
mse.statnuest <- ls.statnuest <- rep(0,length(distance))

# cross-validation for-loop
for (j in 1:length(distance)) {
  print(j)
  # cross-validation of the stationary model
  cv.statnu0.5 <- inla.group.cv(rspde_fit_stat_mean_nu0.5, groups = GROUPS[[j]])
  cv.statnu1.5 <- inla.group.cv(rspde_fit_stat_mean_nu1.5, groups = GROUPS[[j]])
  cv.statnuest <- inla.group.cv(rspde_fit_stat_mean_nuest, groups = GROUPS[[j]])
  # obtain MSE and LS
  mse.statnu0.5[j] <- mean((cv.statnu0.5$mean - data$y)^2)
  mse.statnu1.5[j] <- mean((cv.statnu1.5$mean - data$y)^2)
  mse.statnuest[j] <- mean((cv.statnuest$mean - data$y)^2)
  
  
  ls.statnu0.5[j] <- mean(log(cv.statnu0.5$cv))
  ls.statnu1.5[j] <- mean(log(cv.statnu1.5$cv))
  ls.statnuest[j] <- mean(log(cv.statnuest$cv))
}

# Create data frames
mse_df <- data.frame(
  distance,
  Statnu0.5 = mse.statnu0.5,
  Nonstatnu0.5 = mse.statnu0.5,
  Statnu1.5 = mse.statnu1.5,
  Nonstatnu1.5 = mse.statnu1.5,
  Statnuest = mse.statnuest,
  Nonstatnuest = mse.statnuest
)

ls_df <- data.frame(
  distance,
  Statnu0.5 = -ls.statnu0.5,
  Nonstatnu0.5 = -ls.statnu0.5,
  Statnu1.5 = -ls.statnu1.5,
  Nonstatnu1.5 = -ls.statnu1.5,
  Statnuest = -ls.statnuest,
  Nonstatnuest = -ls.statnuest
)

# Save the results



choose_index <- seq(2, nrow(mse_df), by = 3)
mse_df_red <- mse_df[choose_index,]
ls_df_red <- ls_df[choose_index,]
# Convert to long format
mse_long <- mse_df_red %>%
  pivot_longer(cols = -distance, names_to = "nu", values_to = "MSE")

ls_long <- ls_df_red %>%
  pivot_longer(cols = -distance, names_to = "nu", values_to = "LogScore")


# Update the label mappings with the new legend title
label_mapping <- c(
  "Statnu0.5" = "1", 
  "Nonstatnu0.5" = "1", 
  "Statnu1.5" = "2", 
  "Nonstatnu1.5" = "2", 
  "Statnuest" = paste(round(param_matern["nu", "mean"]+0.5, 3), "(est)"), 
  "Nonstatnuest" = paste(round(param_matern["nu", "mean"]+0.5, 3), "(est)")
)

# Define color and linetype mapping
color_mapping <- c(
  "Statnu0.5" = "blue", 
  "Nonstatnu0.5" = "blue", 
  "Statnu1.5" = "black", 
  "Nonstatnu1.5" = "black", 
  "Statnuest" = "red", 
  "Nonstatnuest" = "red"
)

linetype_mapping <- c(
  "Statnu0.5" = "dotdash", 
  "Nonstatnu0.5" = "solid", 
  "Statnu1.5" = "dotdash", 
  "Nonstatnu1.5" = "solid", 
  "Statnuest" = "dotdash", 
  "Nonstatnuest" = "solid"
)

# Plot MSE
mse_plot <- ggplot(mse_long, aes(x = distance, y = MSE, color = nu, linetype = nu)) +
  geom_line(linewidth = 2) +
  labs(y = "MSE", x = "$\\mbox{Geodesic distance } R\\mbox{ }(\\mbox{km})$") +
  scale_color_manual(values = color_mapping, labels = label_mapping, name = "$\\alpha$") +
  scale_linetype_manual(values = linetype_mapping, labels = label_mapping, name = "$\\alpha$") +
  theme_minimal() +
  theme(text = element_text(family = "Palatino"))

# Plot negative log-score
ls_plot <- ggplot(ls_long, aes(x = distance, y = LogScore, color = nu, linetype = nu)) +
  geom_line(linewidth = 2) +
  labs(y = "Negative Log-Score", x = "$\\mbox{Geodesic distance } R\\mbox{ }(\\mbox{km})$") +
  scale_color_manual(values = color_mapping, labels = label_mapping, name = "$\\alpha$") +
  scale_linetype_manual(values = linetype_mapping, labels = label_mapping, name = "$\\alpha$") +
  theme_minimal() +
  theme(text = element_text(family = "Palatino"))

# Combine plots with a shared legend at the top in a single line
new_combined_plot_pems <- mse_plot + ls_plot + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'right') & 
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(nrow = 1))

new_combined_plot_pems

# Save combined plot
# ggsave(here("data_files/crossval_pems.png"), plot = combined_plot_pems, width = 9.22, height = 4.01, dpi = 500)
#myggsave(new_combined_plot_pems, width = 9.22, height = 4.01, use_mathpazo = FALSE)




