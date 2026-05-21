EXPNUM <- 1

## ---------------------------
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

inla.setOption(num.threads = parallel::detectCores() - 4) #20


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

process_model_results <- function(fit, model) {
  fit_spde <- rspde.result(fit, "field", model, parameterization = "spde")
  fit_matern <- rspde.result(fit, "field", model, parameterization = "matern")
  df_for_plot_spde <- gg_df(fit_spde)
  df_for_plot_matern <- gg_df(fit_matern)
  param_spde <- summary(fit_spde)
  param_matern <- summary(fit_matern)
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
  return(list(allparams = allparams, df_for_plot_spde = df_for_plot_spde, df_for_plot_matern = df_for_plot_matern))
}


# Load the data
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/new_pems_repl1_data.RData"))
load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/tau_from_graphlme.RData"))

load(here::here("data_files/new_pems_repl1_data.RData"))
load(here::here("data_files/tau_from_graphlme.RData"))

# Non-stationary parameters
B.tau =   cbind(0, 1, 0, cov, 0)
B.kappa = cbind(0, 0, 1, 0, cov)


log_tau_from_graphlme <- log(tau_from_graphlme)
log_kappa_from_graphlme <- log(kappa_from_graphlme)

#####################################
#############nu=0.5##################
#####################################

# Build the model
rspde_model_stat <- rspde.metric_graph(graph,
                                       start.ltau = log_tau_from_graphlme,
                                       start.lkappa = log_kappa_from_graphlme,
                                       parameterization = "spde",
                                       nu = 0.5)
# Prepare the data for fitting
data_rspde_bru_stat <- graph_data_rspde(rspde_model_stat,
                                        repl = ".all",
                                        bru = TRUE,
                                        repl_col = "repl")
# Define the component
cmp_stat <- y ~ -1 +
  Intercept(1) +
  mean_value +
  field(cbind(.edge_number, .distance_on_edge), 
        model = rspde_model_stat,
        replicate = repl)
# Fit the model
rspde_fit_stat <-
  bru(cmp_stat,
      data = data_rspde_bru_stat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

output_from_models <-process_model_results(rspde_fit_stat, rspde_model_stat)
parameters_statistics_statnu0.5 <- output_from_models$allparams
rspde_fit_statnu0.5 <- rspde_fit_stat
save(rspde_fit_statnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_statnu0.5.RData"))
save(parameters_statistics_statnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_statnu0.5.RData"))
slackr_msg(
  text = paste(
    capture.output(print(parameters_statistics_statnu0.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_statnu0.5.RData"))
parameters_statistics_statnu0.5

start.theta <- c(log(parameters_statistics_statnu0.5["tau","mean"]), 
                 log(parameters_statistics_statnu0.5["kappa","mean"]), 
                 rep(0, (ncol(B.tau)-3)))



# Build the model
rspde_model_nonstat <- rspde.metric_graph(graph,
                                          start.theta = start.theta,
                                          B.tau = B.tau,
                                          B.kappa =  B.kappa,
                                          parameterization = "spde",
                                          nu = 0.5)
# Prepare the data for fitting
data_rspde_bru_nonstat <- graph_data_rspde(rspde_model_nonstat,
                                           repl = ".all",
                                           bru = TRUE,
                                           repl_col = "repl")
# Define the component
cmp_nonstat <- y ~ -1 +
  Intercept(1) +
  mean_value +
  field(cbind(.edge_number, .distance_on_edge), 
        model = rspde_model_nonstat,
        replicate = repl)
# Fit the model
rspde_fit_nonstat <-
  bru(cmp_nonstat,
      data = data_rspde_bru_nonstat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

output_from_models <- process_model_results(rspde_fit_nonstat, rspde_model_nonstat)
parameters_statistics_nonstatnu0.5 <- output_from_models$allparams
rspde_fit_nonstatnu0.5 <- rspde_fit_nonstat
save(rspde_fit_nonstatnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_nonstatnu0.5.RData"))
save(parameters_statistics_nonstatnu0.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_nonstatnu0.5.RData"))
slackr_msg(
  text = paste(
    capture.output(print(parameters_statistics_nonstatnu0.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_nonstatnu0.5.RData"))
parameters_statistics_nonstatnu0.5


#####################################
#############nu=1.5##################
#####################################


# Build the model
rspde_model_stat <- rspde.metric_graph(graph,
                                       start.ltau = log_tau_from_graphlme,
                                       start.lkappa = log_kappa_from_graphlme,
                                       parameterization = "spde",
                                       nu = 1.5)
# Prepare the data for fitting
data_rspde_bru_stat <- graph_data_rspde(rspde_model_stat,
                                        repl = ".all",
                                        bru = TRUE,
                                        repl_col = "repl")
# Define the component
cmp_stat <- y ~ -1 +
  Intercept(1) +
  mean_value +
  field(cbind(.edge_number, .distance_on_edge), 
        model = rspde_model_stat,
        replicate = repl)
# Fit the model
rspde_fit_stat <-
  bru(cmp_stat,
      data = data_rspde_bru_stat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

output_from_models <-process_model_results(rspde_fit_stat, rspde_model_stat)
parameters_statistics_statnu1.5 <- output_from_models$allparams
rspde_fit_statnu1.5 <- rspde_fit_stat
save(rspde_fit_statnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_statnu1.5.RData"))
save(parameters_statistics_statnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_statnu1.5.RData"))
slackr_msg(
  text = paste(
    capture.output(print(parameters_statistics_statnu1.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_statnu1.5.RData"))
parameters_statistics_statnu1.5

start.theta <- c(log(parameters_statistics_statnu1.5["tau","mean"]), 
                 log(parameters_statistics_statnu1.5["kappa","mean"]), 
                 rep(0, (ncol(B.tau)-3)))


# Build the model
rspde_model_nonstat <- rspde.metric_graph(graph,
                                          start.theta = start.theta,
                                          B.tau = B.tau,
                                          B.kappa =  B.kappa,
                                          parameterization = "spde",
                                          nu = 1.5)
# Prepare the data for fitting
data_rspde_bru_nonstat <- graph_data_rspde(rspde_model_nonstat,
                                           repl = ".all",
                                           bru = TRUE,
                                           repl_col = "repl")
# Define the component
cmp_nonstat <- y ~ -1 +
  Intercept(1) +
  mean_value +
  field(cbind(.edge_number, .distance_on_edge), 
        model = rspde_model_nonstat,
        replicate = repl)
# Fit the model
rspde_fit_nonstat <-
  bru(cmp_nonstat,
      data = data_rspde_bru_nonstat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

output_from_models <- process_model_results(rspde_fit_nonstat, rspde_model_nonstat)
parameters_statistics_nonstatnu1.5 <- output_from_models$allparams
rspde_fit_nonstatnu1.5 <- rspde_fit_nonstat
save(rspde_fit_nonstatnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_nonstatnu1.5.RData"))
save(parameters_statistics_nonstatnu1.5, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_nonstatnu1.5.RData"))
slackr_msg(
  text = paste(
    capture.output(print(parameters_statistics_nonstatnu1.5)),
    collapse = "\n"
  ),
  channel = "#research"
)

load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_nonstatnu1.5.RData"))
parameters_statistics_nonstatnu1.5

#####################################
#############nu=est##################
#####################################



# Build the model
rspde_model_stat <- rspde.metric_graph(graph,
                                       start.ltau = log_tau_from_graphlme,
                                       start.lkappa = log_kappa_from_graphlme,
                                       parameterization = "spde")
# Prepare the data for fitting
data_rspde_bru_stat <- graph_data_rspde(rspde_model_stat,
                                        repl = ".all",
                                        bru = TRUE,
                                        repl_col = "repl")
# Define the component
cmp_stat <- y ~ -1 +
  Intercept(1) +
  mean_value +
  field(cbind(.edge_number, .distance_on_edge), 
        model = rspde_model_stat,
        replicate = repl)
# Fit the model
rspde_fit_stat <-
  bru(cmp_stat,
      data = data_rspde_bru_stat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

output_from_models <-process_model_results(rspde_fit_stat, rspde_model_stat)
parameters_statistics_statnuest <- output_from_models$allparams
rspde_fit_statnuest <- rspde_fit_stat
save(rspde_fit_statnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_statnuest.RData"))
save(parameters_statistics_statnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_statnuest.RData"))
alpha_stat <- parameters_statistics_statnuest[5,1] + 0.5
save(alpha_stat, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/alpha_stat.RData"))
slackr_msg(
  text = paste(
    capture.output(print(parameters_statistics_statnuest)),
    collapse = "\n"
  ),
  channel = "#research"
)

load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_statnuest.RData"))
parameters_statistics_statnuest

start.theta <- c(log(parameters_statistics_statnuest["tau","mean"]), 
                 log(parameters_statistics_statnuest["kappa","mean"]), 
                 rep(0, (ncol(B.tau)-3)))

# Build the model
rspde_model_nonstat <- rspde.metric_graph(graph,
                                          start.nu = parameters_statistics_statnuest[5,1],
                                          start.theta = start.theta,
                                          B.tau = B.tau,
                                          B.kappa =  B.kappa,
                                          parameterization = "spde")
# Prepare the data for fitting
data_rspde_bru_nonstat <- graph_data_rspde(rspde_model_nonstat,
                                           repl = ".all",
                                           bru = TRUE,
                                           repl_col = "repl")
# Define the component
cmp_nonstat <- y ~ -1 +
  Intercept(1) +
  mean_value +
  field(cbind(.edge_number, .distance_on_edge), 
        model = rspde_model_nonstat,
        replicate = repl)
# Fit the model
rspde_fit_nonstat <-
  bru(cmp_nonstat,
      data = data_rspde_bru_nonstat[["data"]],
      family = "gaussian",
      options = list(verbose = FALSE)
  )

output_from_models <- process_model_results(rspde_fit_nonstat, rspde_model_nonstat)
parameters_statistics_nonstatnuest <- output_from_models$allparams
rspde_fit_nonstatnuest <- rspde_fit_nonstat
save(rspde_fit_nonstatnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/rspde_fit_nonstatnuest.RData"))
save(parameters_statistics_nonstatnuest, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_nonstatnuest.RData"))
alpha_nonstat <- parameters_statistics_nonstatnuest[7,1] + 0.5
save(alpha_nonstat, file = paste0("~/Desktop/folder_aux/exp", EXPNUM, "/alpha_nonstat.RData"))
slackr_msg(
  text = paste(
    capture.output(print(parameters_statistics_nonstatnuest)),
    collapse = "\n"
  ),
  channel = "#research"
)

load(paste0("~/Desktop/folder_aux/exp", EXPNUM, "/parameters_statistics_nonstatnuest.RData"))
parameters_statistics_nonstatnuest








