# ## ---- 1. build graph & mesh exactly like the user's example -----------------
# data(pems)
# graph <- metric_graph$new(edges = pems$edges)
# graph$add_observations(data = pems$data, verbose = 0)
# graph$prune_vertices(verbose = 0)
# graph$build_mesh(h = 0.1)
# 
# ## covariate at mesh nodes (used to build kappa/tau at FEM nodes)
# cov_pos_mesh <- as.numeric(
#   (graph$mesh$VtE[, 2] > 0.9) | (graph$mesh$VtE[, 2] < 0.1)
# )
# 
# 
# 
# 
# ## same B matrices as the example
# B.tau <- cbind(0, 1, 0, cov_pos_mesh, 0)
# B.kappa <- cbind(0, 0, 1, 0, cov_pos_mesh)
# 



# Load the data
load(here("data_files/new_pems_repl1_data.RData"))
# graph <- update_graph(graph)
# Extract the data from the graph
data <- graph$get_data()

data_mm <- graph$get_PtE()

B.tau = cbind(0, 1, 0, cov, 0)
B.kappa = cbind(0, 0, 1, 0, cov)



## ---- 2. pick TRUE parameters and build the rSPDE operator ------------------
true_theta <- c(theta1 = log(1.0),   # baseline log-sigma
                theta2 = log(0.5),   # baseline log-range
                theta3 = log(2.0),   # multiplicative bump on sigma at endpoints
                theta4 = log(0.5))   # multiplicative shrink on range at endpoints
true_alpha   <- 1.5            # nu = 1
true_nu      <- true_alpha - 0.5
true_sigma_e <- 0.1
true_beta0   <- 2.0




op_true <- rSPDE::spde.matern.operators(
  graph            = graph,
  B.tau          = B.tau,
  B.kappa         = B.kappa,
  parameterization = "spde",
  alpha               = true_alpha,
  theta            = true_theta,
  m                = 2
)


## ---- 3. simulate the latent field on the mesh and project to obs locs -----
u_mesh <- simulate(op_true, nsim = 13)          # one column

data_priv <- graph$.__enclos_env__$private$data
obs_pte <- data_mm
A_obs <- graph$fem_basis(obs_pte)

u_obs <- as.numeric(A_obs %*% u_mesh)
n_obs <- length(u_obs)

cov_obs <- data$mean_value

y_sim <- as.numeric(true_beta0 + u_obs + cov_obs + true_sigma_e * rnorm(n_obs))

## also build the covariate at observation locations (same definition)
# cov_obs <- data_priv[[".distance_on_edge"]]
# cov_obs <- as.numeric((cov_obs > 0.9) | (cov_obs < 0.1))



## ---- 4. push simulated data into the graph and fit -------------------------
sim_df <- data.frame(y = y_sim, cov_obs = cov_obs)
graph$add_observations(data = graph$mutate(y       = y_sim,
                                                     cov_obs = cov_obs),
                            clear_obs = TRUE,
                            verbose   = 0)

t0 <- Sys.time()

fit <- graph_lme(y ~ cov_obs,
                 graph = graph,
                 parallel = TRUE,
                 which_repl = 2,
                 model = list(type = "WhittleMatern",
                              B.tau = B.tau,
                              B.kappa = B.kappa,
                              #alpha   = true_alpha,
                              fem     = TRUE))

res_exp <- graph_lme(y ~ cov_obs, 
                     graph = graph, 
                     parallel = TRUE,
                     which_repl = 2,
                     model = list(type = "isoCov"))

t0 <- Sys.time()
POST1 <- posterior_crossvalidation(object = fit, mode = "loo", true_CV = FALSE)
t1 <- Sys.time()

t2 <- Sys.time()
POST2 <- posterior_crossvalidation_loo(object = res_exp)
t3 <- Sys.time()

cat("graph_lme fit time:", format(Sys.time() - t0), "\n")

cat("\n================ FIT SUMMARY ================\n")
print(summary(fit))

cat("\n================ COMPARISON ================\n")
re <- fit$coeff$random_effects
fe <- fit$coeff$fixed_effects
me <- fit$coeff$measurement_error

truth_vec <- c(true_theta[1], true_theta[2], true_theta[3], true_theta[4])
names(truth_vec) <- c("theta1", "theta2", "theta3", "theta4")

cat(sprintf("%-12s %12s %12s %12s\n", "param", "true", "estimate", "diff"))
for (nm in names(truth_vec)) {
  if (nm %in% names(re)) {
    est <- as.numeric(re[nm])
    cat(sprintf("%-12s %12.4f %12.4f %12.4f\n", nm, truth_vec[nm],
                est, est - truth_vec[nm]))
  } else {
    cat(sprintf("%-12s %12.4f %12s\n", nm, truth_vec[nm], "<missing>"))
  }
}
# nu (should be fixed at true_nu)
if ("nu" %in% names(re)) {
  est <- as.numeric(re["nu"])
  cat(sprintf("%-12s %12.4f %12.4f %12.4f (log-scale, exp=%.4f)\n",
              "log(nu)", log(true_nu), est, est - log(true_nu), exp(est)))
}
# sigma_e
if (!is.null(me)) {
  cat(sprintf("%-12s %12.4f %12.4f %12.4f\n", "sigma_e",
              true_sigma_e, as.numeric(me), as.numeric(me) - true_sigma_e))
}
cat(sprintf("%-12s %12.4f %12.4f %12.4f\n", "(Intercept)",
            true_beta0, as.numeric(fe), as.numeric(fe) - true_beta0))

