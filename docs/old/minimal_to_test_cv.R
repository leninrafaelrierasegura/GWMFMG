#libraries
library(MetricGraph)
library(rSPDE)
library(inlabru)
library(dplyr)

set.seed(1982)
# Build the graph
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edges = list(edge1, edge2, edge3, edge4)
graph <- metric_graph$new(edges = edges)


# Build the mesh on the graph
graph$build_mesh(h = 0.01)

# Parameters
sigma <- 1.2
range <- 0.2
nu <- 0.8
rspde.order <- 1

kappa <- sqrt(8*nu)/range
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))  #sigma = 1, d = 1

# Build the operator
op <- matern.operators(nu = nu, range = range, sigma = sigma, 
                       parameterization = "matern",
                       m = rspde.order, graph = graph) 

# Set the number of replicates and simulate
n.rep <- 2
u.rep <- simulate(op, nsim = n.rep)

obs.per.edge <- 20
obs.loc <- NULL
for(i in 1:graph$nE) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge)))
}
n.obs <- obs.per.edge*graph$nE
A <- graph$fem_basis(obs.loc)

sigma.e <- 0.1
Y.rep <- A %*% u.rep + sigma.e * matrix(rnorm(n.obs * n.rep), ncol = n.rep)
y_vec <- 1 + as.vector(Y.rep)
repl <- rep(1:n.rep, each = n.obs)                       

df_data_repl <- data.frame(y = y_vec,
                           edge_number = rep(obs.loc[,1], n.rep),
                           distance_on_edge = rep(obs.loc[,2], n.rep), 
                           repl = repl)

graph$add_observations(data = df_data_repl, normalized = TRUE, 
                       group = "repl", clear_obs = TRUE)

############### Case alpha = 1 ###########################
# Build the model
rspde_model_stat <- rspde.metric_graph(graph,
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
rspde_fit_statnu0.5 <- rspde_fit_stat



rSPDE_cv_mse_statnu0.5 <- rSPDE::cross_validation(
  models = rspde_fit_statnu0.5, 
  return_train_test = TRUE,
  true_CV = FALSE)

train_test_aux <- rSPDE_cv_mse_statnu0.5$train_test

MG_cv_mse_statnu0.5 <- MetricGraph::cross_validation(
  models = rspde_fit_statnu0.5,
  train_test_indexes = train_test_aux,
  true_CV = FALSE)


rSPDE_cv_mse_statnu0.5$scores_df$mse
MG_cv_mse_statnu0.5$mse








