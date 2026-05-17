


```{r}
dist_index <- 2

df <- data.frame(
  x = rep(seq_along(GROUPS[[dist_index]]),
          lengths(GROUPS[[dist_index]])),
  y = unlist(GROUPS[[dist_index]])
)

check_neig <- ggplot(df, aes(x, y)) + geom_point()
check_neig <- ggplotly(check_neig)
```

### Using `inla.group.cv()`

```{r, eval = FALSE, echo = FALSE}
j = 2
cv.statnu0.5 <- inla.group.cv(rspde_fit_statnu0.5, groups = GROUPS[[j]])
cv.statnu1.5 <- inla.group.cv(rspde_fit_statnu1.5, groups = GROUPS[[j]])
cv.statnuest <- inla.group.cv(rspde_fit_statnuest, groups = GROUPS[[j]])
cv.nonstatnu0.5 <- inla.group.cv(rspde_fit_nonstatnu0.5, groups = GROUPS[[j]])
cv.nonstatnu1.5 <- inla.group.cv(rspde_fit_nonstatnu1.5, groups = GROUPS[[j]])
cv.nonstatnuest <- inla.group.cv(rspde_fit_nonstatnuest, groups = GROUPS[[j]])

mse.statnu0.5 <- mean((cv.statnu0.5$mean - data$y)^2)
mse.statnu1.5 <- mean((cv.statnu1.5$mean - data$y)^2)
mse.statnuest <- mean((cv.statnuest$mean - data$y)^2)
mse.nonstatnu0.5 <- mean((cv.nonstatnu0.5$mean - data$y)^2)
mse.nonstatnu1.5 <- mean((cv.nonstatnu1.5$mean - data$y)^2)
mse.nonstatnuest <- mean((cv.nonstatnuest$mean - data$y)^2)


MSE_LOO <- data.frame(
  mse.statnu0.5 = mse.statnu0.5,
  mse.statnu1.5 = mse.statnu1.5,
  mse.statnuest = mse.statnuest,
  mse.nonstatnu0.5 = mse.nonstatnu0.5,
  mse.nonstatnu1.5 = mse.nonstatnu1.5,
  mse.nonstatnuest = mse.nonstatnuest
)

save(MSE_LOO, file = here("data_files/new_mse_loo_results.RData"))
```


### Using `rSPDE::cross_validation()`


```{r, eval = FALSE, echo = FALSE}
rSPDE_cv_mse_statnu0.5 <- rSPDE::cross_validation(models = rspde_fit_statnu0.5)
rSPDE_cv_mse_statnu1.5 <- rSPDE::cross_validation(models = rspde_fit_statnu1.5)
rSPDE_cv_mse_statnuest <- rSPDE::cross_validation(models = rspde_fit_statnuest)
rSPDE_cv_mse_nonstatnu0.5 <- rSPDE::cross_validation(models = rspde_fit_nonstatnu0.5)
rSPDE_cv_mse_nonstatnu1.5 <- rSPDE::cross_validation(models = rspde_fit_nonstatnu1.5)
rSPDE_cv_mse_nonstatnuest <- rSPDE::cross_validation(models = rspde_fit_nonstatnuest)
mse_CV <- data.frame(
  cv_mse_statnu0.5 = rSPDE_cv_mse_statnu0.5$mse[1],
  cv_mse_statnu1.5 = rSPDE_cv_mse_statnu1.5$mse[1],
  cv_mse_statnuest = rSPDE_cv_mse_statnuest$mse[1],
  cv_mse_nonstatnu0.5 = rSPDE_cv_mse_nonstatnu0.5$mse[1],
  cv_mse_nonstatnu1.5 = rSPDE_cv_mse_nonstatnu1.5$mse[1],
  cv_mse_nonstatnuest = rSPDE_cv_mse_nonstatnuest$mse[1]
)
save(mse_CV, file = here("data_files/new_mse_cv_results.RData"))
```

```{r, eval = FALSE}
rSPDE_cv_mse_statnu0.5 <- rSPDE::cross_validation(models = rspde_fit_statnu0.5, return_train_test = TRUE)
train_test_aux <- rSPDE_cv_mse_statnu0.5$train_test
MG_cv_mse_statnu0.5 <- MetricGraph::cross_validation(models = rspde_fit_statnu0.5, train_test_indexes = train_test_aux)
```


### Using `MetricGraph::cross_validation()`


```{r, eval = FALSE, echo = FALSE}
library(slackr)
source("keys.R")
slackr_setup(token = token) # token comes from keys.R

output <- capture.output({
  withCallingHandlers({cv_mse_statnu0.5 <- MetricGraph::cross_validation(models = rspde_fit_statnu0.5)},
                      message = function(m) {slackr::slackr_msg(text = conditionMessage(m), channel = "#research")},
                      warning = function(w) {slackr::slackr_msg(text = conditionMessage(w), channel = "#research")})
})
slackr::slackr_msg(text = paste(output, collapse = "\n"),channel = "#research")
save(cv_mse_statnu0.5, file = here("data_files/mse_CV_from_MetricGraph.RData"))

output <- capture.output({
  withCallingHandlers({cv_mse_statnu1.5 <- MetricGraph::cross_validation(models = rspde_fit_statnu1.5)},
                      message = function(m) {slackr::slackr_msg(text = conditionMessage(m), channel = "#research")},
                      warning = function(w) {slackr::slackr_msg(text = conditionMessage(w), channel = "#research")})
})
slackr::slackr_msg(text = paste(output, collapse = "\n"),channel = "#research")
save(cv_mse_statnu1.5, file = here("data_files/mse_CV_from_MetricGraph.RData"))

output <- capture.output({
  withCallingHandlers({cv_mse_statnuest <- MetricGraph::cross_validation(models = rspde_fit_statnuest)},
                      message = function(m) {slackr::slackr_msg(text = conditionMessage(m), channel = "#research")},
                      warning = function(w) {slackr::slackr_msg(text = conditionMessage(w), channel = "#research")})
})
slackr::slackr_msg(text = paste(output, collapse = "\n"),channel = "#research")
save(cv_mse_statnuest, file = here("data_files/mse_CV_from_MetricGraph.RData"))

output <- capture.output({
  withCallingHandlers({cv_mse_nonstatnu0.5 <- MetricGraph::cross_validation(models = rspde_fit_nonstatnu0.5)},
                      message = function(m) {slackr::slackr_msg(text = conditionMessage(m), channel = "#research")},
                      warning = function(w) {slackr::slackr_msg(text = conditionMessage(w), channel = "#research")})
})
slackr::slackr_msg(text = paste(output, collapse = "\n"),channel = "#research")
save(cv_mse_nonstatnu0.5, file = here("data_files/mse_CV_from_MetricGraph.RData"))

output <- capture.output({
  withCallingHandlers({cv_mse_nonstatnu1.5 <- MetricGraph::cross_validation(models = rspde_fit_nonstatnu1.5)},
                      message = function(m) {slackr::slackr_msg(text = conditionMessage(m), channel = "#research")},
                      warning = function(w) {slackr::slackr_msg(text = conditionMessage(w), channel = "#research")})
})
slackr::slackr_msg(text = paste(output, collapse = "\n"),channel = "#research")
save(cv_mse_nonstatnu1.5, file = here("data_files/mse_CV_from_MetricGraph.RData"))

output <- capture.output({
  withCallingHandlers({cv_mse_nonstatnuest <- MetricGraph::cross_validation(models = rspde_fit_nonstatnuest)},
                      message = function(m) {slackr::slackr_msg(text = conditionMessage(m), channel = "#research")},
                      warning = function(w) {slackr::slackr_msg(text = conditionMessage(w), channel = "#research")})
})
slackr::slackr_msg(text = paste(output, collapse = "\n"),channel = "#research")
save(cv_mse_nonstatnuest, file = here("data_files/mse_CV_from_MetricGraph.RData"))

mse_CV_from_MetricGraph <- data.frame(
  cv_mse_statnu0.5 = cv_mse_statnu0.5$mse[1],
  cv_mse_statnu1.5 = cv_mse_statnu1.5$mse[1],
  cv_mse_statnuest = cv_mse_statnuest$mse[1],
  cv_mse_nonstatnu0.5 = cv_mse_nonstatnu0.5$mse[1],
  cv_mse_nonstatnu1.5 = cv_mse_nonstatnu1.5$mse[1],
  cv_mse_nonstatnuest = cv_mse_nonstatnuest$mse[1]
)
save(mse_CV_from_MetricGraph, file = here::here("data_files/mse_CV_from_MetricGraph.RData"))
```



# Linear regression and kNN regression



```{r, eval = FALSE}
# Load the data

load(here::here("data_files/new_Y_mean.RData")) # was created in pems1.Rmd

Y_mu <- apply(Y_raw[1:13,], 2, mean)

data_simple <- data.frame(y = c(t(Y_raw[14:26,])), 
                          mean_value = rep(Y_mu, times = 13), 
                          repl = rep(1:13, each = 314))

data_simple$repl <- factor(data_simple$repl)
```


```{r, eval = FALSE}
library(lme4)
n <- nrow(data_simple)

pred_loocv <- numeric(n)
#pred_loocv_repl <- numeric(n)

for(i in 1:n){
  train_data <- data_simple[-i, ]
  test_data  <- data_simple[i, , drop = FALSE]
  
  model <- lm(y ~ mean_value, data = train_data)
  pred_loocv[i] <- predict(model, newdata = test_data)
  
  # model <- lmer(y ~ mean_value + (1 | repl), data = train_data)
  # pred_loocv_repl[i] <- predict(model, newdata = test_data, re.form = NULL)
  print(paste("Processed observation", i, "out of", n))
}

mse_loocv_lm <- mean((data_simple$y - pred_loocv)^2)
mse_loocv_lm

# mse_loocv_lmer <- mean((data_simple$y - pred_loocv_repl)^2)
# mse_loocv_lmer
```


```{r, eval = FALSE}
# Load the data
load(here::here("data_files/new_pems_repl1_data.RData"))
# Extract the data from the graph

initial_data <- graph$get_data()
data <- initial_data |> as.data.frame() |> select(y, mean_value, repl)

n <- length(data |> filter(repl == 1) |> pull(y))
data$repl <- factor(rep(1:13, each = n))

```

```{r, eval = FALSE}
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


# Another approach
# 
# k_values <- 1:30
# 
# loo_mse <- sapply(k_values, function(k) {
#   n <- nrow(data)
#   pred <- numeric(n)
#   
#   for (i in 1:n) {
#     pred[i] <- knn.reg(
#       train = matrix(data$mean_value[-i]),
#       y     = data$y[-i],
#       test  = matrix(data$mean_value[i]),
#       k     = k
#     )$pred
#   }
#   
#   mean((data$y - pred)^2)
# })
# 
# plot(k_values, loo_mse, type = "b", xlab = "k", ylab = "LOO RMSE")
# 
# best_k    <- k_values[which.min(loo_mse)]
# best_knn_mse <- loo_mse[which.min(loo_mse)]


save(mse_loocv_lm, best_KNN_mse, file = here("data_files/new_simple_linear_regression_results.RData"))
```

Below we show the MSE results for the simple linear regression `mse_loocv_lm` and the kNN regression `best_KNN_mse`.


```{r, eval = TRUE}
load(here::here("data_files/new_simple_linear_regression_results.RData"))
load(here::here("data_files/new_mse_loo_results.RData"))
#load(here::here("data_files/new_mse_cv_results.RData"))
load(here::here("data_files/new_MSE_ISOCOV.RData"))
#load(here::here("data_files/mse_CV_from_MetricGraph.RData"))
#mse_CV
#mse_CV_from_MetricGraph
MSE_LOO
data.frame(LM = mse_loocv_lm, 
           kNNdistMAT = best_KNN_mse,
           MSE_ISOCOV = MSE_ISOCOV)
```

