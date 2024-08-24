library(glmnet)
library(igraph)
library(MASS)
library(lars)
library(Rcpp)
library(RcppArmadillo)


sourceCpp("main_functions.cpp")
source("auxiliary_functions.R")




task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
seed <- task_id
print(seed)
set.seed(seed)


n0 <- 200
K <- 90
n_vec <- sample(1:5, size = K, replace = T) * n0
N <- sum(n_vec)
p <- 100
r <- 0.5



intercept1 <- 1.5
intercept2 <- 1
intercept3 <- 0.5
intercept4 <- -0.5
intercept5 <- -1
intercept6 <- -1.5


coe1 <- c(intercept1, c(c(rnorm(5, mean = 0.4, sd = 0.1), rnorm(5, mean = 0.8, sd = 0.1)) * sign(runif(10, min = -0.5, max = 0.5)), rep(0, p-11)))
coe2 <- c(intercept2, c(c(rnorm(5, mean = 0.4, sd = 0.1), rnorm(5, mean = 0.8, sd = 0.1)) * sign(runif(10, min = -0.5, max = 0.5)), rep(0, p-11)))
coe3 <- c(intercept3, c(c(rnorm(5, mean = 0.4, sd = 0.1), rnorm(5, mean = 0.8, sd = 0.1)) * sign(runif(10, min = -0.5, max = 0.5)), rep(0, p-11)))
coe4 <- c(intercept4, c(c(rnorm(5, mean = 0.4, sd = 0.1), rnorm(5, mean = 0.8, sd = 0.1)) * sign(runif(10, min = -0.5, max = 0.5)), rep(0, p-11)))
coe5 <- c(intercept5, c(c(rnorm(5, mean = 0.4, sd = 0.1), rnorm(5, mean = 0.8, sd = 0.1)) * sign(runif(10, min = -0.5, max = 0.5)), rep(0, p-11)))
coe6 <- c(intercept6, c(c(rnorm(5, mean = 0.4, sd = 0.1), rnorm(5, mean = 0.8, sd = 0.1)) * sign(runif(10, min = -0.5, max = 0.5)), rep(0, p-11)))

theta_true <- matrix(0, nrow = p, ncol = K)
theta_true[1:11,1:(K/9)] <- coe1[1:11]
theta_true[1:11,(K/9+1):(2*K/9)] <- coe2[1:11]
theta_true[1:11,(2*K/9+1):(7*K/18)] <- coe3[1:11]
theta_true[1:11,(7*K/18+1):(5*K/9)] <- coe4[1:11]
theta_true[1:11,(5*K/9+1):(7*K/9)] <- coe5[1:11]
theta_true[1:11,(7*K/9+1):K] <- coe6[1:11]

mu <- rep(0, p-1)
Sigma <- matrix(0, nrow = p-1, ncol = p-1)


for (i in 1:(p-1)) {
  for (j in 1:(p-1)) {
    Sigma[i,j] <- r^(abs(i-j))
  }
}


X_lst <- vector(mode = "list", length = K)
Y_lst <- vector(mode = "list", length = K)


for (k in 1:(K/9)) {
  X <- mvrnorm(n_vec[k], mu = mu, Sigma = Sigma)
  odds <- exp(cbind(rep(1, n_vec[k]), X) %*% coe1)
  p_vec <- odds / (1 + odds)
  Y <- rbinom(n_vec[k], 1, p_vec)
  X_lst[[k]] <- X
  Y_lst[[k]] <- Y
  print(sum(Y_lst[[k]] == 1))
}


for (k in (K/9+1):(2*K/9)) {
  X <- mvrnorm(n_vec[k], mu = mu, Sigma = Sigma)
  odds <- exp(cbind(rep(1, n_vec[k]), X) %*% coe2)
  p_vec <- odds / (1 + odds)
  Y <- rbinom(n_vec[k], 1, p_vec)
  X_lst[[k]] <- X
  Y_lst[[k]] <- Y
  print(sum(Y_lst[[k]] == 1))
}




for (k in (2*K/9+1):(7*K/18)) {
  X <- mvrnorm(n_vec[k], mu = mu, Sigma = Sigma)
  odds <- exp(cbind(rep(1, n_vec[k]), X) %*% coe3)
  p_vec <- odds / (1 + odds)
  Y <- rbinom(n_vec[k], 1, p_vec)
  X_lst[[k]] <- X
  Y_lst[[k]] <- Y
  print(sum(Y_lst[[k]] == 1))
}



for (k in (7*K/18+1):(5*K/9)) {
  X <- mvrnorm(n_vec[k], mu = mu, Sigma = Sigma)
  odds <- exp(cbind(rep(1, n_vec[k]), X) %*% coe4)
  p_vec <- odds / (1 + odds)
  Y <- rbinom(n_vec[k], 1, p_vec)
  X_lst[[k]] <- X
  Y_lst[[k]] <- Y
  print(sum(Y_lst[[k]] == 1))
}



for (k in (5*K/9+1):(7*K/9)) {
  X <- mvrnorm(n_vec[k], mu = mu, Sigma = Sigma)
  odds <- exp(cbind(rep(1, n_vec[k]), X) %*% coe5)
  p_vec <- odds / (1 + odds)
  Y <- rbinom(n_vec[k], 1, p_vec)
  X_lst[[k]] <- X
  Y_lst[[k]] <- Y
  print(sum(Y_lst[[k]] == 1))
}



for (k in (7*K/9+1):K) {
  X <- mvrnorm(n_vec[k], mu = mu, Sigma = Sigma)
  odds <- exp(cbind(rep(1, n_vec[k]), X) %*% coe6)
  p_vec <- odds / (1 + odds)
  Y <- rbinom(n_vec[k], 1, p_vec)
  X_lst[[k]] <- X
  Y_lst[[k]] <- Y
  print(sum(Y_lst[[k]] == 1))
}



################################################################################
##################                  fit the ICR.                ################
################################################################################
system.time(local_summary <- Local_fit(Y_lst, X_lst, family = "binomial")) 
print(n_vec)


tau_lst <- 0.1 * 5:15 * sqrt((log(p)) / mean(n_vec))
length_lst <- c()
for (k in 1:K){
  length_lst <- c(length_lst, length(X_lst[[k]][,1]))
}

system.time(debiased_beta <- debias_fit(local_summary$I, local_summary$U, beta_fit = local_summary$beta, length_lst, tau_lst)) 


lambda1_cand <- seq(0.2, 0.6, 0.05) 
lambda2_cand <- seq(0.2, 0.6, 0.1) 

run.time <- system.time(optimal_tuning <- tuning_selection_cpp(theta_ini = debiased_beta, local_summary$I, local_summary$U, 
                                                               K, p, N, n_vec, lambda1_cand, lambda2_cand, con = 1)) 

print(run.time)


V_tilde_mat <- do.call("cbind", local_summary$I)
zeta_tilde_mat <- do.call("cbind", local_summary$U)
index <- t(combn(K,2))

system.time(result_icr <- ICR_pg_cpp(theta_ini = debiased_beta, V_tilde_mat, zeta_tilde_mat, n_vec, index,
                                     nu = 1, tau = 3, lambda1 = optimal_tuning[1], lambda2 = optimal_tuning[2], 
                                     eps_outer = 1e-3, max_iter = 100)) 







write.csv(lambda1_cand, file = paste('Example5_K90_n200_p500_ICR_lambda1_cand_', seed, '.csv', sep = ""))
write.csv(lambda2_cand, file = paste('Example5_K90_n200_p500_ICR_lambda2_cand_', seed, '.csv', sep = ""))
write.csv(run.time[3], file = paste('Example5_K90_n200_p500_ICR_run.time_', seed, '.csv', sep = ""))
write.csv(result_icr$theta, file = paste('Example5_K90_n200_p500_ICR_theta_', seed, '.csv', sep = ""))
write.csv(result_icr$alpha, file = paste('Example5_K90_n200_p500_ICR_alpha_', seed, '.csv', sep = ""))
write.csv(theta_true, file = paste('Example5_K90_n200_p500_ICR_thetatrue_', seed, '.csv', sep = ""))


