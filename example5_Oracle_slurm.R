library(glmnet)
library(igraph)
library(MASS)
library(lars)
library(Rcpp)
library(RcppArmadillo)




source("auxiliary_functions.R")


task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
seed <- task_id
print(seed)
set.seed(seed)


n0 <- 400
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
##################                  fit the Oracle.             ################
################################################################################
system.time(local_summary <- Local_fit(Y_lst, X_lst, family = "binomial")) 
print(n_vec)


term11 <- matrix(0, nrow = 11, ncol = 11)
for (k in 1:(K/9)) {
  term11 <- term11 + local_summary$I[[k]][1:11, 1:11]
}

term12 <- rep(0, 11)
for (k in 1:(K/9)) {
  term12 <- term12 + local_summary$U[[k]][1:11]
}

coeff1 <- as.numeric(solve(term11) %*% term12)



term21 <- matrix(0, nrow = 11, ncol = 11)
for (k in (K/9+1):(2*K/9)) {
  term21 <- term21 + local_summary$I[[k]][1:11, 1:11]
}
term22 <- rep(0, 11)
for (k in (K/9+1):(2*K/9)) {
  term22 <- term22 + local_summary$U[[k]][1:11]
}

coeff2 <- as.numeric(solve(term21) %*% term22)



term31 <- matrix(0, nrow = 11, ncol = 11)
for (k in (2*K/9+1):(7*K/18)) {
  term31 <- term31 + local_summary$I[[k]][1:11, 1:11]
}

term32 <- rep(0, 11)
for (k in (2*K/9+1):(7*K/18)) {
  term32 <- term32 + local_summary$U[[k]][1:11]
}

coeff3 <- as.numeric(solve(term31) %*% term32)


term41 <- matrix(0, nrow = 11, ncol = 11)
for (k in (7*K/18+1):(5*K/9)) {
  term41 <- term41 + local_summary$I[[k]][1:11, 1:11]
}

term42 <- rep(0, 11)
for (k in (7*K/18+1):(5*K/9)) {
  term42 <- term42 + local_summary$U[[k]][1:11]
}

coeff4 <- as.numeric(solve(term41) %*% term42)



term51 <- matrix(0, nrow = 11, ncol = 11)
for (k in (5*K/9+1):(7*K/9)) {
  term51 <- term51 + local_summary$I[[k]][1:11, 1:11]
}

term52 <- rep(0, 11)
for (k in (5*K/9+1):(7*K/9)) {
  term52 <- term52 + local_summary$U[[k]][1:11]
}

coeff5 <- as.numeric(solve(term51) %*% term52)




term61 <- matrix(0, nrow = 11, ncol = 11)
for (k in (7*K/9+1):K) {
  term61 <- term61 + local_summary$I[[k]][1:11, 1:11]
}

term62 <- rep(0, 11)
for (k in (7*K/9+1):K) {
  term62 <- term62 + local_summary$U[[k]][1:11]
}

coeff6 <- as.numeric(solve(term61) %*% term62)



theta_ora <- matrix(c(rep(coeff1, K/9), rep(coeff2, K/9), rep(coeff3, K/6), rep(coeff4, K/6), rep(coeff5, 2*K/9), rep(coeff6, 2*K/9)), nrow = 11)



write.csv(theta_ora, file = paste('Example5_K90_n400_p100_ora_theta_', seed, '.csv', sep = ""))
write.csv(theta_true, file = paste('Example5_K90_n400_p100_ora_thetatrue_', seed, '.csv', sep = ""))



