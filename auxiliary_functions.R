



############### Proximal ADMM ###############


create_adjacency <- function(V,n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0);
  index = t(combn(n,2));
  i <- index[connected_ix,1]
  j <- index[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}


Local_fit <- function(Y_lst, X_lst, family = 'binomial', lambda_lst = NULL){
  
  K <- length(X_lst)
  I_lst <- vector('list', K)
  U_lst <- vector('list', K)
  beta_lst <- c()
  
  if (family == "binomial") {
    for (k in 1:K){
      Y <- Y_lst[[k]]
      X <- X_lst[[k]]
      
      if (length(lambda_lst) == 1){
        lambda.cv <- lambda_lst[1]
      }
      else{
        cv.result <- cv.glmnet(X, Y, family = family, lambda = lambda_lst)
        lambda.cv <- cv.result$lambda.min
      }
      
      model <- glmnet(X, Y, family = family, lambda = lambda.cv)
      beta_fit <- c(as.vector(model$a0), as.vector(model$beta))
      
      n <- length(Y)
      X_all <- cbind(rep(1, n), X)
      
      pi_vec <- as.vector(1 / (1 + exp(- X_all %*% beta_fit)))
      grad <- t(X_all) %*% (pi_vec - Y)
      I_mat <- t(X_all) %*% diag(pi_vec * (1- pi_vec)) %*% X_all
      U_mat <- I_mat %*% beta_fit - grad
      I_lst[[k]] <- I_mat
      U_lst[[k]] <- U_mat
      beta_lst <- cbind(beta_lst, beta_fit)
    }
  }
  else if (family == "gaussian"){
    for (k in 1:K){
      Y <- Y_lst[[k]]
      X <- X_lst[[k]]
      
      if (length(lambda_lst) == 1){
        lambda.cv <- lambda_lst[1]
      }
      else{
        cv.result <- cv.glmnet(X, Y, family = family, lambda = lambda_lst)
        lambda.cv <- cv.result$lambda.min
      }
      
      model <- glmnet(X, Y, family = family, lambda = lambda.cv)
      beta_fit <- c(as.vector(model$a0), as.vector(model$beta))
      
      n <- length(Y)
      X_all <- cbind(rep(1, n), X)
      grad <- t(X_all) %*% (as.numeric(X_all %*% beta_fit) - Y)
      I_mat <- t(X_all) %*% X_all
      U_mat <- I_mat %*% beta_fit - grad
      
      I_lst[[k]] <- I_mat
      U_lst[[k]] <- U_mat
      beta_lst <- cbind(beta_lst, beta_fit)
    }
  }
  
  
  return(list(I = I_lst, U = U_lst, beta = beta_lst))
}

debias_fit <- function(H_lst, g_lst, beta_fit, n_lst, tau_lst = NULL){
  
  
  M <- length(H_lst)
  n <- sum(n_lst) / M
  p <- length(H_lst[[1]][1, ]) - 1
  
  if (is.null(tau_lst)){
    tau_lst <- sqrt(log(p) / n)
  }
  
  
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  
  for (i in 1:M){
    H <- H_lst[[i]]
    d <- g_lst[[i]]
    mat_all <- cbind(H, d)
    mat_all <- rbind(mat_all, t(c(d, max(mat_all) + n)))
    
    svd_result <- svd(mat_all)
    s_value <- svd_result$d
    s_mat <- diag(sqrt(s_value))[ ,1:(min(p + 1, n_lst[[i]]) + 1)]
    data_all <- svd_result$u %*% s_mat
    
    X <- t(data_all[-length(data_all[ ,1]),])
    Y <- data_all[length(data_all[ ,1]),]
    X_lst[[i]] <- X
    Y_lst[[i]] <- Y
  }
  
  X_all <- c()
  for (i in 1:M){
    X_all <- rbind(X_all, X_lst[[i]])
  }
  
  X_alpha <- X_lst[[1]]
  for (i in 2:M){
    X_alpha <- bdiag(X_alpha, X_lst[[i]])
  }
  X_all <- cbind(X_all, X_alpha)
  
  enlarge_mat <- matrix(0, p + 1, p + 1)
  for (i in 1:M) {
    enlarge_mat <- cbind(enlarge_mat, 10000 * diag(rep(1, p + 1)))
  }
  X_enlarge <- rbind(X_all, enlarge_mat)
  
  
  beta_deb_mat <- matrix(0, p + 1, M)
  for (j in 1:(p + 1)) {
    ej <- rep(0, p + 1)
    ej[j] <- 1
    
    for (m in 1:M){
      X_mj <- X_lst[[m]][, -j]
      X_j <- X_lst[[m]][,j]
      
      if (length(tau_lst) >= 2){
        model_fit <- cv.glmnet(X_mj, X_j, lambda = tau_lst * sqrt(n / n_lst[m]), intercept = F,
                               standardize = F)
        tau <- model_fit$lambda.min
        #print(tau)
        model_fit <- glmnet(X_mj, X_j, lambda = tau, intercept = F, standardize = F)
      }else{
        tau <- tau_lst[1]
        model_fit <- glmnet(X_mj, X_j, lambda = tau, intercept = F, standardize = F)
      }
      
      u_j <- rep(1, p + 1)
      u_j[-j] <- - model_fit$beta
      sigma2 <- mean((X_j - X_mj %*% model_fit$beta)^2) + tau * sqrt(n / n_lst[m]) * sum(abs(model_fit$beta))
      u_m <- as.vector(u_j / sigma2)
      
      beta_deb <- beta_fit[j, m] + 
        t(u_m) %*% (g_lst[[m]] - H_lst[[m]] %*% beta_fit[,m]) / n_lst[m]
      beta_deb_mat[j, m] <- beta_deb
    }
  }
  
  return(beta_deb_mat)
}



BIC_cal <- function(theta_mat, V_tilde_lst, zeta_tilde_lst, N, con, K_hat){
  
  p <- nrow(theta_mat)
  K <- ncol(theta_mat)
  
  loss <- rep(0, K)
  for (k in 1:K) {
    loss[k] <- as.numeric(t(theta_mat[,k]) %*% V_tilde_lst[[k]] %*% theta_mat[,k] - 2*t(theta_mat[,k]) %*% zeta_tilde_lst[[k]])
  }
  
  #print(loss)
  loss_val <- sum(loss)/N
  qhat <- sum(theta_mat[,1] != 0)*K_hat
  penal <- con*log(log(K*p))*log(N)*(qhat)/N
  #penal <- 2*(qhat)/N
  BIC_val <- loss_val + penal
  return(BIC_val)
}


tuning_selection_cpp <- function(theta_ini, V_tilde_lst, zeta_tilde_lst, K, p, N, n_vec,
                                lambda1_candi, lambda2_candi, con){
  
  index <- t(combn(K,2))
  m <- nrow(index)
  V_tilde_mat <- do.call("cbind", V_tilde_lst)
  zeta_tilde_mat <- do.call("cbind", zeta_tilde_lst)
  
  nlambda1 <- length(lambda1_candi)
  nlambda2 <- length(lambda2_candi)
  bic_mat <- matrix(0, nrow = nlambda1, ncol = nlambda2)
  iter <- 0
  for (l1 in 1:nlambda1) {
    for (l2 in 1:nlambda2) {
      result <- ICR_pg_cpp(theta_ini, V_tilde_mat, zeta_tilde_mat, n_vec, index,
                           nu = 1, tau = 3, lambda1 = lambda1_candi[l1], lambda2 = lambda2_candi[l2],  
                           eps_outer = 1e-3, max_iter = 100)
      ## take some average and sparse procedures to get the final estimators theta_average
      theta_average <- matrix(0, nrow = p, ncol = K)
      
      
      Ad_final <- create_adjacency(result$alpha, K)
      G_final <- graph.adjacency(Ad_final, mode = 'upper')
      cls_final_fea <- components(G_final);
      #print(theta_average)
      for (g in 1:cls_final_fea$no) {
        theta_average[,cls_final_fea$membership == g] <- apply(as.matrix(result$theta[,cls_final_fea$membership == g]), MARGIN = 1, mean)
      }
      
      
      bic_mat[l1,l2] <- BIC_cal(theta_average, V_tilde_lst, zeta_tilde_lst, N, con = con, K_hat = cls_final_fea$no)
      iter <- iter + 1
      print(iter)
    }
  }
  
  print(bic_mat)
  id <- which(bic_mat == min(bic_mat), arr.ind = T)[1,]
  lambda1_final <- lambda1_candi[id[1]]
  lambda2_final <- lambda2_candi[id[2]]
  lambda_final <- c(lambda1_final, lambda2_final)
  return(lambda_final)
}





############### Proximal ADMM for IP ###############

BIC_logistic_cal <- function(theta_mat, X_intercept_lst, Y_lst, n_vec, con, K_hat){
  
  p <- nrow(theta_mat)
  K <- ncol(theta_mat)
  N <- sum(n_vec)
  

  loss <- rep(0, K)
  for (k in 1:K) {
    loss[k] <- sum(as.numeric(X_intercept_lst[[k]] %*% theta_mat[,k]) * Y_lst[[k]]) - sum(log(1+exp(X_intercept_lst[[k]] %*% theta_mat[,k])))
  }
  
  #print(loss)
  loss_val <- -2 * sum(loss)/N
  qhat <- sum(theta_mat[,1] != 0)*K_hat
  penal <- con*log(log(K*p))*log(N)*(qhat)/N
  #penal <- 2*(qhat)/N
  BIC_val <- loss_val + penal
  return(BIC_val)
}


tuning_selection_logistic_cpp <- function(theta_ini, X_intercept_lst, Y_lst, K, p, N, n_vec, lambda1_candi, lambda2_candi, con){
  
  index <- t(combn(K,2))
  m <- nrow(index)
  X_intercept_mat <- do.call("rbind", X_intercept_lst)
  Y_vec <- do.call("c", Y_lst)
  
  nlambda1 <- length(lambda1_candi)
  nlambda2 <- length(lambda2_candi)
  bic_mat <- matrix(0, nrow = nlambda1, ncol = nlambda2)
  iter <- 0
  for (l1 in 1:nlambda1) {
    for (l2 in 1:nlambda2) {
      result <- ICR_pg_logistic_cpp(theta_ini, X_intercept_mat, Y_vec, n_vec, index, nu = 1, tau = 3, 
                                    lambda1 = lambda1_candi[l1], lambda2 = lambda2_candi[l2],  eps_outer = 1e-5, max_iter = 500)
      ## take some average and sparse procedures to get the final estimators theta_average
      theta_average <- matrix(0, nrow = p, ncol = K)
      
      
      Ad_final <- create_adjacency(result$alpha, K)
      G_final <- graph.adjacency(Ad_final, mode = 'upper')
      cls_final_fea <- components(G_final);
      #print(theta_average)
      for (g in 1:cls_final_fea$no) {
        theta_average[,cls_final_fea$membership == g] <- apply(as.matrix(result$theta[,cls_final_fea$membership == g]), MARGIN = 1, mean)
      }
      
      
      bic_mat[l1,l2] <- BIC_logistic_cal(theta_average, X_intercept_lst, Y_lst, n_vec, con = con, K_hat = cls_final_fea$no)
      iter <- iter + 1
      print(iter)
    }
  }
  
  print(bic_mat)
  id <- which(bic_mat == min(bic_mat), arr.ind = T)[1,]
  lambda1_final <- lambda1_candi[id[1]]
  lambda2_final <- lambda2_candi[id[2]]
  lambda_final <- c(lambda1_final, lambda2_final)
  return(lambda_final)
}



BIC_linear_cal <- function(theta_mat, X_intercept_lst, Y_lst, n_vec, con, K_hat){
  
  p <- nrow(theta_mat)
  K <- ncol(theta_mat)
  N <- sum(n_vec)
  
  
  loss <- rep(0, K)
  for (k in 1:K) {
    loss[k] <- sum((as.numeric(X_intercept_lst[[k]] %*% theta_mat[,k]) - Y_lst[[k]])^2)
  }
  
  #print(loss)
  #loss_val <- log(sum(loss)/N)
  loss_val <- sum(loss)/N
  qhat <- sum(theta_mat[,1] != 0)*K_hat
  penal <- con*log(log(K*p))*log(N)*(qhat)/N
  #penal <- 2*(qhat)/N
  BIC_val <- loss_val + penal
  return(BIC_val)
}




tuning_selection_linear_cpp <- function(theta_ini, X_intercept_lst, Y_lst, K, p, N, n_vec, lambda1_candi, lambda2_candi, con){
  
  index <- t(combn(K,2))
  m <- nrow(index)
  X_intercept_mat <- do.call("rbind", X_intercept_lst)
  Y_vec <- do.call("c", Y_lst)
  
  nlambda1 <- length(lambda1_candi)
  nlambda2 <- length(lambda2_candi)
  bic_mat <- matrix(0, nrow = nlambda1, ncol = nlambda2)
  iter <- 0
  for (l1 in 1:nlambda1) {
    for (l2 in 1:nlambda2) {
      result <- ICR_pg_linear_cpp(theta_ini, X_intercept_mat, Y_vec, n_vec, index, nu = 1, tau = 3, 
                                    lambda1 = lambda1_candi[l1], lambda2 = lambda2_candi[l2],  eps_outer = 1e-5, max_iter = 500)
      ## take some average and sparse procedures to get the final estimators theta_average
      theta_average <- matrix(0, nrow = p, ncol = K)
      
      
      Ad_final <- create_adjacency(result$alpha, K)
      G_final <- graph.adjacency(Ad_final, mode = 'upper')
      cls_final_fea <- components(G_final);
      #print(theta_average)
      for (g in 1:cls_final_fea$no) {
        theta_average[,cls_final_fea$membership == g] <- apply(as.matrix(result$theta[,cls_final_fea$membership == g]), MARGIN = 1, mean)
      }
      
      
      bic_mat[l1,l2] <- BIC_linear_cal(theta_average, X_intercept_lst, Y_lst, n_vec, con = con, K_hat = cls_final_fea$no)
      iter <- iter + 1
      print(iter)
    }
  }
  
  print(bic_mat)
  id <- which(bic_mat == min(bic_mat), arr.ind = T)[1,]
  lambda1_final <- lambda1_candi[id[1]]
  lambda2_final <- lambda2_candi[id[2]]
  lambda_final <- c(lambda1_final, lambda2_final)
  return(lambda_final)
}


