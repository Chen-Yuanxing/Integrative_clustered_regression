
#include <RcppArmadillo.h>
#include <string.h>
#include <stdio.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// calculate prox_L2

arma::vec prox_L2(arma::vec x, double sigma){
  int n = x.n_elem;
  double lv;
  arma::vec px(n);
  
  lv = norm(x, 2);
  
  if (lv == 0.)
    px = x;
  else
    px = fmax(0., 1. - (sigma/lv))*x;
  
  return px;
}

arma::vec prox_mcp(arma::vec x, double eta, double tau, double lambda){
  int n = x.n_elem;
  double lv;
  arma::vec px(n);
  
  lv = norm(x, 2);
  
  if (lv <= lambda*eta)
    px = arma::zeros<arma::vec>(n);
  else if (lv <= lambda*tau)
    px = (1 - lambda*eta/lv)/(1 - eta/tau) * x;
  else
    px = x;
  
  return px;
}

double mcp_dx(double x, double tau, double lambda){

  double deriv;

  if (x <= lambda*tau)
    deriv = lambda*(1.-x/(lambda*tau));
  else
    deriv = 0.;
  
  return deriv;
}

arma::vec proj_L2(arma::vec x, double sigma){
  int n = x.n_elem;
  double lv;
  arma::vec px(n);
  
  lv = norm(x, 2);
  
  if (lv > sigma)
    px = x * sigma/lv;
  else
    px = x;
  
  return px;
}


arma::mat update_theta_cpp(arma::mat theta, arma::mat Lambda, arma::mat V_tilde_mat, arma::mat zeta_tilde_mat, 
                       arma::mat A, arma::mat index, arma::vec n_vec, double eta, double nu, double tau, 
                       double lambda1, double lambda2, int iter_outer, int max_iter_inner){
  
  int p = theta.n_rows;
  int K = theta.n_cols;
  int m = index.n_rows;
  
  int N = 0;
  for (int l = 0; l < K; l++) {
    N = N + n_vec(l);
  }
  
  arma::mat theta_middle  = theta;
  arma::mat theta_old  = theta;
  arma::mat theta_new  = theta;
  arma::mat f_gra(p, K);
  arma::mat prox_term1(p, m);
  
  arma::vec omega(m);
 
  
  double res_inner;
  double eps_inner = 0.01/(nu*iter_outer);
  double xi_old = 1.;
  double xi_new = 1.;
  

  for (int iter_inner = 0; iter_inner < max_iter_inner; iter_inner++) {
    
    for (int l=0; l<m; l++) {
      int i = index(l, 0);
      int j = index(l, 1);
      
      double l2_norm = norm(theta_middle.col(i-1) - theta_middle.col(j-1), 2);
      omega(l) = mcp_dx(l2_norm, lambda2, tau);
      
    }
    
    for (int k=0; k<K; k++) {
      int r1 = k*p;
      int r2 = (k+1)*p - 1;
      
      f_gra.col(k) = 2.0/N * (V_tilde_mat.cols(r1, r2) * theta_middle.col(k) - zeta_tilde_mat.col(k));

    }
    
    arma::mat theta_middle_A(p, m);
    for (int l=0; l<m; l++) {
      int i = index(l, 0);
      int j = index(l, 1);
      
      theta_middle_A.col(l) = theta_middle.col(i-1) - theta_middle.col(j-1);
      
    }
    
    
    arma::mat term1 = nu * theta_middle_A + Lambda;
    for (int l=0; l<m; l++) {
      prox_term1.col(l) = proj_L2(term1.col(l), omega(l));
      
    }
    
    arma::mat prox_term1_At(p, K);
    for(int k=0;k<K;k++){
      arma::vec v1 = arma::zeros<arma::vec>(p);
      arma::vec v2 = arma::zeros<arma::vec>(p);
      for(int ix=0;ix<m;ix++){
        int l1 = index(ix,0);
        if(l1 == k + 1)
          v1 = v1 + prox_term1.col(ix);
        else
          v1 = v1 + arma::zeros<arma::vec>(p);
        
        int l2 = index(ix,1);
        if(l2 == k + 1)
          v2 = v2 + prox_term1.col(ix);
        else
          v2 = v2 + arma::zeros<arma::vec>(p);
      }
      prox_term1_At.col(k) = v1 - v2;
      
    }
    
    arma::mat phi_gra = f_gra + prox_term1_At;
    
    arma::mat uu1 = arma::trans(theta_middle);
    arma::mat uu2 = arma::trans(phi_gra);
    arma::mat uu3(K, p);
    
    uu3.col(0) = uu1.col(0) - eta*uu2.col(0);
    for (int j=1; j<p; j++) {
      uu3.col(j) = prox_mcp(uu1.col(j) - eta*uu2.col(j), eta, tau, lambda1);

    }
    theta_new = arma::trans(uu3);
    xi_new = (1 + sqrt(1+4*pow(xi_old, 2)))/2;
    theta_middle = theta_new - (xi_old-1)/xi_new * (theta_new - theta_old);
    
    res_inner = arma::norm(theta_new - theta_old, "fro");
    if (res_inner <= eps_inner) break;
    
    theta_old = theta_new;
    xi_old = xi_new;
    
  }

  return theta_new;
  
}


arma::mat update_Lambda_cpp(arma::mat theta, arma::mat Lambda_old, arma::mat A, arma::mat index,
                        double nu, double tau, double lambda2){
  
  int p = theta.n_rows;
  int m = index.n_rows;
  

  arma::mat Lambda_new(p, m);
  arma::mat theta_A(p, m);
  arma::vec omega(m);
  
  for (int l=0; l<m; l++) {
    int i = index(l, 0);
    int j = index(l, 1);
    
    double l2_norm = norm(theta.col(i-1) - theta.col(j-1), 2);
    omega(l) = mcp_dx(l2_norm, lambda2, tau);
    theta_A.col(l) = theta.col(i-1) - theta.col(j-1);
  }
  

  arma::mat term1 = nu * theta_A + Lambda_old;

  for (int l=0; l<m; l++) {
    Lambda_new.col(l) = proj_L2(term1.col(l), omega(l));
    
  }
  return Lambda_new;
  
}


arma::mat update_alpha_cpp(arma::mat theta, arma::mat Lambda, arma::mat A, arma::mat index,
                        double nu, double tau, double lambda2){
  
  int p = theta.n_rows;
  int m = index.n_rows;
  
  arma::mat alpha(p, m);
  arma::mat theta_A(p, m);
  arma::vec omega(m);
  
  for (int l=0; l<m; l++) {
    int i = index(l, 0);
    int j = index(l, 1);
    
    double l2_norm = norm(theta.col(i-1) - theta.col(j-1), 2);
    omega(l) = mcp_dx(l2_norm, lambda2, tau);
    theta_A.col(l) = theta.col(i-1) - theta.col(j-1);
  }
  
  arma::mat term2 = theta_A + Lambda/nu;
  
  for (int l=0; l<m; l++) {
    alpha.col(l) = prox_L2(term2.col(l), omega(l)/nu);
    
  }
  return alpha;
  
}


// [[Rcpp::export]]
List ICR_pg_cpp(arma::mat theta_ini, arma::mat V_tilde_mat, arma::mat zeta_tilde_mat, arma::vec n_vec, 
            arma::mat index, double nu, double tau, double lambda1, double lambda2, 
            double eps_outer, int max_iter) {
  
  int p = theta_ini.n_rows;
  int K = theta_ini.n_cols;
  int m = index.n_rows;
  
  arma::mat theta_old = theta_ini;
  arma::mat Lambda_old = arma::zeros<arma::mat>(p, m);
  arma::mat A = arma::zeros<arma::mat>(K, m);
  arma::mat theta_new(p, K);
  arma::mat Lambda_new(p, m);
  arma::mat alpha(p, m);
  
  for (int l=0; l<m; l++) {
    int i = index(l, 0);
    int j = index(l, 1);
    
    A(i-1, l) = 1;
    A(j-1, l) = -1;
    
  }
  
  arma::mat G = A * arma::trans(A);
  double eta = 1.0/(1 + 2*nu*max(G.diag(0)));
  
  int its = 1;

  for (int iter = 0; iter < max_iter; iter++) {
    theta_new = update_theta_cpp(theta_old, Lambda_old, V_tilde_mat, zeta_tilde_mat, A, index, n_vec, eta,
                             nu, tau, lambda1, lambda2, its, 100);
    Lambda_new = update_Lambda_cpp(theta_new, Lambda_old, A, index, nu, tau, lambda2);
    
    double res_outer = fmax(arma::norm(Lambda_new - Lambda_old, "fro"), arma::norm(theta_new - theta_old, "fro"));
    its = its + 1;
    
    if (res_outer <= eps_outer) break;
      
    theta_old = theta_new;
    Lambda_old = Lambda_new;

  }

  alpha = update_alpha_cpp(theta_new, Lambda_new, A, index, nu, tau, lambda2);
  
  return  List::create(Named("theta") = theta_new,
                       Named("alpha") = alpha,
                       Named("Lambda") = Lambda_new,
                       Named("iterations") = its);
  
}




// [[Rcpp::export]]
List ICR_pg_logistic_cpp(arma::mat theta_ini, arma::mat X_all, arma::vec Y_all, arma::vec n_vec, 
                         arma::mat index, double nu, double tau, double lambda1, double lambda2, 
                         double eps_outer, int max_iter) {
  
  int p = theta_ini.n_rows;
  int K = theta_ini.n_cols;
  int m = index.n_rows;
  int N = arma::sum(n_vec);
  

  
  arma::mat theta_old = theta_ini;
  arma::mat Lambda_old = arma::zeros<arma::mat>(p, m);
  arma::mat A = arma::zeros<arma::mat>(K, m);
  arma::mat theta_new(p, K);
  arma::mat Lambda_new(p, m);
  arma::mat alpha(p, m);
  
  for (int l=0; l<m; l++) {
    int i = index(l, 0);
    int j = index(l, 1);
    
    A(i-1, l) = 1;
    A(j-1, l) = -1;
  }
  
  arma::mat G = A * arma::trans(A);
  
  
  int its = 1;
  
  arma::mat hessian_mat(p, p*K);
  arma::mat zeta_mat(p, K); 
  
  for (int iter = 0; iter < max_iter; iter++) {
    
    int n_1 = 0;
    arma::vec eigen_vec(K);
    
    for (int k=0; k<K; k++) {
      n_1 = n_1 + n_vec(k);
      int r1 = k*p;
      int r2 = (k+1)*p - 1;
      
      arma::vec pi_vec = 1/(1 + exp(-1 * X_all.rows(n_1-n_vec(k), n_1-1) * theta_old.col(k)));
      arma::mat diag_mat = arma::diagmat(pi_vec % (1-pi_vec)); 
      hessian_mat.cols(r1, r2) = arma::trans(X_all.rows(n_1-n_vec(k), n_1-1)) * diag_mat * X_all.rows(n_1-n_vec(k), n_1-1);

      arma::vec grad = arma::trans(X_all.rows(n_1-n_vec(k), n_1-1)) * (pi_vec - Y_all.subvec(n_1-n_vec(k), n_1-1));
      zeta_mat.col(k) = hessian_mat.cols(r1, r2) * theta_old.col(k) - grad;
      
      arma::vec eigenvalues = arma::eig_sym(hessian_mat.cols(r1, r2));
      eigen_vec(k) = arma::max(eigenvalues);
    }
    
    double eta = 1/(2*arma::max(eigen_vec)/N + 2*nu*max(G.diag(0)));
    
    theta_new = update_theta_cpp(theta_old, Lambda_old, hessian_mat, zeta_mat, A, index, n_vec, eta,
                                 nu, tau, lambda1, lambda2, its, 100);
    
   // theta_new = update_theta_logistic_cpp(X_all, Y_all, theta_old, Lambda_old, A, index, n_vec, eta,
   //                                       nu, tau, lambda1, lambda2, its, 100);
    
    Lambda_new = update_Lambda_cpp(theta_new, Lambda_old, A, index, nu, tau, lambda2);
    
    double res_outer = fmax(arma::norm(Lambda_new - Lambda_old, "fro"), arma::norm(theta_new - theta_old, "fro"));
    its = its + 1;
    
    if (res_outer <= eps_outer) break;
    
    theta_old = theta_new;
    Lambda_old = Lambda_new;
    
  }
  
  alpha = update_alpha_cpp(theta_new, Lambda_new, A, index, nu, tau, lambda2);
  
  return  List::create(Named("theta") = theta_new,
                       Named("alpha") = alpha,
                       Named("Lambda") = Lambda_new,
                       Named("iterations") = its);
  
}


// [[Rcpp::export]]
List ICR_pg_linear_cpp(arma::mat theta_ini, arma::mat X_all, arma::vec Y_all, arma::vec n_vec, 
                         arma::mat index, double nu, double tau, double lambda1, double lambda2, 
                         double eps_outer, int max_iter) {
  
  int p = theta_ini.n_rows;
  int K = theta_ini.n_cols;
  int m = index.n_rows;
  int N = arma::sum(n_vec);
  
  
  arma::mat theta_old = theta_ini;
  arma::mat Lambda_old = arma::zeros<arma::mat>(p, m);
  arma::mat A = arma::zeros<arma::mat>(K, m);
  arma::mat theta_new(p, K);
  arma::mat Lambda_new(p, m);
  arma::mat alpha(p, m);
  
  for (int l=0; l<m; l++) {
    int i = index(l, 0);
    int j = index(l, 1);
    
    A(i-1, l) = 1;
    A(j-1, l) = -1;
  }
  
  arma::mat G = A * arma::trans(A);
  
  
  int its = 1;
  
  arma::mat hessian_mat(p, p*K);
  arma::mat zeta_mat(p, K); 
  
  for (int iter = 0; iter < max_iter; iter++) {
    
    int n_1 = 0;
    arma::vec eigen_vec(K);
    
    for (int k=0; k<K; k++) {
      n_1 = n_1 + n_vec(k);
      int r1 = k*p;
      int r2 = (k+1)*p - 1;
      
      //arma::vec pi_vec = 1/(1 + exp(-1 * X_all.rows(n_1-n_vec(k), n_1-1) * theta_old.col(k)));
      //arma::mat diag_mat = arma::diagmat(pi_vec % (1-pi_vec)); 
      hessian_mat.cols(r1, r2) = arma::trans(X_all.rows(n_1-n_vec(k), n_1-1)) * X_all.rows(n_1-n_vec(k), n_1-1);
      
      arma::vec grad = arma::trans(X_all.rows(n_1-n_vec(k), n_1-1)) * (X_all.rows(n_1-n_vec(k), n_1-1) * theta_old.col(k) - Y_all.subvec(n_1-n_vec(k), n_1-1));
      zeta_mat.col(k) = hessian_mat.cols(r1, r2) * theta_old.col(k) - grad;
      
      arma::vec eigenvalues = arma::eig_sym(hessian_mat.cols(r1, r2));
      eigen_vec(k) = arma::max(eigenvalues);
    }
    
    double eta = 1/(2*arma::max(eigen_vec)/N + 2*nu*max(G.diag(0)));
    
    theta_new = update_theta_cpp(theta_old, Lambda_old, hessian_mat, zeta_mat, A, index, n_vec, eta,
                                 nu, tau, lambda1, lambda2, its, 100);
    
    // theta_new = update_theta_logistic_cpp(X_all, Y_all, theta_old, Lambda_old, A, index, n_vec, eta,
    //                                       nu, tau, lambda1, lambda2, its, 100);
    
    Lambda_new = update_Lambda_cpp(theta_new, Lambda_old, A, index, nu, tau, lambda2);
    
    double res_outer = fmax(arma::norm(Lambda_new - Lambda_old, "fro"), arma::norm(theta_new - theta_old, "fro"));
    its = its + 1;
    
    if (res_outer <= eps_outer) break;
    
    theta_old = theta_new;
    Lambda_old = Lambda_new;
    
  }
  
  alpha = update_alpha_cpp(theta_new, Lambda_new, A, index, nu, tau, lambda2);
  
  return  List::create(Named("theta") = theta_new,
                       Named("alpha") = alpha,
                       Named("Lambda") = Lambda_new,
                       Named("iterations") = its);
  
}



