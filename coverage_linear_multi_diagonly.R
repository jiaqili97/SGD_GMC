data_gen_linear <- function(n,d){
  x <- matrix(rnorm(n*d, mean = 0, sd = 1), ncol=d)
  theta <- seq(0, 1, length.out = d + 1)
  theta <- theta[2:length(theta)]
  epsilon <- rnorm(n, mean = 0, sd = 1)
  y <- x %*% theta + epsilon
  dat <- data.frame(y=y, x=x)
  return(list(dat,theta))
}

loss_grad_linear <- function(x, y, theta){
  loss_grad <- 2*x*as.numeric(t(x) %*% theta - y) # transpose x
}

N_rep <- 200 # number of repetitions
n <- 100000 # number of observations
d <- 2 # number of predictors

gamma_step <- 0.05 # constant learning rate

################### non-overlapping blocks {B_m} ###############
M <- round(sqrt(n))
eta <- rep(0,M)
for (m in 1:M){
  eta[m] <- m^2
}
B <- list()
for (m in 1:(M-1)){
  B[[m]] <- seq.int(eta[m], eta[m+1]-1, by = 1)
}

count <- numeric(n)
true.idx <- matrix(0, nrow = N_rep, ncol = n)
unit_v <- rep(1,d)/sqrt(d) # unit-length vector

for (s in 1:N_rep) {
  
  print(paste(s, "iterations have been processed"))
  
  ############################### generate random samples #############################
  data_linear <- data.matrix(data_gen_linear(n,d)[[1]])
  X <- data_linear[, -1]
  y <- data_linear[, 1]
  theta <- data_gen_linear(n,d)[[2]]
  
  
  ############################### ASGD estimation #############################
  theta_sgd <- matrix(0, nrow = d, ncol = n)
  theta_sgd_bar <- matrix(0, nrow = d, ncol = n)
  
  theta_sgd[,1] <- 0 # can be changed to other initial points
  theta_sgd_bar[,1] <- 0
  
  for (k in 2:n) {
    theta_sgd[,k] <- theta_sgd[,k-1] - gamma_step*loss_grad_linear(X[k-1,], y[k-1], theta_sgd[,k-1])
    theta_sgd_bar[,k] <- ((k-1)*theta_sgd_bar[,k-1] + theta_sgd[,k])/k
  }
  
  # # (uncomment for sanity check) Plot convergence trace for each coordinate of theta
  # matplot(1:n, t(theta_sgd_bar), type = 'l', col = 1:d,
  #         xlab = "n steps",
  #         ylab = "estimated theta",
  #         main = "Convergence of ASGD Iterates")
  # legend("topright", legend = paste("theta", 1:d), col = 1:d, lty = 1)


  ################### Online estimation of long-run covariance #############################
  psi <- numeric(n)
  R_n <- matrix(0, nrow = d, ncol = n)
  delta <- numeric(n)
  K_n <- numeric(n)
  H_n <- matrix(0, nrow = d, ncol = n)
  Vcur_n <- matrix(0, nrow = d, ncol = d)
  V_n <- matrix(0, nrow = d, ncol = d)
  Sigma2_hat <- matrix(0, nrow = d, ncol = d)
  Sigma2_hat_diag <- array(0, dim = c(d,d,n))
  
  psi[1] <- 1
  R_n[,1] <- theta_sgd[,1]
  delta[1] <- 1
  K_n[1] <- 1
  H_n[,1] <- 1
  Vcur_n <- theta_sgd[,1] %o% theta_sgd[,1]
  V_n <- Vcur_n + K_n[1]*theta_sgd_bar[,1]%o%theta_sgd_bar[,1] - 2*H_n[,1]%o%theta_sgd_bar[,1]
  Sigma2_hat <- V_n
  
  for (k in 2:n){
    phi = (psi[k-1] + 1)^2
    
    if (k < phi){
      R_n[,k] <- R_n[,k-1] + theta_sgd[,k]
      psi[k] <- psi[k-1] 
      delta[k] <- delta[k-1] + 1
      K_n[k] <- K_n[k-1] - delta[k-1]^2 + delta[k]^2
      H_n[,k] <- H_n[,k-1] - delta[k-1]*R_n[,k-1] + delta[k]*R_n[,k]
      Vcur_n <- Vcur_n - R_n[,k-1]%o%R_n[,k-1] + R_n[,k]%o%R_n[,k]
    } else {
      R_n[,k] <- theta_sgd[,k]
      delta[k] <- 1
      psi[k] <- psi[k-1] + 1
      K_n[k] <- K_n[k-1] + 1
      H_n[,k] <- H_n[,k-1] + R_n[,k]
      Vcur_n <- Vcur_n + sum(R_n[,k]%o%R_n[,k])
    }
    
    V_n <- Vcur_n + K_n[k]*theta_sgd_bar[,k]%o%theta_sgd_bar[,k] - 2*H_n[,k]%o%theta_sgd_bar[,k]
    Sigma2_hat <- V_n/k
    
    Sigma2_hat_diag[,,k] <- diag(diag(Sigma2_hat)) #only use the estimated diagonals for large d
  }
  
  # Extract diagonals of long-run covariance matrices, i.e., long-run variances
  sigma2_hat <- matrix(0, nrow = d, ncol = n) 
  for (j in 1:d) {
    sigma2_hat[j,] <- Sigma2_hat[j,j]
  }
  # inner_sqrt[s,] <- apply(Sigma2_hat, 3, function(x) unit_v%*%x%*%unit_v)
  inner_sqrt_diag <- apply(Sigma2_hat_diag, 3, function(x) unit_v%*%x%*%unit_v)
  
  # # (uncomment for sanity check) Plot convergence trace for long-run variances of coordinate of theta
  # matplot(1:n, t(sigma2_hat), type = 'l', col = 1:d,
  #         xlab = "n steps",
  #         ylab = "estimated long-run var",
  #         main = "Estimated long-run variances of ASGD")
  # legend("topright", legend = paste("Variance", 1:d), col = 1:d, lty = 1)

  ################### 95% CI with standardized error ###########################################
  CI_left <- unit_v%*%theta_sgd_bar - 1.96*sqrt(abs(inner_sqrt_diag)/n)
  CI_right <- unit_v%*%theta_sgd_bar + 1.96*sqrt(abs(inner_sqrt_diag)/n)
  for (i in 1:n) {
    if (unit_v%*%theta<CI_right[i] & unit_v%*%theta>CI_left[i]) {
      true.idx[s,i] <- 1
      count[i] <- count[i] + 1
    }
  }
}

############################# Plots #############################
coverage_prob <- count/N_rep
coverage_prob[c(n-10000,n-5000,n-1000,n-500,n-100,n)]

coverage_prob[c(5000,10000,50000,80000,90000,n)]

plot(coverage_prob, type = 'l', ylab = 'Coverage Probability', ylim = c(0,1),
     xlab = 'n steps', main = paste('Empirical coverage rate of 95% CI,','d=',d), col = 'black')
abline(h=0.95,col='red', lty=2)

sd.f <- function(x){sd(true.idx[,x])/sqrt(length(true.idx[,x]))}
coverage_prob_sd <- sapply(1:length(true.idx[1,]),sd.f)

coverage_prob_sd[c(n-10000,n-5000,n-1000,n-500,n-100,n)]

coverage_prob_sd[c(5000,10000,50000,80000,90000,n)]


# coverage <- rbind(coverage_prob,coverage_prob_sd)
# 
# coverage[,82600:82700]
# coverage[,c(5000,10000,40000,60000,80000)]






