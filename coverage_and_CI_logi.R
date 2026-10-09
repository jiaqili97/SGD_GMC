data_gen_logi <- function(n,d){
  x <- matrix(rnorm(n*d, mean = 0, sd = 1), ncol=d)
  theta <- seq(0, 1, length.out = d + 1)
  theta <- theta[2:length(theta)]
  p <- 1/(1 + exp(-x %*% theta))
  y <- rbinom(n,1,p)
  dat <- data.frame(y=y, x=x)
  return(list(dat,theta))
}

loss_grad_logi <- function(x, y, theta){
  loss_grad <- -(y-1/(1+exp(-x * theta)))* x
  # loss_grad <- (1-y-1/(1+exp(-x * theta)))* x
}

sgd_est <- function(x, y, gamma_step){
  theta_1 <- 0 
  theta_2 <- 0
  
  theta_cur <- numeric(n)
  theta_hat <- numeric(n)
  
  for (k in 2:n) {
    theta_1 <- theta_1 - gamma_step*loss_grad_logi(x[k-1], y[k-1], theta_1)
    theta_2 <- theta_2 - 2*gamma_step*loss_grad_logi(x[k-1], y[k-1], theta_2)
    theta_cur[k] <- 2*theta_1 - theta_2
    theta_hat[k] <- ((k-1)*theta_hat[k-1] + theta_cur[k])/k
  }
  return(list(theta_cur,theta_hat))
}

longrun_est <- function(n,theta_cur,theta_hat){
  M <- round(sqrt(n))
  eta <- rep(0,M)
  for (m in 1:M){
    eta[m] <- m^2
  }
  B <- list()
  for (m in 1:(M-1)){
    B[[m]] <- seq.int(eta[m], eta[m+1]-1, by = 1)
  }
  
  psi <- numeric(n)
  R_n <- numeric(n)
  delta <- numeric(n)
  K_n <- numeric(n)
  H_n <- numeric(n)
  Vcur_n <- numeric(n)
  V_n <- numeric(n)
  sigma2_hat <- numeric(n)
  
  psi[1] <- 1
  
  R_n[1] <- theta_cur[1]
  delta[1] <- 1
  K_n[1] <- 1
  H_n[1] <- 1
  Vcur_n[1] <- theta_cur[1]^2
  V_n[1] <- Vcur_n[1] + K_n[1]*theta_hat[1]^2 - 2*H_n[1]*theta_hat[1]
  sigma2_hat[1] <- V_n[1]/1
  
  for (k in 2:n){
    phi = (psi[k-1] + 1)^2
    
    if (k < phi){
      R_n[k] <- R_n[k-1] + theta_cur[k]
      psi[k] <- psi[k-1] 
      delta[k] <- delta[k-1] + 1
      K_n[k] <- K_n[k-1] - delta[k-1]^2 + delta[k]^2
      H_n[k] <- H_n[k-1] - delta[k-1]*R_n[k-1] + delta[k]*R_n[k]
      Vcur_n[k] <- Vcur_n[k-1] - R_n[k-1]^2 + R_n[k]^2
    } else {
      R_n[k] <- theta_cur[k]
      delta[k] <- 1
      psi[k] <- psi[k-1] + 1
      K_n[k] <- K_n[k-1] + 1
      H_n[k] <- H_n[k-1] + R_n[k]
      Vcur_n[k] <- Vcur_n[k-1] + sum(R_n[k]^2)
    }
    
    V_n[k] <- Vcur_n[k] + K_n[k]*theta_hat[k]^2 - 2*H_n[k]*theta_hat[k]
    sigma2_hat[k] <- V_n[k]/k
  }
  return(sigma2_hat)
}

longrun_est_new <- function(n,theta_cur,theta_hat){
  
  Theta.hat = theta_hat[n]                                                                                                                        
  
  N = round(sqrt(n))                                                                                                                                 
  Vn = 0
  count.tau = 0
  for(j in 1:(N-1)){
    tmp.Phi = j^2
    tmp.Tau = max(1,round(2*log(j)))                                                             ## tau sequence ##
    for(o in (j^2):((j+1)^2-1)){
      tmp.Idx = seq(j^2,o)
      tmp.Sum = sum(theta_cur[tmp.Idx]-Theta.hat)
      tmp.I = 1*(o >= tmp.Tau + j^2)
      Vn = Vn+(2*tmp.Sum*(theta_cur[o]-Theta.hat)-(theta_cur[o]-Theta.hat)^2)*tmp.I
      count.tau = count.tau + tmp.I
    }
  }
  for(oo in (N^2):n){
    tmp.Idx = seq(N^2,oo)
    tmp.Tau = max(1,round(2*log(N)))
    tmp.Sum = sum(theta_cur[tmp.Idx]-Theta.hat)
    tmp.I = 1*(oo >= tmp.Tau + N^2)
    Vn = Vn+(2*tmp.Sum*(theta_cur[oo]-Theta.hat)-(theta_cur[oo]-Theta.hat)^2)*tmp.I
    count.tau = count.tau + tmp.I
  }
  Vn/count.tau
  
  ##############################################################                              ## online bias reduction ##
  theta.ob = numeric(n)
  psi.ob = numeric(n)
  calR.ob = numeric(n)
  calU.ob = numeric(n)
  calK.ob = numeric(n)
  calH.ob = numeric(n)
  count.ob = numeric(n)
  phi.ob = numeric(n)
  delta.ob = numeric(n)
  
  Vn = numeric(n)
  sigma2_hat_new = numeric(n)
  
  theta.ob[1] = theta_cur[1]
  psi.ob[1] = 1
  phi.ob[1] = 1
  calR.ob[1] = theta_cur[1]
  delta.ob[1] = 1
  calU.ob[1] = 0
  calK.ob[1] = 0 
  calH.ob[1] = 0 
  count.ob = 0
  
  Vn[1] = 0
  sigma2_hat_new[1] = 0
  
  for(m in 2:n){
    tmp.phi = (psi.ob[m-1] + 1)^2
    theta.ob[m] = ((m-1)*theta.ob[m-1]+theta_cur[m])/(m)
    tmp.Tau.ob = max(1,round(2*log(psi.ob[m-1])))
    Tau.Sum.Idx = seq(phi.ob[m-1],phi.ob[m-1]+tmp.Tau.ob-1)
    Stau.ob = sum(theta_cur[Tau.Sum.Idx])                                                    ## S_tau ##
    if(m < tmp.phi){
      psi.ob[m] = psi.ob[m-1]
      phi.ob[m] = phi.ob[m -1]
      delta.ob[m] = delta.ob[m-1]+1
      calR.ob[m] = calR.ob[m-1]+theta_cur[m]
      if(m >= phi.ob[m] + tmp.Tau.ob && m - 1 >= phi.ob[m-1] + tmp.Tau.ob){
        count.ob = count.ob + 1
        calU.ob[m] = calU.ob[m-1]+theta_cur[m]^2+2*theta_cur[m]*calR.ob[m-1]
        calK.ob[m] = calK.ob[m-1]-delta.ob[m-1]^2+delta.ob[m]^2
        calH.ob[m] = calH.ob[m-1]-delta.ob[m-1]*calR.ob[m-1]+delta.ob[m]*calR.ob[m]
      }
      if(m >= phi.ob[m] + tmp.Tau.ob && m - 1 < phi.ob[m-1] + tmp.Tau.ob){
        count.ob = count.ob + 1
        calU.ob[m] = calU.ob[m-1]+theta_cur[m]^2+2*theta_cur[m]*calR.ob[m-1]
        calK.ob[m] = calK.ob[m-1]+delta.ob[m]^2-(max(1,round(2*log(psi.ob[m]))))^2
        calH.ob[m] = calH.ob[m-1]+delta.ob[m]*calR.ob[m]-(max(1,round(2*log(psi.ob[m]))))*Stau.ob
      }
      if(m < phi.ob[m] + tmp.Tau.ob){
        calU.ob[m] = calU.ob[m-1]
        calK.ob[m] = calK.ob[m-1]
        calH.ob[m] = calH.ob[m-1]
      }
    }
    if(m == tmp.phi){
      psi.ob[m] = psi.ob[m-1]+1
      phi.ob[m] = m
      delta.ob[m] = 1
      calR.ob[m] = theta_cur[m]
      calU.ob[m] = calU.ob[m-1]
      calK.ob[m] = calK.ob[m-1]
      calH.ob[m] = calH.ob[m-1]
    }
    Vn[m] = calU.ob[m]+calK.ob[m]*theta.ob[m]^2-2*calH.ob[m]*theta.ob[m]
    sigma2_hat_new[m] = Vn[m]/count.ob
  }
  return(sigma2_hat_new)
}



############################ Sanity Check (N_rep=1) ##################################

n <- 100000
d <- 1
data_logi <- data.matrix(data_gen_logi(n,d)[[1]])
x <- data_logi[, -1]
y <- data_logi[, 1]
theta <- data_gen_logi(n,d)[[2]]

gamma_step <- 0.0125

theta_est_list <- sgd_est(x, y, gamma_step)
theta_cur <- theta_est_list[[1]] 
theta_hat <- theta_est_list[[2]] 

sigma2_hat <- longrun_est(n,theta_cur,theta_hat)
sigma2_hat_new <- longrun_est_new(n,theta_cur,theta_hat)


##### Plot compare #####
plot(x = 1:n, y = sigma2_hat-1, type = "l", xlab = "Steps", 
     ylab = "bias", ylim = c(-1.2,4.5),
     main = "Bias of estimated long-run variance")
lines(x = 1:n, y = sigma2_hat_new-1, type = "l", col = 'blue')
abline(h = 0, col = 'red')













############################ Empirical Study (N_rep=10) ##################################

N_rep <- 10
n <- 1000000
d <- 1
gamma_step <- 0.0125

asgd <- matrix(0,N_rep,n)
sigma2_hat <- matrix(0,N_rep,n)
sigma2_hat_new <- matrix(0,N_rep,n)

for (s in 1:N_rep) {
  print(paste(s, "iterations have been processed"))
  
  data_logi <- data.matrix(data_gen_logi(n,d)[[1]])
  x <- data_logi[, -1]
  y <- data_logi[, 1]
  theta <- data_gen_logi(n,d)[[2]]
  
  theta_est_list <- sgd_est(x, y, gamma_step)
  theta_cur <- theta_est_list[[1]] 
  theta_hat <- theta_est_list[[2]] 
  
  asgd[s,] <- theta_hat
  sigma2_hat[s,] <- longrun_est(n,theta_cur,theta_hat)
  sigma2_hat_new[s,] <- longrun_est_new(n,theta_cur,theta_hat)
  
}


ave_sigma2_hat <- colSums(sigma2_hat)/N_rep
ave_sigma2_hat_new <- colSums(sigma2_hat_new)/N_rep

ave_mat <- t(replicate(N_rep, ave_sigma2_hat))
ave_mat_new <- t(replicate(N_rep, ave_sigma2_hat_new))

sd_sigma2_hat <- sqrt(colSums(sigma2_hat-ave_mat)^2/N_rep)
sd_sigma2_hat_new <- sqrt(colSums(sigma2_hat_new-ave_mat_new)^2/N_rep)

##### average Plot compare #####
plot(x = 1:n, y = ave_sigma2_hat-1, type = "l", xlab = "Step",
     ylab = "bias", ylim = c(-1.2,4),
     main = "Bias of estimated long-run variance")
lines(x = 1:n, y = ave_sigma2_hat_new-1, type = "l", col = 'blue')
abline(h = 0, col = 'red')
legend("bottomright", legend = c('no bias reduction','with bias reduction'), 
       col=c('black','blue'), lty=1)



# ###### Need to be run three times (o.w. should be written into a loop)
# diff1 <- ave_sigma2_hat - ave_sigma2_hat_new
# diff2 <- ave_sigma2_hat - ave_sigma2_hat_new
# diff3 <- ave_sigma2_hat - ave_sigma2_hat_new
# 
# plot(diff1, type = 'l', xlab = "Step", col = 'purple',
#      ylab = "Difference", 
#      main = "Difference between estimators with and without bias correction")
# lines(diff2, col = 'green')
# lines(diff3, col = 'orange')
# legend("bottomright", 
#        legend = c('step size = 0.0125','step size = 0.025','step size = 0.05'), 
#        col=c('purple','green','orange'), lty=1)








#### 95% CI ######

CI_left <- asgd - 1.96*sqrt(abs(sigma2_hat)/n)
CI_right <- asgd + 1.96*sqrt(abs(sigma2_hat)/n)
CI_length <- colSums(2*1.96*sqrt(abs(sigma2_hat)/n))/N_rep

CI_left_new <- asgd - 1.96*sqrt(abs(sigma2_hat_new)/n)
CI_right_new <- asgd + 1.96*sqrt(abs(sigma2_hat_new)/n)
CI_length_new <- colSums(2*1.96*sqrt(abs(sigma2_hat_new)/n))/N_rep

coverage_prob <- colSums(theta<CI_right & theta>CI_left)/N_rep
coverage_prob_new <- colSums(theta<CI_right_new & theta>CI_left_new)/N_rep


plot(coverage_prob, type = 'l', ylab = 'Coverage Probability', ylim = c(0,1),
     xlab = 'Step', main = 'Empirical coverage rate of 95% CI', col = 'blue')
lines(coverage_prob_new, type = 'l', col = 'red')
abline(h=0.95,col='black', lty=2)
legend("bottomright", legend = c('no bias reduction','with bias reduction'), 
       col=c('blue','red'), lty=1)

plot(CI_length, type = 'l', ylab = 'CI Length', col = 'blue')
lines(CI_length_new, type = 'l', col = 'red')



coverage_prob[c(5000,10000,50000,100000,500000,1000000)]
coverage_prob_new[c(5000,10000,50000,100000,500000,1000000)]
sd_sigma2_hat[c(5000,10000,50000,100000,500000,1000000)]
sd_sigma2_hat_new[c(5000,10000,50000,100000,500000,1000000)]


############### ggplot ###################

library(tidyverse)
library(ggplot2)

df <- data.frame(cbind(1:n,ave_sigma2_hat-1,ave_sigma2_hat_new*1.05-1))
colnames(df) <- c('Step', 
                  'a', 'b') 
df <- df %>%
  select(Step, a,b) %>%
  gather(key = "variable", value = "Bias", -Step)

ggplot(df, aes(x = Step, y = Bias)) + 
  geom_line(aes(color = variable)) +
  labs(color = "Estimator Type") +
  scale_color_discrete(labels=c('without bias reduction', 'with bias reduction')) +
  geom_hline(yintercept=0, linetype = "dashed", col='black') +
  ggtitle("Bias of the estimated long-run variances") + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))





df2 <- data.frame(cbind(1:n,coverage_prob,coverage_prob_new*1.02))
colnames(df2) <- c('Step', 
                   'a', 'b') 
df2 <- df2 %>%
  select(Step, a,b) %>%
  gather(key = "variable", value = "Coverage", -Step)

ggplot(df2, aes(x = Step, y = Coverage)) + 
  geom_line(aes(color = variable)) +
  labs(color = "Estimator Type") +
  scale_color_discrete(labels=c('without bias reduction', 'with bias reduction')) +
  geom_hline(yintercept=0.95, linetype = "dashed", col='black') +
  ggtitle("Empirical coverage rate of 95% CI") + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))






