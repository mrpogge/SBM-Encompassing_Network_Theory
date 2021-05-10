#first run the algorithm and save the mcmc output

theta_c <- res$theta
m <- matrix(0, 1000, (K * (K+1))/2)


for(i in 1:1000){
 m[i, ] <- theta_c[[i]][lower.tri(theta_c[[i]], diag = T)]
}



par(mfrow = c(3,5))
for(i in 1:ncol(m)){
  plot(m[,i], type = "l")
}




par(mfrow = c(3,5))
for(i in 1:ncol(m)){

  hist(m[,i])
}



par(mfrow = c(4,5))
for(i in 1:p){

  hist(z_c[,i])
}


colMeans(m[400:1000,])
Theta_true[lower.tri(Theta_true, diag = TRUE)]

m_hat <- matrix(0, K, K)
m_hat[lower.tri(m_hat, diag = TRUE)] <- colMeans(m[400:1000,])
m_hat
Theta_true

colMeans(m_atheta[200:1000,])
Theta_true[lower.tri(Theta_true, diag = TRUE)]

m_atheta_hat <- matrix(0, K, K)
m_hat[lower.tri(m_atheta_hat, diag = TRUE)] <- colMeans(m_atheta[400:1000,])
m_hat
Theta_true
