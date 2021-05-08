#first run the algorithm and save the mcmc output

m <- matrix(0, iter, (K * (K+1))/2)


for(i in 1:1000){
 m[i, ] <- theta_c[[i]][lower.tri(theta_c[[i]], diag = T)]
}



par(mfrow = c(2,5))
for(i in 1:ncol(m)){
  plot(m[,i], type = "l")
}




par(mfrow = c(2,5))
for(i in 1:ncol(m)){

  hist(m[,i])
}



par(mfrow = c(2,5))
for(i in 1:p){

  hist(z_c[,i])
}


colMeans(m[400:1000,])
Theta_true[lower.tri(Theta_true, diag = TRUE)]

m_hat <- matrix(0, K, K)
m_hat[lower.tri(m_hat, diag = TRUE)] <- colMeans(m[400:1000,])
m_hat
Theta_true

##########################↓ artificial identifiable constrain

#simple IC 1>2>3>4>5>6

m_reordered <- matrix(0, iter, (K * (K+1))/2)
for(reorder in 1:nrow(m)){
  m_reordered[reorder,] <- sort(m[reorder,], decreasing = TRUE)
}



par(mfrow = c(2,5))
for(i in 1:ncol(m)){
  plot(m_reordered[,i], type = "l")
}




par(mfrow = c(3, 5))
for(i in 1:ncol(m)){

  hist(m_reordered[,i])
}


colMeans(m_reordered)
sort(Theta_true[lower.tri(Theta_true, diag = TRUE)], decreasing = TRUE)



