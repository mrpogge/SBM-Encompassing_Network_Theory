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

