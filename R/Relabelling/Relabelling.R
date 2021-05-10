library(label.switching)

relabelR <- function(p_c, theta_c, K, iter){
  require(label.switching)

  #detect relabelling
  run <- stephens(p_c)

  #save the permutations
  perm <- run$permutations

  #select the iterations with relabelling
  switch_logic <- logical(length = iter)
  for(switch in 1:iter){
    switch_logic[switch] <- all(perm[switch, -(K+1)] == 1:K)
  }

  #permuting theta matrix
  m_atheta <- matrix(0, iter, (K * (K+1))/2)
  for(i in 1:iter){

    if(switch_logic[i] == TRUE){
      m_atheta[i, ] <- theta_c[[i]][lower.tri(theta_c[[i]], diag = T)]
    } else {
      counter <- 0
      for(s in 1:K) {
        for(t in s:K) {
          counter <- counter + 1
          m_atheta[i, ][counter] <- theta_c[[i]][perm[i,][s], perm[i,][t]]
        }
      }
    }
  }
  l <- list(m_atheta, perm)
  return(l)

}

m_atheta <- relabelR(res$class_prob, theta_c = theta_c, K = K, iter = 1000)

par(mfrow = c(3,5))
for(i in 1:ncol(m)){
  plot(m_atheta[[1]][,i], type = "l")
}


par(mfrow = c(3,5))
for(i in 1:ncol(m)){

  hist(m_atheta[[1]][,i])
}


