library(label.switching)

relabelR <- function(obj){
  require(label.switching)

  #detect relabelling
  run <- stephens(obj$class_prob)

  #save the permutations
  perm <- run$permutations

  iter <- obj$iter
  #select the iterations with relabelling
  switch_logic <- logical(length = iter)
  for(switch in 1:iter){
    switch_logic[switch] <- all(perm[switch, -(K+1)] == 1:K)
  }

  #permuting theta matrix
  m_atheta <- matrix(0, iter, (K * (K+1))/2)
  for(i in 1:iter){

    if(switch_logic[i] == TRUE){
      m_atheta[i, ] <- obj$theta[[i]][lower.tri(obj$theta[[i]], diag = T)]
    } else {
      counter <- 0
      for(s in 1:K) {
        for(t in s:K) {
          counter <- counter + 1
          m_atheta[i, ][counter] <- obj$theta[[i]][perm[i,][s], perm[i,][t]]
        }
      }
    }
  }
  l <- list(theta = m_atheta,
            z = obj$z,
            class_prob = obj$class_prob,
            iter = obj$iter,
            N_Blocks = obj$N_Blocks,
            N_nodes = obj$N_nodes,
            Sample_size = obj$Sample_size,
            N_observation = obj$N_observation,
            Obs_Topology = obj$Obs_Topology,
            permutations = perm)
  class(l) <- "MCMC_SBM"
  return(l)

}

m_atheta <- relabelR(res)

par(mfrow = c(3,5))
for(i in 1:ncol(m)){
  plot(m_atheta[[1]][,i], type = "l")
}


par(mfrow = c(3,5))
for(i in 1:ncol(m)){

  hist(m_atheta[[1]][,i])
}


