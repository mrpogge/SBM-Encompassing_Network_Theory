library(label.switching)

relabelR <- function(obj){
  require(label.switching)

  #detect relabelling
  run <- stephens(obj$class_prob)

  #save the permutations
  perm <- run$permutations

  iter <- obj$iter
  K <- obj$N_Blocks
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

  #permuting z vector
  z_relab <- matrix(0, iter, obj$N_nodes)
  for(perm_rows in 1:nrow(perm)){
    for(zs in 1:obj$N_nodes){
      z_relab[perm_rows, zs] <- which(perm[perm_rows, ] == obj$z[perm_rows, zs])
    }
  }

  l <- list(theta = m_atheta,
            z = z_relab,
            class_prob = obj$class_prob,
            iter = obj$iter,
            N_Blocks = obj$N_Blocks,
            N_nodes = obj$N_nodes,
            Sample_size = obj$Sample_size,
            N_observation = obj$N_observation,
            Obs_Topology = obj$Obs_Topology,
            run = run)
  class(l) <- "MCMC_SBM"
  return(l)

}







