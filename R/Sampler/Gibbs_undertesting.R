SBM_Gibbs <- function(aw, p, K, N, m = 1, iter){
  
  ####indexing#################################################
  indexing <- function(p) {
    index <- matrix(0, nrow = p * (p - 1) / 2, ncol = 3)
    counter <- 0
    for(i in 1:(p - 1)) {
      for(j in (i + 1):p) {
        counter <- counter + 1
        index[counter, 1] <- counter
        index[counter, 2] <- i
        index[counter, 3] <- j
      }
    }
    return(index)
  }
  
  Indexing <- indexing(p)
  
  ####### setting beta hyperparameters ######################
  alpha <- beta <- 1
  
  ####### creating a random allocation of the nodes to blocks ##############
  az <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
  az <- as.vector(az)
  
  ######## creating a matrix for the edge probabilites (0-s) ###############
  aTheta <- matrix(0, nrow = K, ncol = K)
  for(b in 1:(K-1)) {
    for(bb in (b+1):K) {
      aTheta[b, bb] <- runif(n = 1, min = 0, max = 0.5)
      aTheta[bb, b] <- aTheta[b, bb]
    }
  }
  diag(aTheta) <-  runif(n = K, min = 0, max = 0.5)
  
  ##############Gibbs-sampler ################################################################
  
  #################containers############################
  theta_c <- list()
  z_c <- matrix(0, nrow = iter, ncol = p)
  p_c <- array(0, dim = c(iter, p, K))
  
  for(iter in 1:iter){
    
    counter <- 0
    for(b in 1:K) {
      for(bb in b:K) {
        w_plus <- 0
        w_min <- 0
        if(b == bb) { #within block 
          for(s in 1:(p-1)) {
            for(t in (s+1):p) {
              if(az[s] == b & az[t] == b) {
                counter <- Indexing[Indexing[, 2] == s & Indexing[, 3] == t, 1]
                w_plus <- w_plus + sum(aw[,counter]) #number of connections across all individual in the given i,j node pair for every node pair witihn cluster
                w_min <- w_min + N - sum(aw[,counter]) 
              }
            }
          }
        } else { #between block (upper triangular)
          for(s in 1:(p-1)) {
            for(t in (s+1):p) {
              if(az[s] == b & az[t] == bb) {
                counter <- Indexing[Indexing[, 2] == s & Indexing[, 3] == t, 1]
                w_plus <- w_plus + sum(aw[,counter])
                w_min <- w_min + N - sum(aw[,counter])
              }
              if(az[s] == bb & az[t] == b) {
                counter <- Indexing[Indexing[, 2] == s & Indexing[, 3] == t, 1]
                w_plus <- w_plus + sum(aw[,counter])
                w_min <- w_min + N - sum(aw[,counter])
              }
            }
          }        
        }
        aTheta[b, bb] <- rbeta(n = 1, 
                               shape1 = alpha + w_plus,
                               shape2 = beta + w_min) # is it right? we need to rename that to alpha and beta
        aTheta[bb, b] <- aTheta[b, bb] #simulate theta matrix from the conditional, since we used beta-binomial conjugacy the conditional is Beta(a + w_plus, b + w_min)
      }
    }
    
    #save the current theta matrix to the container ######
    theta_c[[iter]] <- aTheta
    
    ##################Z conditional for details see. SBM-Encompassing_Network_Theory/R/Sampler/README############################
    w_plus <- colSums(aw)
    p_c_m <- matrix(0, p, K)
    
    for(nodes in 1:p){
      index1 <- Indexing[Indexing[,2] == nodes, -2]
      index2 <- Indexing[Indexing[,3] == nodes, -3]
      index <- rbind(index1, index2)
      
      prob <- numeric()
      prob_parts <- numeric()
      for(blocks in 1:K){
        
        prob[blocks] <- log(aTheta[blocks, az[index[, 2]]]) %*% w_plus[index[ ,1]]  + log(1-aTheta[blocks, az[index[,2]]]) %*% (N - w_plus[index[, 1]])
        
      }
      prob <- prob - max(prob) - log(sum(exp(prob - max(prob))))
      az[nodes] <- 1:K %*% rmultinom(1, 1, exp(prob))
      
      #save clustering probabilites
      p_c_m[nodes, ] <- exp(prob) / sum(exp(prob))
    }
    
    p_c[iter,,] <- p_c_m
    z_c[iter, ] <- az
    
  }
  #creating MCMC_SBM S3 class
  l <- list(theta = theta_c,
            z = z_c,
            class_prob = p_c,
            iter = iter,
            N_Blocks = K,
            N_nodes = p,
            Sample_size = N,
            N_observation = m,
            Obs_Topology = aw)
  
  class(l) <- "MCMC_SBM"
  return(l)
}
