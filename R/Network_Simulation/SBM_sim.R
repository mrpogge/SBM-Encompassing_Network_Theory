######################################################################
#This is an R script to simulate a network under the Stochastic Block Model with Ludkin (2020) parameterisation. 
#The function requires the following parameters:
# p = number of nodes
# K = number of blocks
# N = number of observations
# alpha, beta = hyperparameters of a Beta distribution
#The function gives a list as an output with the following structure:
# nodes = number of nodes
# K = number of blocks
# N = number of observation
# alpha, beta = hyperparameters of a Beta distribution
#Theta = matrix of edge probabilities in every block-pair
#theta = vector form of the edge probabilites in every block-pair
#W = adjacency matrix of the resulting network
#w = a dataframe with binary entries; the rows of the dataframe is the observations; the columns of the dataframe are the observed edges (order of edges given by Index)
#Index = indexing of the edges, first column the unique number of the given edge, second and third column which node is connected/disconnected from which node
#z = allocation vector; which node belongs to which block
########################################################################



######SBM simulation##############################################################
SBM_sim <- function(p, K, N = 1, alpha, beta){
 
######creating indexing###########################################################
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
  ########################variable assignment######################################
  Index <- indexing(p)
  
  theta <- vector(length = choose(p, 2)) #edge prob vector form
  Theta <- matrix(0, K, K) #edge probability matrix
  w <- matrix(0, N, choose(p, 2)) #individual topology
  W <- matrix(0, p, p) #adj matrix
  
  
  ###### sampling cluster membership form a multinomial distribution ###############
  z <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
  z <- as.vector(z)
  
  ############## sampling edge probabilities from a beta distribution with hyperparameters alpha beta ######################
  diag(Theta) <- rbeta(length(diag(Theta)), alpha, beta)
  Theta[lower.tri(Theta)] <- Theta[upper.tri(Theta)] <- rep(rbeta(1, beta, alpha), length(Theta[lower.tri(Theta)])) 
  
  counter <- 0 
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      counter <- counter + 1 
      w[, counter] <- rbinom(n = N,
                             size = 1, 
                             prob = Theta[z[i], z[j]]) 
      theta[counter] <- Theta[z[i], z[j]] 
    }
  }
  
  for(v in 1:N) {
    W <- matrix(0, nrow = p, ncol = p)
    counter <- 0
    for(i in 1:(p-1)) {
      for(j in (i+1):p) {
        counter <- counter + 1
        W[i, j] <- W[j, i] <- w[v, counter]
      }
    }
  }
	
  ###### saving results in a list ########################################################
  res <- list(nodes = p, K = K, N = N, alpha = alpha, beta = beta,  Theta = Theta, theta = Theta, W = W, w = w, Index = Index, z = z)
  return(res)
}
