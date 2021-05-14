########################################################################################
#Gibbs sampler for the Stochastic Block Model
#This R script contains an example of the sampler (currently not in function form)
#First we create a random SBM network with SBM_sim for further info see SBM-Encompassing_Network_Theory/R/Network_simulation/SBM_sim
#Then we set the initial values for the sampler, see. 'initialisation'

#################SBM simulation######################################################
SBM_sim <- function(p, K, N = 1, alpha, beta){

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
  ########################variable assignment##########################
  Index <- indexing(p)

  theta <- vector(length = choose(p, 2)) #edge prob vector form
  Theta <- matrix(0, K, K) #edge probability matrix
  w <- matrix(0, N, choose(p, 2)) #individual topology
  W <- matrix(0, p, p) #adj matrix


  ###### sampling cluster membership form a multinomial distribution ###############
  z <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
  z <- as.vector(z)

  ############## sampling for edge probabilities from a beta distribution with hyperparameters alpha beta ######################
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

  res <- list(nodes = p, K = K, N = N, alpha = alpha, beta = beta,  Theta = Theta, theta = Theta, W = W, w = w, Index = Index, z = z)
  return(res)
}
#### initialisation################################################################

t_network <- SBM_sim(10, 3, 10, 2, 1)
K <- t_network$K
p <- t_network$nodes
aw <- t_network$w
N <- t_network$N

Theta_true <- t_network$Theta
aTheta <- Theta_true
z_true <- t_network$z
#az <- t_network$z
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

###### creating a random allocation of the nodes to blocks ##############
az <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
az <- as.vector(az)

######## creating a matrix for the edge probabilites (0-s) ###############
# aTheta <- matrix(0, nrow = K, ncol = K)
# for(b in 1:(K-1)) {
#   for(bb in (b+1):K) {
#     aTheta[b, bb] <- runif(n = 1, min = 0, max = 0.5)
#     aTheta[bb, b] <- aTheta[b, bb]
#   }
# }
# diag(aTheta) <-  runif(n = K, min = 0, max = 0.5)

##############Gibbs-sampler ################################################################

#################containers############################
theta_c <- matrix(0, 100000, K*(K+1)/2)
z_c <- matrix(0, nrow = 100000, ncol = p)
p_c <- array(0, dim = c(100000, p, K))

for(iter in 1:100000){

  # ########theta conditional see. SBM-Encompassing_Network_Theory/R/Sampler/README
  # for(b in 1:K){
  #   for(bb in b:K){
  #     w_plus <- 0
  #     w_min <- 0  #setting the power to zero for each blockpair
  #     if(b == bb){
  #       nodes <- which(az == b) #select the nodes from the block
  #       within_nodes <- length(nodes)
  #       for(i in 1:within_nodes){ #loop over every unique pair of nodes in the block
  #         for(j in i:within_nodes){
  #           if(i != j){
  #             place <- Indexing[Indexing[,2] == nodes[i] & Indexing[,3] == nodes[j], 1] #find their place in aw matrix
  #             w_plus <- w_plus + sum(aw[, place]) #calculate the number of edges across observation
  #             w_min <- w_min + N - sum(aw[,place]) #calculate the number of non-edges across observations
  #           } else{
  #             next
  #           }
  #         }
  #       }
  #     } else { #between blocks
  #       nodes_b <- which(az == b) #nodes from block b
  #       nodes_bb <- which(az == bb)#nodes from block bb
  #       for(i in nodes_b){
  #         for(j in nodes_bb){ #ADD AN IF STATEMENT / I've added an if statement because you can't subset Indexing if the first value is higher
  #           if(i > j){
  #             place <- Indexing[Indexing[,2] == i & Indexing[,3] == j, 1] #find their place in aw matrix
  #             w_plus <- w_plus + sum(aw[, place]) #calculate the number of edges across observation
  #             w_min <- w_min + N - sum(aw[,place]) #calculate the number of non-edges across observation
  #           } else {
  #             place <- Indexing[Indexing[,2] == j & Indexing[,3] == i, 1] #find their place in aw matrix
  #             w_plus <- w_plus + sum(aw[, place]) #calculate the number of edges across observation
  #             w_min <- w_min + N - sum(aw[,place]) #calculate the number of non-edges across observation
  #           }
  #         }
  #       }
  #     }
  #     aTheta[b, bb] <- rbeta(n = 1,
  #                            shape1 = alpha + w_plus,
  #                            shape2 = beta + w_min)
  #     aTheta[bb, b] <- aTheta[b, bb]
  #   }
  # }
  #
  # #save the current theta matrix to the container ######
  # theta_c[iter, ] <- aTheta[lower.tri(aTheta, diag = TRUE)]

  #################Z conditional for details see. SBM-Encompassing_Network_Theory/R/Sampler/README############################
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

theta_hat <- colMeans(theta_c)
theta_hat
Theta_true

z_true
az

par(mfrow = c(2,5))
for(i in 1:10){
  hist(z_c[,i])
}

z_hat <- numeric()
for(col in 1:10){
  z_hat[col] <- median(z_c[,col])
}

