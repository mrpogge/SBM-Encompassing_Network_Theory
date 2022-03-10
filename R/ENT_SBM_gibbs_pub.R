library(igraph)
N <- 500 #number of persons
p <- 15 #number of nodes
m <- 10 #number of observations on node states
B <- K <- 5 #number of stochastic blocks
iteration <- 10

alpha <- beta <- 1 #hyperparameters

#indexing of the network
indexing <- function(p) {
  index <- matrix(0, nrow = p * (p - 1) / 2, ncol = 3)
  counter <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      counter <- counter + 1
      index[counter, 1] <- counter
      index[counter, 2] <- s
      index[counter, 3] <- t
    }
  }
  return(index)
}
Indexing <- indexing(p)

#containers for the true network
person <- rep(1:N, each = m)
X <- matrix(0, nrow = length(person), ncol = p)
w <- matrix(0, nrow = N, ncol = choose(p, 2))

#model
z <- 1:B %*% rmultinom(n = p, size = 1, prob = rep(1, B))
z <- as.vector(z)

theta <- matrix(0, nrow = B, ncol = B)
for(b in 1:(B-1)) {
  for(bb in (b+1):B) {
    theta[b, bb] <- rbeta(n = 1, shape1 = 1, shape2 =  2)
    theta[bb, b] <- theta[b, bb]
  }
}

diag(theta) <- rbeta(n = B, shape1 = 3, shape2 =  1)

mu <- rnorm(n = p, 0, 1) / sqrt(p)

counter <- 0
Theta <- vector(length = choose(p, 2))
for(s in 1:(p-1)) {
  for(t in (s+1):p) {
    counter <- counter + 1
    w[, counter] <- rbinom(n = N,
                           size = 1, 
                           prob = theta[z[s], z[t]])    
    Theta[counter] <- theta[z[s], z[t]]
  }
}


for(v in 1:N) {
  omega <- matrix(0, nrow = p, ncol = p)
  counter <- 0
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      counter <- counter + 1
      omega[s, t] <- omega[t, s] <- w[v, counter]
    }
  }
  
  G <- graph_from_adjacency_matrix(omega, 
                                   mode = "undirected")
  network <- igraph::components(G)
  
  x <- matrix(0, nrow = m, ncol = p)
  for(rep in 1:m) {
    for (c in 1:network$no) {
      tmp <- sum(mu[network$membership == c])
      prob <- exp(2* tmp) / (1 + exp(2 * tmp))        
      x[rep, network$membership == c] <- - 1 + 2 *  rbinom(n = 1,
                                                           size = 1, 
                                                           prob = prob)
    }
  }
  X[person == v, ] <- x
}
plot(mu, colMeans(X));abline(lm(colMeans(X)~mu))
#Score distribution
plot.ecdf(rowSums(X))

source("gibbs_functions(1).R")

####################################################################################################

#################INITIAL VALUES####################################################################
az <- 1:B %*% rmultinom(n = p, size = 1, prob = rep(1, B))
az <- as.vector(az)

aTheta <- matrix(0, nrow = B, ncol = B)
for(b in 1:(B-1)) {
  for(bb in (b+1):B) {
    aTheta[b, bb] <- runif(n = 1, min = 0, max = 0.5)
    aTheta[bb, b] <- aTheta[b, bb]
  }
}
diag(aTheta) <-  runif(n = B, min = 0, max = 0.5)

amu <- rnorm(p)
aw <- matrix(0, nrow = N, ncol = choose(p, 2))

counter <- 0
atheta <- vector(length = choose(p, 2))
for(s in 1:(p-1)) {
  for(t in (s+1):p) {
    counter <- counter + 1
    aw[, counter] <- rbinom(n = N,
                            size = 1, 
                            prob = aTheta[az[s], az[t]])    
    atheta[counter] <- aTheta[az[s], az[t]]
  }
}

at <- 0 * atheta
am <- 0 * amu

#####creating containers####################################################################
theta_c <- matrix(0, nrow = iteration, choose(p, 2))
z_c <- matrix(0, nrow = iteration, ncol = p)
mu_c <- matrix(0, nrow = iteration, ncol = p)

##########################################################################################♣###

##############################SAMPLER#########################################################
for(iter in 1:iteration){
  

  #updating the vector form of Theta
  atheta <- vector(length = choose(p, 2))
  
  
  counter <- 0
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      counter <- counter + 1
      atheta[counter] <- aTheta[az[s], az[t]]
    }
  } #using the inital blockmembership simulate the Theta edge weight vector
  
  theta_c[iter, ] <- atheta
  
  #sampling for the indiviudal topology
  for(v in 1:N) {
    aw[v, ] <- sample_person_topology_dac(x = X[person == v, ], #observed
                                          w = aw[v, ], #empty
                                          node_threshold = amu, #simulated from normal
                                          edge_probability = atheta)  #for every person we simulate an individual topology based on the previously simulated aTheta and the initual values
  }
  
  
  #sampling for theta
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
  
  
  ##################Z conditional for details see. SBM-Encompassing_Network_Theory/R/Sampler/README############################
  w_plus <- colSums(aw)
  
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
    
  }
  
  z_c[iter, ] <- az
  # if(iter < burn_in){
  #   z_c[iter, ] <- az
  # } else {
  #   new_az <- online_relabelling(z_c[1:burn_in, ], az, K)
  #   z_c[iter, ] <- new_az

  old <- amu
  for(i in 1:p) {
    amu[i] <- sample_node_threshold_dac(X = X,
                                        person = person,
                                        w = aw, 
                                        node_threshold = amu,
                                        i = i)
  }
  
  mu_c[iter, ] <- amu
  
}

