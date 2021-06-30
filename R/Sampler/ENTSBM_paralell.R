library(igraph)
library(parallel)
library(svMisc)
#Simulation settings ----------------------------------------------------------
set.seed(2)   #for reproducibility
N <- 500     #number of persons
p <- 10       #number of nodes
m <- 1        #number of observations on node states
alpha <- beta <- 1   #hyperparameters theta
nIter <- 1  #number of iterations
K <- 3        #number of blocks

#Declare data vectors ---------------------------------------------------------
person <- rep(1:N, each = m)
X <- matrix(0, nrow = length(person), ncol = p)
w <- matrix(0, nrow = N, ncol = choose(p, 2))

#Generate encompassing network model parameters -------------------------------
Theta <- matrix(0, nrow = K, ncol = K)
for(b in 1:(K-1)) {
  for(bb in (b+1):K) {
    Theta[b, bb] <- rbeta(1, 1, 8)
    Theta[bb, b] <- Theta[b, bb]
  }
}
diag(Theta) <-  rbeta(K, 5, 2)

z <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
z <- as.vector(z)

mu <- rnorm(n = p, 0, 1)

#Generate idiographic topologies ----------------------------------------------
theta <- numeric(length = choose(p, 2))
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

#Generate node states (color the graph) ---------------------------------------
delta <- matrix(0, nrow = N, ncol = choose(p, 2))
y_plus <- matrix(0, nrow = N, ncol = p)
nr_obs_person <-matrix(m, nrow = 1, ncol = N)
for(v in 1:N) {
  #Build network from vector w[v, ] -------------------------------------------
  omega <- matrix(0, nrow = p, ncol = p)
  counter <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      counter <- counter + 1
      omega[s, t] <- omega[t, s] <- w[v, counter]
    }
  }
  
  #Use igraph to determine clusters -------------------------------------------
  G <- graph_from_adjacency_matrix(omega, 
                                   mode = "undirected")
  network <- igraph::components(G)
  
  #Color the graph ------------------------------------------------------------
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
  #Use "person" to index the observations of an individual v ------------------
  X[person == v, ] <- x
  counter <- 0
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      counter <- counter + 1
      if(sum(x[, s] * x[, t]) == m)
        delta[v, counter] <- 1
    }
  }
  y_plus[v, ] <- colSums(1/2 + x/2)
}

##########################################################################################???
# -----------------------------------------------------------------------------
# The Gibbs sampler -----------------------------------------------------------
# -----------------------------------------------------------------------------

#Generate starting values -----------------------------------------------------
aTheta <- matrix(0, nrow = K, ncol = K)
for(b in 1:(K-1)) {
  for(bb in (b+1):K) {
    aTheta[b, bb] <- rbeta(1, 1, 20)
    aTheta[bb, b] <- aTheta[b, bb]
  }
}
diag(aTheta) <-  rbeta(K, 1, 10)

az <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
az <- as.vector(az)
amu <- rnorm(p)
aw <- matrix(0, nrow = N, ncol = choose(p, 2))
atheta <- numeric(length = choose(p, 2))
counter <- 0
for(i in 1:(p-1)) {
  for(j in (i+1):p) {
    counter <- counter + 1
    aw[, counter] <- rbinom(n = N,
                           size = 1,
                           prob = Theta[az[i], az[j]])
    atheta[counter] <- aTheta[az[i], az[j]]
  }
}
#Declare vectors for the posterior means and output ---------------------------
eap_t <- rep(0, length = (K * (K+1))/2)
eap_m <- rep(0, length = p)
acceptance_rate <- rep(0, length = p)

Theta_c <- matrix(0, nrow = nIter, ncol = (K * (K+1))/2)
Mu_c <- matrix(0, nrow = nIter, ncol = p)
z_c <- matrix(0, nrow = nIter, ncol = p)
w_c <- w_c <- array(data = 0, dim = c(N, choose(p, 2), nIter))

#utility
numCores <- detectCores() - 1

sample.topology <- function(v) {
  topology <- sample_person_topology_dac_delta(delta[v, ],
                                               w = aw[v, ], 
                                               node_threshold = amu,
                                               edge_probability = atheta, 
                                               m = k) 
  return(list(topology = topology, person = v)) 
}


for(iter in 1:nIter) {
  #Step 0: progressbar and vectorform of the theta from the previous step
  progress(iter, max.value = nIter)
  #updating the vector form of Theta
  atheta <- vector(length = choose(p, 2))
  counter <- 0
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      counter <- counter + 1
      atheta[counter] <- aTheta[az[s], az[t]]
    }
  } #using the inital blockmembership simulate the Theta edge weight vector
  #Step 1: sample the topologies ----------------------------------------------
  
  list_topologies <- mclapply(X = 1:N, 
                              FUN = sample.topology, 
                              mc.cores = numCores)
  for(v in 1:N){
    aw[list_topologies[[v]]$person, ] <- list_topologies[[v]]$topology
  }
  
  #Step 2: sample the main effects  -------------------------------------------
  #First, assess topologies
  topologies <- topology_clustering(N = N, p = p, w = aw)
  #topologies <- topology_clustering_parallel(N = N, 
  #                                           p = p, 
  #                                           w = aw, 
  #                                           numCores = numCores)
  #Now, update the thresholds
  current_state <- amu
  for(i in 1:p) {
    amu[i] <- sample_node_threshold_dac(nr_obs_person = nr_obs_person, 
                                        y_plus = y_plus,
                                        p = p,
                                        clustering = topologies$clustering, 
                                        clustersize = topologies$clustersize,
                                        node_threshold = amu,
                                        node = i)
    if(current_state[i] != amu[i]) {
      acceptance_rate[i] <- (iter - 1) * acceptance_rate[i] /iter + 1 / iter
    } else {
      acceptance_rate[i] <- (iter - 1) * acceptance_rate[i] /iter 
    }
  }
  
  #Step 3: sample the wiring probabilities ------------------------------------
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
  
  #Step 4: sample z------------------------------------------------------------
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
  
  if(iter < burn_in){
    z_c[iter, ] <- az
  } else {
    new_az <- online_relabelling(z_c[1:burn_in, ], az, K)
    z_c[iter, ] <- new_az
  }
  
  #Compute posterior mean wiring probabilities and main effects --------------- 
  eap_t <- eap_t * (iter - 1) / iter + atheta / iter
  eap_m <- eap_m * (iter - 1) / iter + amu / iter
  
  Theta[iter, ] <- atheta
  Mu[iter, ] <- amu
  
  #Screen output --------------------------------------------------------------
  if((iter - 1)%%5 == 0) {
    par(mfrow = c(2, 1))
    
    par(cex.main = 1.5, mar = c(5, 6, 1, 1) + 0.1, mgp = c(3.5, 1, 0), 
        cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
    plot(eap_t, theta, col = "black", pch = 21, bg = "grey", cex = 2,
         xlim = c(0, 1), ylim = c(0, 1), ylab = "", xlab = "", axes = FALSE)
    abline(0, 1, lwd = 2, lty =2 , col = "gray")
    par(las = 0)
    mtext(expression(paste("Posterior Mean ", theta)), side = 1, line = 2.5, cex = 1.5)
    mtext(expression(paste("True Value ", theta)), side = 2, line = 3.7, cex = 1.5)
    
    
    lim <- c(min (c(eap_m, mu))-.01, max (c(eap_m, mu))+.01)
    par(cex.main = 1.5, mar = c(5, 6, 3, 1) + 0.1, mgp = c(3.5, 1, 0), 
        cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
    plot(eap_m, mu, col = "black", pch = 21, bg = "grey", cex = 2,
         xlim = lim, ylim = lim, ylab = "", xlab = "", axes = FALSE)
    abline(0, 1, lwd = 2, lty =2 , col = "gray")
    par(las = 0)
    mtext(expression(paste("Posterior Mean ", mu)), side = 1, line = 2.5, cex = 1.5)
    mtext(expression(paste("True Value ", mu)), side = 2, line = 3.7, cex = 1.5)
    
  }
  print(paste(iter, "  ", round(100 * mean(acceptance_rate))))  
}
