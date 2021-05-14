#########################BASIC 3 node 2 cluster network with 10 observations
set.seed(1)
#The settings
number_nodes <- 3
number_blocks <- 2
N <- 10

#The model.
theta <- matrix(0, nrow = 2, ncol = 2)
theta[1, 1] <- 0.5
theta[2, 2] <- 0.4
theta[1, 2] <- theta[2, 1] <- 0.1

#Node assignment.
Z <- matrix(0, nrow = 3, ncol = 1)
Z[1] <- 1
Z[2] <- 1
Z[3] <- 2

#Networks
networks <- matrix (0, nrow = N, ncol = choose(number_nodes, 2))
for(person in 1:N) {
  counter <- 0
  for(node_1 in 1:(number_nodes - 1)) {
    for(node_2 in (node_1 + 1):number_nodes) {
      counter <- counter + 1
      networks[person, counter] <- rbinom(n = 1,
                                          size = 1,
                                          prob = theta[Z[node_1], Z[node_2]])
    }
  }
}
colnames(networks) <- c("1-2", "1-3", "2-3")
head(networks)

Indexing <- indexing(3)
alpha <- beta <- 1

aTheta <- matrix(0, 2, 2)
aTheta[1, 1] <- 0.5
aTheta[2, 2] <- 0.5
aTheta[1, 2] <- aTheta[2, 1] <- 0.5

###############test for the Theta part of the conditional

##############MY CODE

theta_c <- matrix(0, 100000, 3)
for(iter in 1:100000){

  ########theta conditional see. SBM-Encompassing_Network_Theory/R/Sampler/README
  for(b in 1:2){
    for(bb in b:2){
      w_plus <- 0
      w_min <- 0  #setting the power to zero for each blockpair
      if(b == bb){
        nodes <- which(Z == b) #select the nodes from the block
        within_nodes <- length(nodes)
        for(i in 1:within_nodes){ #loop over every unique pair of nodes in the block
          for(j in i:within_nodes){
            if(i != j){
              place <- Indexing[Indexing[,2] == nodes[i] & Indexing[,3] == nodes[j], 1] #find their place in aw matrix
              w_plus <- w_plus + sum(networks[, place]) #calculate the number of edges across observation
              w_min <- w_min + N - sum(networks[,place]) #calculate the number of non-edges across observations
            } else{
              next
            }
          }
        }
      } else { #between blocks
        nodes_b <- which(Z == b) #nodes from block b
        nodes_bb <- which(Z == bb)#nodes from block bb
        for(i in nodes_b){
          for(j in nodes_bb){ #ADD AN IF STATEMENT / I've added an if statement because you can't subset Indexing if the first value is higher
            if(i > j){
              place <- Indexing[Indexing[,2] == i & Indexing[,3] == j, 1] #find their place in aw matrix
              w_plus <- w_plus + sum(networks[, place]) #calculate the number of edges across observation
              w_min <- w_min + N - sum(networks[,place]) #calculate the number of non-edges across observation
            } else {
              place <- Indexing[Indexing[,2] == j & Indexing[,3] == i, 1] #find their place in aw matrix
              w_plus <- w_plus + sum(networks[, place]) #calculate the number of edges across observation
              w_min <- w_min + N - sum(networks[,place]) #calculate the number of non-edges across observation
            }
          }
        }
      }
      aTheta[b, bb] <- rbeta(n = 1,
                             shape1 = alpha + w_plus,
                             shape2 = beta + w_min)
      aTheta[bb, b] <- aTheta[b, bb]
    }
  }

  theta_c[iter, ] <- aTheta[lower.tri(aTheta, diag = TRUE)]

}
################ MY TEST RESULTS
aTheta_hat <- colMeans(theta_c)
aTheta_hat
theta

################MAARTEN'S CODE
for(iter in 1:100000){
  counter <- 0
  for(b in 1:2) {
    for(bb in b:2) {
      w_plus <- 0
      w_min <- 0
      if(b == bb) { #within block
        for(s in 1:(3-1)) {
          for(t in (s+1):3) {
            if(Z[s] == b & Z[t] == b) {
              counter <- Indexing[Indexing[, 2] == s & Indexing[, 3] == t, 1]
              w_plus <- w_plus + sum(networks[,counter]) #number of connections across all individual in the given i,j node pair for every node pair witihn cluster
              w_min <- w_min + N - sum(networks[,counter])
            }
          }
        }
      } else { #between block (upper triangular)
        for(s in 1:(3-1)) {
          for(t in (s+1):3) {
            if(Z[s] == b & Z[t] == bb) {
              counter <- Indexing[Indexing[, 2] == s & Indexing[, 3] == t, 1]
              w_plus <- w_plus + sum(networks[,counter])
              w_min <- w_min + N - sum(networks[,counter])
            }
            if(Z[s] == bb & Z[t] == b) {
              counter <- Indexing[Indexing[, 2] == s & Indexing[, 3] == t, 1]
              w_plus <- w_plus + sum(networks[,counter])
              w_min <- w_min + N - sum(networks[,counter])
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

  theta_c[iter, ] <- aTheta[lower.tri(aTheta, diag = TRUE)]

}
############ MAARTEN'S TEST RESULTS
aTheta_hat <- colMeans(theta_c[10000:100000, ])
aTheta_hat
theta







