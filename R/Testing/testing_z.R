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

az <- c(2, 1, 2)


#################TESTING MY CODE
p_c <- array(0, dim = c(100000, 3, 2))
z_c <- matrix(0, 100000, 3)
for(iter in 100000){
  w_plus <- colSums(networks)
  p_c_m <- matrix(0, 3, 2)

  for(nodes in 1:3){
    index1 <- Indexing[Indexing[,2] == nodes, -2]
    index2 <- Indexing[Indexing[,3] == nodes, -3]
    index <- rbind(index1, index2)

    prob <- numeric()
    prob_parts <- numeric()
    for(blocks in 1:2){

      prob[blocks] <- log(theta[blocks, az[index[, 2]]]) %*% w_plus[index[ ,1]]  + log(1-theta[blocks, az[index[,2]]]) %*% (N - w_plus[index[, 1]])

    }
    prob <- prob - max(prob) - log(sum(exp(prob - max(prob))))
    az[nodes] <- 1:2 %*% rmultinom(1, 1, exp(prob))

    #save clustering probabilites
    p_c_m[nodes, ] <- exp(prob) / sum(exp(prob))
  }

  p_c[iter,,] <- p_c_m
  z_c[iter, ] <- az

}

par(mfrow = c(1,3))
for(i in 1:3){
  hist(z_c[,i])
}

