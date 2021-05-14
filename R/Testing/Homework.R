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


#Compute p(Z[1] | Z[2], Z[3]) -------------------------------------------------
w_12 <- sum(networks[, 1])
w_13 <- sum(networks[, 2])
w_23 <- sum(networks[, 3])

prior_1 <- 1/number_blocks
prior_2 <- 1/number_blocks

#First part, assigning node 1 to block 1 --------------------------------------
Z[1] <- 1
p1 <- theta[Z[1], Z[2]] ^ w_12 * (1 - theta[Z[1], Z[2]])^{N - w_12} *
  theta[Z[1], Z[3]] ^ w_13 * (1 - theta[Z[1], Z[3]])^{N - w_13} *
  theta[Z[2], Z[3]] ^ w_23 * (1 - theta[Z[2], Z[3]])^{N - w_23} *
  prior_1

#Second part, assigning node 1 to block 2 -------------------------------------
Z[1] <- 2
p2 <- theta[Z[1], Z[2]] ^ w_12 * (1 - theta[Z[1], Z[2]])^{N - w_12} *
  theta[Z[1], Z[3]] ^ w_13 * (1 - theta[Z[1], Z[3]])^{N - w_13} *
  theta[Z[2], Z[3]] ^ w_23 * (1 - theta[Z[2], Z[3]])^{N - w_23} *
  prior_2

# The probability p(Z[1] = 1 | Z[2], Z[3]) ------------------------------------
p1 / (p1  + p2)
# The part "theta[Z[2], Z[3]] ^ w_23 * (1 - theta[Z[2], Z[3]])^{N - w_23}"
# cancels here. Assignment 1: Check this by hand!

#Compute p(Z[1] | Z[2], Z[3]) now using logs ----------------------------------
#First part, assigning node 1 to block 1 --------------------------------------
Z[1] <- 1
l1 <- w_12 * log(theta[Z[1], Z[2]]) +
  {N - w_12} * log(1 - theta[Z[1], Z[2]]) +
  w_13 * log(theta[Z[1], Z[3]]) + {N - w_13} * log(1 - theta[Z[1], Z[3]]) +
  w_23 * log(theta[Z[2], Z[3]]) + {N - w_23} * log(1 - theta[Z[2], Z[3]]) +
  log(prior_1)

#Second part, assigning node 1 to block 2 -------------------------------------
Z[1] <- 2
l2 <- w_12 * log(theta[Z[1], Z[2]]) +
  {N - w_12} * log(1 - theta[Z[1], Z[2]]) +
  w_13 * log(theta[Z[1], Z[3]]) + {N - w_13} * log(1 - theta[Z[1], Z[3]]) +
  w_23 * log(theta[Z[2], Z[3]]) + {N - w_23} * log(1 - theta[Z[2], Z[3]]) +
  log(prior_2)



# The probability p(Z[1] = 1 | Z[2], Z[3]) ------------------------------------
exp(l1) / (exp(l1)  + exp(l2))

exp(l2) / (exp(l1) + exp(l2))
#Check this probability (simulating) also by setting Z[1] to 2 in simulating --
# data ------------------------------------------------------------------------

Indexing <- indexing(3)
index1 <- Indexing[Indexing[,2] == 1, -2]
index2 <- Indexing[Indexing[,3] == 1, -3]
index <- rbind(index1, index2)
w_plus <- c(w_12, w_13, w_23)

prob <- numeric()
prob_parts <- numeric()
for(blocks in 1:2){

  prob[blocks] <- log(theta[blocks, Z[index[, 2]]]) %*% w_plus[index[ ,1]]  + log(1-theta[blocks, Z[index[,2]]]) %*% (N - w_plus[index[, 1]])

}
prob <- prob - max(prob) - log(sum(exp(prob - max(prob))))
az <- numeric(length = 3)
az[1] <- 1:2 %*% rmultinom(1, 1, exp(prob))


#### it works !!!!


#Now, a Gibbs sampler for this one case is trivial ----------------------------
aZ <- Z
estimateZ1 <- NULL
for(iteration in 1:1e5) {
  aZ[1] <- 2 - rbinom(n = 1, size = 1, prob = exp(l1) / (exp(l1)  + exp(l2)))
  estimateZ1[length(estimateZ1) + 1] <- aZ[1]
}
rbind(round(table(estimateZ1) / iteration, 4),
      round(c(exp(l1) / (exp(l1)  + exp(l2)), exp(l2) / (exp(l1)  + exp(l2))), 4))

#It works! Now we need to expand this to the situation that we need to sample -
# both Z[1] and Z[2]. That is, you need to complete the following skeleton ----

Z[2] <- 1
l1_n2 <- w_12 * log(theta[Z[1], Z[2]]) +
  {N - w_12} * log(1 - theta[Z[1], Z[2]]) +
  w_13 * log(theta[Z[1], Z[3]]) + {N - w_13} * log(1 - theta[Z[1], Z[3]]) +
  w_23 * log(theta[Z[2], Z[3]]) + {N - w_23} * log(1 - theta[Z[2], Z[3]]) +
  log(prior_1)

#Second part, assigning node 2 to block 2 -------------------------------------
Z[2] <- 2
l2_n2 <- w_12 * log(theta[Z[1], Z[2]]) +
  {N - w_12} * log(1 - theta[Z[1], Z[2]]) +
  w_13 * log(theta[Z[1], Z[3]]) + {N - w_13} * log(1 - theta[Z[1], Z[3]]) +
  w_23 * log(theta[Z[2], Z[3]]) + {N - w_23} * log(1 - theta[Z[2], Z[3]]) +
  log(prior_2)

exp(l1_n2) / (exp(l2_n2) + exp(l1_n2))
exp(l2_n2) / (exp(l2_n2) + exp(l1_n2))

#Check with the algorithm
Indexing <- indexing(3)
index1 <- Indexing[Indexing[,2] == 2, -2]
index2 <- Indexing[Indexing[,3] == 2, -3]
index <- rbind(index1, index2)
w_plus <- c(w_12, w_13, w_23)

prob <- numeric()
prob_parts <- numeric()
for(blocks in 1:2){

  prob[blocks] <- log(theta[blocks, Z[index[, 2]]]) %*% w_plus[index[ ,1]]  + log(1-theta[blocks, Z[index[,2]]]) %*% (N - w_plus[index[, 1]])

}
prob <- prob - max(prob) - log(sum(exp(prob - max(prob))))
az <- numeric(length = 3)
az[1] <- 1:2 %*% rmultinom(1, 1, exp(prob))


aZ <- Z
estimateZ1 <- NULL
estimateZ2 <- NULL
for(iteration in 1:1e5) {
  aZ[1] <- 2 - rbinom(n = 1, size = 1, prob = exp(l1) / (exp(l1)  + exp(l2)))
  aZ[2] <- 2 - rbinom(n = 1, size = 1, prob = exp(l1_n2) / (exp(l1_n2)  + exp(l2_n2)))
  estimateZ1[length(estimateZ1) + 1] <- aZ[1]
  estimateZ2[length(estimateZ2) + 1] <- aZ[2]
}
#Observe that here both Z[1] and Z[2] are unknown. Their full-conditionals need
# to be derived in a manner similar to the one we followed above. -------------

# Bonus assignment (not sure if my hunch is correct) --------------------------
# A maths question: Can you derive the probability p(Z[1] | Z[3]) by hand? Note
# that it is the conditional based on the joint distribution p(Z[1], Z[3]) with
# Z[2] "integrated out". The distribution from "estimateZ1" should align with -
# p(Z[1] | Z[3]). Do you also know why? ---------------------------------------
