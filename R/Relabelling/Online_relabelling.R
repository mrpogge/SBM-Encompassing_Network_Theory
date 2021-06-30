# #saving trial mcmc for z
# z_MCMC <- sim_results[[10]]$estimate$z
# plot(z_MCMC[,9])
#
# #easy example
# z_proba <- c(1,1,3,2,2,1,3, 2,2,1,3,3,2,1, 2,2,1,3,3,2,1, 3,3,2,1,1,3,2, 2,2,1,3,3,2,1)
# z_proba <- matrix(data = z_proba, nrow = 5, ncol = 7, byrow = TRUE)
# vec <- c(1,1,3,2,2,1,3)
# perm_0 <- c(1,2)
# labs <- c(2,1)
# new_vec <- matrix(0, nrow = 1, ncol = length(vec))
# new_vec[which(vec == perm_0[1])] <- labs[1]
#
#
# z_proba <- z_MCMC[2000:4000, ]
#
# z_new <- z_MCMC[7000, ]
#
# ######################################################################################################################
#
# K <- 2
# all_perm <- permut(1:K)
# relabelR_online <- function(z_prev, z_new, all_perm, K){
#
#   #calculating the distance function of two steps of the sampler
#   class_dist <- function(z_u, z_t){
#
#     delta <- 0
#     delta <- sum(z_u - z_t != 0)
#
#     return(delta)
#   }
#
#   #saving all of the permutations
#   class_perm_prop <- function(z_u, all_perm, K){
#
#     require(randtests)
#     perm_0 <- 1:K
#     all_perm <- all_perm
#
#     new_z_u <- matrix(0, nrow = nrow(all_perm), ncol = length(z_u))
#     for(prm in 1:nrow(all_perm)){
#
#       for(lbl in 1:length(perm_0)){
#
#         new_z_u[prm, which(z_u == perm_0[lbl])] <- all_perm[prm, lbl]
#
#       }
#     }
#     return(new_z_u)
#   }
#
#   #calculating the cost matrix per permutation
#   all_alloc <- class_perm_prop(z_new, all_perm, K)
#   cost <- matrix(0, nrow = nrow(z_prev), ncol = nrow(all_alloc))
#   for(steps in 1:nrow(z_prev)){
#     for(allocs in 1:nrow(all_alloc)){
#
online_relabelling <- function(z_prev, z_new, K){

  #requirements
  require(RcppHungarian)

  #creating full contignency tables
  full_table <- function(z_prev_i, z_new, K) {
    m <- matrix(0, K, K)
    for(i in 1:K){
      for(j in 1:K){
        m[i,j] <- sum(z_new != i & z_prev_i == j)
      }
    }
    return(m)
  }

  #creating the cost function
  cost_m <- matrix(0, nrow = K, ncol = K)

  for(rw in 1:nrow(z_prev)){
    cost_m <- cost_m + full_table(z_prev[rw, ], z_new, K)
  }

  #searching for the best match with the hungarian algorithm
  solved <- HungarianSolver(cost_m)

  new_lab <- solved$pairs[,2]

  #relabelling the current z proposal
  perm0 <- solved$pairs[,1]

  new_z_u <- numeric(length = length(z_new))
  for(lbl in 1:length(perm0)){

    new_z_u[which(z_new == perm0[lbl])] <- new_lab[lbl]

  }

  return(new_z_u)

}

