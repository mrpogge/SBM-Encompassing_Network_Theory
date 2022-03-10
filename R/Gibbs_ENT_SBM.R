ENTSBM_Gibbs <- function(X, p, K, N, m = 1, iteration, burn_in, mu){

  require(svMisc)
  ########Maarten Gibbs functions####################################################################################
  sample_person_topology_dac <- function(x,
                                         w,
                                         node_threshold,
                                         edge_probability) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package `igraph' is needed for handling the topological structures.
         Please install it.", call. = FALSE)
    }
    if(is.matrix(x)) {
      p <- ncol(x)
      m <- nrow(x)
    } else {
      p <- length(x)
      m <- 1
    }

    #if the response matrix is a matrix than p will be the number of columns (nodes), m will be the number of rows (observations).

    # build the adjacency matrix omega ------------------------------------------
    omega <- matrix(0, nrow = p, ncol = p)
    counter <- 0
    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        counter <- counter + 1
        omega[s, t] <- w[counter]
        omega[t, s] <- omega[s, t]
      }
    }
    #first it creates an empty matrix and a counter variable. Than it fills the upper triangle with the values of w

    # update the edge states ----------------------------------------------------
    counter <- 0                    #set the counter to 0
    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        counter <- counter + 1

        if(is.matrix(x)) {
          if(sum(x[, s] * x[, t]) < m) {
            delta = FALSE
          } else {
            delta = TRUE
          }
        } else {
          if(x[s] * x[t] < 1) {
            delta = FALSE
          } else {
            delta = TRUE
          }
        }
        #create a matrix of logicals coming from the indicator function delta
        # node states do not match (delta[s,t] = 0) -----------------------------
        if(delta == FALSE) {
          sampled_value <- 0

          # node states do match (delta[s,t] = 1) ---------------------------------
        } else {
          # check if absence of omega[s, t] cuts the cluster in two -------------
          # determine network for edge inclusion --------------------------------
          omega[s, t] <- omega[t, s] <- 1
          G <- graph_from_adjacency_matrix(omega, mode = "undirected")
          network_1 <- igraph::components(G)

          # determine network for edge exclusion --------------------------------
          omega[s, t] <- omega[t, s] <- 0
          G <- graph_from_adjacency_matrix(omega, mode = "undirected")
          network_0 <- igraph::components(G)

          # determine the Bernoulli probability --------------------------------- #try to understand the maths of this part
          if(network_0$no == network_1$no) {
            bernoulli_probability <- edge_probability[counter]
          } else {
            indexS <- which(network_0$membership == network_0$membership[s])
            indexT <- which(network_0$membership == network_0$membership[t])
            log_lambda0 <- log(2 * cosh(sum(node_threshold[indexS])))
            log_lambda0 <- log_lambda0 +
              log(2 * cosh(sum(node_threshold[indexT])))
            indexST <- which(network_1$membership == network_1$membership[s])
            log_lambda1 <- log(2 * cosh(sum(node_threshold[indexST])))

            log_tmp1 <- -m * log_lambda1 + log(edge_probability[counter])
            log_tmp0 <- -m * log_lambda0 + log(1 - edge_probability[counter])

            bernoulli_probability <- exp(log_tmp1) /
              (exp(log_tmp1) +  exp(log_tmp0))
          }
          sampled_value <- rbinom(n = 1,
                                  size = 1,
                                  prob = bernoulli_probability)
        }
        omega[s, t] <- sampled_value
        omega[t, s] <- sampled_value
        w[counter] <- sampled_value
      }
    }
    return(w)
  }

  sample_node_threshold_dac <- function(X,
                                        person,
                                        w,
                                        node_threshold,
                                        i,
                                        a = 1,
                                        b = 1) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package `igraph' is needed for handling the topological structures.
         Please install it.", call. = FALSE)
    }

    N <- max(person)
    p <- ncol(X)
    q <- y_plus <- vector(length = N)
    current_state <- node_threshold[i]

    # extracting the likelihood parameters --------------------------------------
    for(v in 1:N) {
      # extract the data for person v -------------------------------------------
      m <- sum(person == v)
      y_plus[v] <- sum(1 / 2 + X[person == v, i] / 2)

      # build the adjacency matrix omega ----------------------------------------
      omega <- matrix(0, nrow = p, ncol = p)
      counter <- 0
      for(s in 1:(p - 1)) {
        for(t in (s + 1):p) {
          counter <- counter + 1
          omega[s, t] <- w[v, counter]
          omega[t, s] <- omega[s, t]
        }
      }

      # extract the likelihood's parameters -------------------------------------
      G <- graph_from_adjacency_matrix(omega, mode = "undirected")
      network <- igraph::components(G)
      cluster <- network$membership[i]
      if(network$csize[cluster] > 1) {
        tmp <- sum(node_threshold[network$membership == cluster]) - current_state
        q[v] <- exp(2 * tmp)
      } else {
        q[v] = 1
      }
    }

    # determine constant d of the proposal distribution -------------------------
    beta_current <- exp(2 * current_state)
    tmp1 <- (a + b) / (1 + beta_current)
    for(v in 1:N) {
      m <- sum(person == v)
      tmp1 <- tmp1 + m *  q[v] / (1 + beta_current * q[v])
    }
    tmp2 <- nrow(X) + a + b
    d <- tmp1 / (tmp2 - tmp1* beta_current)

    # proposal ------------------------------------------------------------------
    Z <- rbeta(n = 1,
               shape1 = sum(y_plus) + a,
               shape2 = nrow(X) - sum(y_plus) + b)
    proposal <- log(Z / (d * (1 - Z))) / 2

    # compute log acceptance probability ----------------------------------------
    beta_proposal <- exp(2 * proposal)
    ln_alpha <- 0
    # likelihoods ---------------------------------------------------------------
    for(v in 1:N) {
      m <- sum(person == v)
      ln_alpha <- ln_alpha + y_plus[v] * log(beta_proposal * q[v])
      ln_alpha <- ln_alpha - m * log(1 + beta_proposal * q[v])
      ln_alpha <- ln_alpha - y_plus[v] * log(beta_current * q[v])
      ln_alpha <- ln_alpha + m * log(1 + beta_current * q[v])
    }
    # priors --------------------------------------------------------------------
    ln_alpha <- ln_alpha + a * log(beta_proposal)
    ln_alpha <- ln_alpha - (a + b) * log(1 + beta_proposal)
    ln_alpha <- ln_alpha - a * log(beta_current)
    ln_alpha <- ln_alpha + (a + b) * log(1 + beta_current)

    # proposals -----------------------------------------------------------------
    ln_alpha <- ln_alpha - (sum(y_plus) + a) * log(d * beta_proposal)
    ln_alpha <- ln_alpha + (nrow(X) + a + b) * log(1 + d * beta_proposal)

    # recompute d ---------------------------------------------------------------
    #  tmp1 <- (a + b) / (1 + beta_proposal)
    #  for(v in 1:N) {
    #    m <- sum(person == v)
    #    tmp1 <- tmp1 + m *  q[v] / (1 + beta_proposal * q[v])
    #  }
    #  tmp2 <- nrow(X) + a + b
    #  d <- tmp1 / (tmp2 - tmp1* beta_proposal)


    ln_alpha <- ln_alpha + (sum(y_plus) + a) * log(d * beta_current)
    ln_alpha <- ln_alpha - (nrow(X) + a + b) * log(1 + d * beta_current)

    # Metropolis step -----------------------------------------------------------
    if(log(runif(1)) <= ln_alpha) {
      output <- proposal
    } else {
      output <- current_state
    }
    return(output)
  }
#####################################################

  ###############Online relabelling function##############################################
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
  ####indexing############################################################################
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

  #number of people
  person <- rep(1:N, each = m)
  ####### setting beta hyperparameters ######################################################
  alpha <- beta <- 1

  ####### creating a random allocation of the nodes to blocks ###############################
  az <- c(1:K) %*% rmultinom(n = p, size = 1, prob = rep(1, K))
  az <- as.vector(az)

  ######## creating a matrix for the edge probabilites (0-s) #################################
  aTheta <- matrix(0, nrow = K, ncol = K)
  for(b in 1:(K-1)) {
    for(bb in (b+1):K) {
      aTheta[b, bb] <- runif(n = 1, min = 0, max = 0.5)
      aTheta[bb, b] <- aTheta[b, bb]
    }
  }
  diag(aTheta) <-  runif(n = K, min = 0, max = 0.5)
  ######## creating node thresholds ##########################################################
  #amu <- rnorm(p)
  amu <- mu
  #####creating the individual network matrix aw and the vector version of theta##############
  counter <- 0
  aw <- matrix(0, nrow = N, ncol = choose(p, 2))
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

  #####creating containers####################################################################

  w_c <- array(data = 0, dim = c(N, choose(p, 2), iteration))
  theta_c <- list()
  z_c <- matrix(0, nrow = iteration, ncol = p)
  mu_c <- matrix(0, nrow = iteration, ncol = p)
  ##############Gibbs-sampler ################################################################

  for(iter in 1:iteration){

    progress(iter, max.value = iteration)
    #updating the vector form of Theta
    atheta <- vector(length = choose(p, 2))
    counter <- 0
    for(s in 1:(p-1)) {
      for(t in (s+1):p) {
        counter <- counter + 1
        atheta[counter] <- aTheta[az[s], az[t]]
      }
    } #using the inital blockmembership simulate the Theta edge weight vector

    #sampling for the indiviudal topology
    for(v in 1:N) {
      aw[v, ] <- sample_person_topology_dac(x = X[person == v, ], #observed
                                            w = aw[v, ], #empty
                                            node_threshold = amu, #simulated from normal
                                            edge_probability = atheta)  #for every person we simulate an individual topology based on the previously simulated aTheta and the initual values
    }

    #save the topology
    w_c[,,iter] <- aw

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

    #save the current theta matrix to the container ######
    theta_c[[iter]] <- aTheta

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

    if(iter < burn_in){
      z_c[iter, ] <- az
    } else {
      new_az <- online_relabelling(z_c[1:burn_in, ], az, K)
      z_c[iter, ] <- new_az
    }

  #   #mu conditional ########################################
  #   for(i in 1:p) {
  #     amu[i] <- sample_node_threshold_dac(X = X,
  #                                         person = person,
  #                                         w = aw,
  #                                         node_threshold = amu,
  #                                         i = i)
  #   } #sample node threshold for every node
  #
  #   mu_c[iter, ] <- amu
  #
   }

  #creating MCMC_ENT_SBM S3 class
  l <- list(theta = theta_c,
            z = z_c,
            w = w_c,
            # mu = mu_c,
            iter = iteration,
            N_Blocks = K,
            N_nodes = p,
            Sample_size = N,
            N_observation = m,
            Obs_Responses = X)

  class(l) <- "MCMC_ENT_SBM"
  return(l)
}
