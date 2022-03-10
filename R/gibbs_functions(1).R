#' Draw from the posterior distribution of an idiographic topology. 
#'
#' The function \code{sample_person_topology_dac} is used to sample from the 
#' full-conditional posterior distribution of the network structure for a
#' single person. It is used by \code{TOBEDETERMINED} to summarize an 
#' idiographic topology and by \code{MORETOBEDETERMINED} to estimate Divide and
#' Color models with the Gibbs sampler. 
#'
#' The function \code{sample_person_topology_dac} samples from the 
#' full-conditional posterior distribution of a person's topology given the 
#' topology's current state, the states of the node thresholds and edge 
#' probabilities. The stipulated prior distribution on the topology is an
#' simple random graph model (i.e., an Erd\H{o}s-R\'{e}nyi model or Bernoulli 
#' model).
#'
#' @param x An \code{m} by \code{p} matrix containing effect coded
#'   responses (i.e., coded \code{-1, +1}) for \code{m} independent 
#'   observations on the \code{p} variables of the network or graph.
#'
#' @param w An \code{1} by \code{p \choose 2} matrix containing binary 
#'   variables indicating absence (coded \code{0}) and presence 
#'   (coded \code{1}) of edges in the idiographic network.
#'
#' @param node_threshold A vector of length \code{p} containing the thresholds
#'   of the network's nodes. 
#'   
#' @param edge_probabilities A \code{p} by \code{p} matrix containing the prior
#'   probability of edge inclusion. 
#'   
#' @return An \code{1} by \code{p \choose 2} matrix containing updated binary 
#'   variables indicating absence (coded \code{0}) and presence 
#'   (coded \code{1}) of edges in the idiographic network.
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
  
  # update the edge states ----------------------------------------------------
  counter <- 0
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
        
        # determine the Bernoulli probability ---------------------------------
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
          max.log <- max(c(log_tmp0, log_tmp1))
          
          bernoulli_probability <- exp(log_tmp1 - max.log) / 
            (exp(log_tmp1 - max.log) +  exp(log_tmp0 - max.log))
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


#' Draw from the posterior distribution of an idiographic topology. 
#'
#' The function \code{sample_person_topology_dac_delta} is used to sample from the 
#' full-conditional posterior distribution of the network structure for a
#' single person. It is used by \code{TOBEDETERMINED} to summarize an 
#' idiographic topology and by \code{MORETOBEDETERMINED} to estimate Divide and
#' Color models with the Gibbs sampler. 
#'
#' The function \code{sample_person_topology_dac} samples from the 
#' full-conditional posterior distribution of a person's topology given the 
#' topology's current state, the states of the node thresholds and edge 
#' probabilities. The stipulated prior distribution on the topology is an
#' simple random graph model (i.e., an Erd\H{o}s-R\'{e}nyi model or Bernoulli 
#' model).
#'
#' @param delta An \code{1} by \code{p \choose 2} matrix containing indicators 
#'   for nodes being in same state across waves.
#'
#' @param w An \code{1} by \code{p \choose 2} matrix containing binary 
#'   variables indicating absence (coded \code{0}) and presence 
#'   (coded \code{1}) of edges in the idiographic network.
#'
#' @param node_threshold A vector of length \code{p} containing the thresholds
#'   of the network's nodes. 
#'   
#' @param edge_probabilities A \code{p} by \code{p} matrix containing the prior
#'   probability of edge inclusion. 
#'   
#' @param m A positive integer, the number of observations on the topology. 
#'   
#' @return An \code{1} by \code{p \choose 2} matrix containing updated binary 
#'   variables indicating absence (coded \code{0}) and presence 
#'   (coded \code{1}) of edges in the idiographic network.
sample_person_topology_dac_delta <- function(delta,
                                       w, 
                                       node_threshold,
                                       edge_probability,
                                       m) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph' is needed for handling the topological structures. 
         Please install it.", call. = FALSE)
  }
  p <- length(node_threshold)
  
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
  
  # update the edge states ----------------------------------------------------
  counter <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      counter <- counter + 1
      
      # node states do not match (delta[s,t] = 0) -----------------------------
      if(delta[counter] == FALSE) {
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
        
        # determine the Bernoulli probability ---------------------------------
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
          max.log <- max(c(log_tmp0, log_tmp1))
          
          bernoulli_probability <- exp(log_tmp1 - max.log) / 
            (exp(log_tmp1 - max.log) +  exp(log_tmp0 - max.log))
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




#' Draw from the posterior distribution of a node threshold. 
#'
#' The function \code{sample_node_threshold_dac_old} is used to sample from the 
#' full-conditional posterior distribution of the node threshold. It is used 
#' for backward compatibility and replaced by \code{sample_node_threshold_dac}.
#'
#' @param X An \code{N} by \code{p} matrix containing effect coded
#'   responses (i.e., coded \code{-1, +1}) for \code{m[v]} independent 
#'   observations for each of \code{N} persons on the \code{p} variables 
#'   of the network or graph.
#'   
#' @param person A vector of length \code{O} indicating which row in \code{X} 
#'   belongs to a person \code{1}, \code{2}, ..., \code{N}. 
#'
#' @param w An \code{N} by \code{p \choose 2} matrix containing binary 
#'   variables indicating absence (coded \code{0}) and presence 
#'   (coded \code{1}) of edges in the idiographic network.
#'
#' @param node_threshold A vector of length \code{p} containing the thresholds
#'   of the network's nodes. 
#'   
#' @param i A positive integer indicating which threshold to update.
#'   
#' @param a,b The hyperparameters of the Beta-Prime distribution.
#' 
#' @return The outcome value of the Metropolis algorithm. 
sample_node_threshold_dac_old <- function(X,
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


#' Draw from the posterior distribution of a node threshold. 
#'
#' The function \code{sample_node_threshold_dac} is used to sample from the 
#' full-conditional posterior distribution of the node threshold. It is used by 
#' \code{TOBEDETERMINED} to estimate Divide and Color models with the Gibbs 
#' sampler. 
#'
#' The function \code{sample_node_threshold_dac} uses a Metropolis algorithm to 
#' update the current of a node's threshold in the Markov chain given the 
#' states of the other node thresholds and edges. The stipulated prior 
#' distribution on the threshold is a generalized Beta-Prime distribution.
#'
#' @param nr_obs_person An \code{1} by \code{N} matrix containing the number of
#'   observations per person.
#'   
#' @param y_plus The number of "+1" observations per node. That is, the sum of 
#'   (0, 1)-coded node states.
#'
#' @param p number of nodes in the network.
#'    
#' @param clustering A \code{N} by \code{p} matrix containing the cluster 
#'   assignments of the \code{p} nodes for each of the \code{N} topologies. 
#'   Computed by \code{topology_clustering}.
#'   
#' @param clustersize A list of length \code{N} containing the cluster sizes
#'   for each of the \code{N} topologies. Computed by 
#'   \code{topology_clustering}.
#'
#' @param node_threshold A vector of length \code{p} containing the current 
#'   states of the node thresholds.
#'   
#' @param node A positive integer indicating which threshold to update.
#'   
#' @param a,b The hyperparameters of the Beta-Prime distribution.
#' 
#' @return The outcome value of the Metropolis algorithm. 
sample_node_threshold_dac <- function(nr_obs_person,
                                      y_plus,
                                      p,
                                      clustering,
                                      clustersize,
                                      node_threshold,
                                      node,
                                      a = 1,
                                      b = 1) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph' is needed for handling the topological structures. 
         Please install it.", call. = FALSE)
  }
  
  total_obs <- sum(nr_obs_person)
  N <- length(nr_obs_person)
  q <- vector(length = N)
  current_state <- node_threshold[node]
  
  # extracting the likelihood parameters --------------------------------------
  for(v in 1:N) {
    cluster <- clustering[v, node]
    if(clustersize[[v]][cluster] > 1) {
      tmp <- sum(node_threshold[clustering[v, ] == cluster]) - current_state
      q[v] <- exp(2 * tmp)
    } else {
      q[v] = 1
    }
  }
  
  # determine constant d of the proposal distribution -------------------------
  beta_current <- exp(2 * current_state)
  tmp1 <- (a + b) / (1 + beta_current)
  for(v in 1:N) {
    tmp1 <- tmp1 + nr_obs_person[v] *  q[v] / (1 + beta_current * q[v])
  }
  tmp2 <- total_obs + a + b
  d <- tmp1 / (tmp2 - tmp1* beta_current)
  
  # proposal ------------------------------------------------------------------
  Z <- rbeta(n = 1, 
             shape1 = sum(y_plus[, node]) + a, 
             shape2 = total_obs - sum(y_plus[, node]) + b)
  proposal <- log(Z / (d * (1 - Z))) / 2
  
  # compute log acceptance probability ----------------------------------------
  beta_proposal <- exp(2 * proposal) 
  ln_alpha <- 0
  # likelihoods ---------------------------------------------------------------
  for(v in 1:N) {
    ln_alpha <- ln_alpha + y_plus[v, node] * log(beta_proposal * q[v])
    ln_alpha <- ln_alpha - nr_obs_person[v] * log(1 + beta_proposal * q[v])
    ln_alpha <- ln_alpha - y_plus[v, node] * log(beta_current * q[v])
    ln_alpha <- ln_alpha + nr_obs_person[v] * log(1 + beta_current * q[v])
  }
  # priors --------------------------------------------------------------------
  ln_alpha <- ln_alpha + a * log(beta_proposal)
  ln_alpha <- ln_alpha - (a + b) * log(1 + beta_proposal)
  ln_alpha <- ln_alpha - a * log(beta_current)
  ln_alpha <- ln_alpha + (a + b) * log(1 + beta_current)
  
  # proposals -----------------------------------------------------------------
  ln_alpha <- ln_alpha - (sum(y_plus[, node]) + a) * log(d * beta_proposal)
  ln_alpha <- ln_alpha + (total_obs + a + b) * log(1 + d * beta_proposal)
  
  ln_alpha <- ln_alpha + (sum(y_plus[, node]) + a) * log(d * beta_current)
  ln_alpha <- ln_alpha - (total_obs + a + b) * log(1 + d * beta_current)
  
  # Metropolis step -----------------------------------------------------------
  if(log(runif(1)) <= ln_alpha) {
    output <- proposal
  } else {
    output <- current_state
  }
  return(output)
}

