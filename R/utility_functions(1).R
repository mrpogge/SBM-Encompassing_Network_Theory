topology_clustering <- function(N, p, w) {
  clustering <- matrix(0, nrow = N, ncol = p)
  clustersize <- list()
  for(v in 1:N) {
    # Build the topology of person v ------------------------------------------
    topology <- matrix(0, nrow = p, ncol = p)
    counter <- 0
    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        counter <- counter + 1
        topology[s, t] <- w[v, counter]
        topology[t, s] <- topology[s, t]
      }
    }
    
    # Extract the clusters ----------------------------------------------------
    G <- graph_from_adjacency_matrix(topology, mode = "undirected")
    tmp <- igraph::components(G)
    clustering[v, ] <- tmp$membership
    clustersize[[v]] <- tmp$csize
  }
  return(list(clustering = clustering, clustersize = clustersize))
}

topology_clustering_parallel <- function(N, p, w, numCores, mc.preschedule = TRUE) {
  #The following function is parallelized -------------------------------------
  topology_v <- function(v) {
    # Build the topology of person v ------------------------------------------
    topology <- matrix(0, nrow = p, ncol = p)
    counter <- 0
    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        counter <- counter + 1
        topology[s, t] <- w[v, counter]
        topology[t, s] <- topology[s, t]
      }
    }
    
    # Extract the clusters ----------------------------------------------------
    G <- graph_from_adjacency_matrix(topology, mode = "undirected")
    tmp <- igraph::components(G) 
    return(list(clustering = tmp$membership, 
                clustersize = tmp$csize, 
                person = v)) 
  }
  
  clustering <- matrix(0, nrow = N, ncol = p)
  clustersize <- list()
  topologies <- mclapply(X = 1:N, 
                  FUN = topology_v,
                  mc.cores = numCores,
                  mc.preschedule = mc.preschedule)
  # Assign output of mclapply in the correct order ----------------------------
  for(v in 1:N) {
    clustering[topologies[[v]]$person, ] <- topologies[[v]]$clustering
    clustersize[[topologies[[v]]$person]] <- topologies[[v]]$clustersize
  }
  return(list(clustering = clustering, clustersize = clustersize))
}
