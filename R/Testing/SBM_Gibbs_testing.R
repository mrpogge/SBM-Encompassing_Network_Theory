#LARGE SCALE TESTING
clusters <- c(2,3,5,10)
nodes <- c(3, 5, 10)
observation <- c(10, 100, 1000)

sim_results <- list()

counter <- 0
for(clus in 1:length(clusters)){
  for(node in 1:length(nodes)){
    for(obs in 1:length(observation)){
      counter <- counter + 1

      net <- SBM_sim(p = nodes[node] * clusters[clus], K = clusters[clus], N = observation[obs], alpha = 8, beta = 3)
      estimate <- SBM_Gibbs(net$w, net$nodes, net$K, net$N, iter = 1e4)

      m <- matrix(0, 1e+4, (net$K * (net$K + 1)) / 2)
      for(i in 1:1e+4){
        m[i, ] <- estimate$theta[[i]][lower.tri(estimate$theta[[i]], diag = T)]
      }

      theta_hat <- colMeans(m[1e3:1e4, ])

      estimate_relab <- relabelR(estimate)

      theta_hat_relab <- colMeans(estimate_relab$theta[1e3:1e4, ])

      z_hat <- numeric()
      for(zs in 1:net$K){
        z_hat[zs] <- median(estimate_relab$z[, zs])
      }


      l <- list(network = net,
                estimate = estimate,
                theta_hat = theta_hat,
                estimate_relabelled = estimate_relab,
                theta_hat_relabelled = theta_hat_relab,
                theta_true = net$theta,
                z_hat = z_hat,
                z_true = net$z)

      sim_results[[counter]] <- l

    }
  }
}







