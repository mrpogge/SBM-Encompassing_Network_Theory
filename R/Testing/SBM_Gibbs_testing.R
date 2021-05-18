###################
#testing
set.seed(1)
net <- SBM_sim(4, 2, N = 10, 3, 2)
net$Theta

proba <- SBM_Gibbs(net$w, net$nodes, net$K, net$N, iter = 1e5)

m <- matrix(0, 1e+5, (net$K * (net$K + 1)) / 2)
for(i in 1:1e+5){
  m[i, ] <- proba$theta[[i]][lower.tri(proba$theta[[i]], diag = T)]
}

theta_hat <- colMeans(m[1e4:1e5, ])
theta_hat
net$Theta

proba_relab <- relabelR(proba)

theta_hat_relab <- colMeans(proba_relab$theta[1e4:1e5, ])
net$Theta
theta_hat_relab


