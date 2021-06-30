#LARGE SCALE TESTING
clusters <- c(2,3,5,10)
nodes <- c(3, 5, 10)
observation <- c(10, 100, 1000)

#sim_results <- list()

counter <- 0
for(clus in 1:length(clusters)){
  for(node in 1:length(nodes)){
    for(obs in 1:length(observation)){
      counter <- counter + 1

      net <- SBM_sim(p = nodes[node] * clusters[clus], K = clusters[clus], N = observation[obs], alpha = 8, beta = 3)
      estimate <- SBM_Gibbs_rl(net$w, net$nodes, net$K, net$N, iter = 1e4, burn_in = 3000)

      m <- matrix(0, 1e+4, (net$K * (net$K + 1)) / 2)
      for(i in 1:1e+4){
        m[i, ] <- estimate$theta[[i]][lower.tri(estimate$theta[[i]], diag = T)]
      }
      theta_true <- net$theta
      theta_hat <- colMeans(m[3000:1e4, ])

      z_true <- net$z
      z_MCMC <- estimate$z
      #calculating z_hat based on the MCMC run
      z_hat <- numeric(length = length(z_true))
      for(cols in 1:ncol(z_MCMC[3000:10000,])){
        z_est[cols] <- round(median(z_MCMC[, cols]))
      }

      poccc <- P_correctly_classified(z_true, z_hat, net$K)


      l <- list(network = net,
                estimate = estimate,
                theta_true = theta_true,
                theta_hat = theta_hat,
                z_true = z_true,
                z_hat = z_hat,
                correctly_classified = poccc)

      sim_results[[counter]] <- l

    }
  }
}

#extracting test results
poccc_list <- numeric(length = length(sim_results))
z_est_list <- list()
z_true_list <- list()
for(l_e in 1:length(sim_results)){

  z_true <- sim_results[[l_e]]$network$z
  z_MCMC <- sim_results[[l_e]]$estimate$z
  z_est <- numeric(length = length(z_true))
  for(cols in 1:ncol(z_MCMC[3000:10000,])){
    z_est[cols] <- round(median(z_MCMC[, cols]))
  }

  z_est_list[[l_e]] <- z_est
  z_true_list[[l_e]] <- z_true
  poccc_list[l_e] <- P_correctly_classified(z_true, z_est, sim_results[[1]]$network$nodes)

}

#checking why network 12 produced an NA (because the true network only contains 2 clusters)
z_true <- sim_results[[12]]$network$z
z_MCMC <- sim_results[[12]]$estimate$z
z_est <- numeric(length = length(z_true))
for(cols in 1:ncol(z_MCMC[3000:10000,])){
  z_est[cols] <- round(median(z_MCMC[, cols]))
}
z_12 <- online_relabelling(matrix(z_true), z_est, 3)
poccc_list[12] <- sum(z_true - z_12 == 0) / length(z_true)


sim_names <- character(length = length(sim_results))
counter <- 0
for(clus in 1:length(clusters)){
  for(node in 1:length(nodes)){
    for(obs in 1:length(observation)){
      counter <- counter + 1
      sim_names[counter] <- paste("net", clusters[clus], nodes[node]*clusters[clus], observation[obs], sep = "_")
    }
  }
}

poccc_matrix <- matrix(poccc_list, nrow = 1, ncol = 36, byrow = TRUE, dimnames = list("" , sim_names))
median(poccc_matrix[1,])
plot(poccc_matrix[1,], type = "l")


#theta values
residual_list <- list()
relab_eap_list_true <- list()
relab_eap_list_est <- list()
for(l_e in 1:length(sim_results)){
  theta_est <- sim_results[[l_e]]$theta_hat
  theta_true <- sim_results[[l_e]]$theta_true
  theta_true <- theta_true[lower.tri(theta_true, diag = TRUE)]

  residual_list[[l_e]] <- abs(sort(theta_true) - sort(theta_est))
  relab_eap_list_true[[l_e]] <- sort(theta_true)
  relab_eap_list_est[[l_e]] <- sort(theta_est)
}



#matrix format of the results of THETA
theta_MCMC_m <- matrix(0, nrow = 10000, (sim_results[[15]]$network$K * (sim_results[[15]]$network$K + 1))/2)
for(iteration in 1:length(sim_results[[15]]$estimate$theta)){

  theta_MCMC <- sim_results[[15]]$estimate$theta[[iteration]]
  theta_MCMC <- theta_MCMC[lower.tri(theta_MCMC, diag = TRUE)]
  theta_MCMC_m[iteration, ] <-theta_MCMC
}

library(LaplacesDemon)
#integrated autocorrelation time
theta_MCMC_m <- theta_MCMC_m

IAT_vec <- numeric(length = ncol(theta_MCMC_m))
for(el in 1:ncol(theta_MCMC_m)){

  IAT_vec[el] <- IAT(theta_MCMC_m[, el])
  library(tseries)
  acf(theta_MCMC_m[,el])
}

mean_residual <- numeric(length = 36)
for(l_e in 1:length(residual_list)){
  mean_residual[l_e] <- mean(unlist(residual_list[[l_e]]))
}

plot(mean_residual, type = "l")

mean(mean_residual)
sd(mean_residual)

#######################################################################################################################
#making the plot
library(gt)
library(tidyverse)

#first row
r_2_3 <- numeric()
for(r in 1:3){
  r_2_3 <- c(r_2_3, unlist(residual_list[r]))
}
r_2_3 <- mean(r_2_3)
z_2_3 <- median(poccc_matrix[, 1:3])

r_2_5 <- numeric()
for(r in 4:6){
  r_2_5 <- c(r_2_5, unlist(residual_list[r]))
}
r_2_5 <- mean(r_2_5)
z_2_5 <- median(poccc_matrix[, 4:6])

r_2_10 <- numeric()
for(r in 7:9){
  r_2_10 <- c(r_2_10, unlist(residual_list[r]))
}
r_2_10 <- mean(r_2_10)
z_2_10 <- median(poccc_matrix[, 7:9])

r_3_3 <- numeric()
for(r in 10:11){
  r_3_3 <- c(r_3_3, unlist(residual_list[r]))
}
r_3_3 <- mean(r_3_3)
z_3_3 <- median(poccc_matrix[, 10:11])

r_3_5 <- numeric()
for(r in 13:15){
  r_3_5 <- c(r_3_5, unlist(residual_list[r]))
}
r_3_5 <- mean(r_3_5)
z_3_5 <- median(poccc_matrix[, 13:15])

r_3_10 <- numeric()
for(r in 16:18){
  r_3_10 <- c(r_3_10, unlist(residual_list[r]))
}
r_3_10 <- mean(r_3_10)
z_3_10 <- median(poccc_matrix[, 16:18])

r_5_3 <- numeric()
for(r in 19:21){
  r_5_3 <- c(r_5_3, unlist(residual_list[r]))
}
r_5_3 <- mean(r_5_3)
z_5_3 <- median(poccc_matrix[, 19:21])

r_5_5 <- numeric()
for(r in 22:24){
  r_5_5 <- c(r_5_5, unlist(residual_list[r]))
}
r_5_5 <- mean(r_5_5)
z_5_5 <- median(poccc_matrix[, 22:24])

r_5_10 <- numeric()
for(r in 25:27){
  r_5_10 <- c(r_5_10, unlist(residual_list[r]))
}
r_5_10 <- mean(r_5_10)
z_5_10 <- median(poccc_matrix[, 25:27])

r_10_3 <- numeric()
for(r in 28:30){
  r_10_3 <- c(r_10_3, unlist(residual_list[r]))
}
r_10_3 <- mean(r_10_3)
z_10_3 <- median(poccc_matrix[, 28:30])

r_10_5 <- numeric()
for(r in 31:33){
  r_10_5 <- c(r_10_5, unlist(residual_list[r]))
}
r_10_5 <- mean(r_10_5)
z_10_5 <- median(poccc_matrix[, 31:33])

r_10_10 <- numeric()
for(r in 34:36){
  r_10_10 <- c(r_10_10, unlist(residual_list[r]))
}
r_10_10 <- mean(r_10_10)
z_10_10 <- median(poccc_matrix[, 34:36])

#creating the dataframe for the table
vec <- c(r_2_3, z_2_3,
         r_2_5, z_2_5,
         r_2_10, z_2_10,
         r_3_3, z_3_3,
         r_3_5, z_3_5,
         r_3_10, z_3_10,
         r_5_3, z_5_3,
         r_5_5, z_5_5,
         r_5_10, z_5_10,
         r_10_3, z_10_3,
         r_10_5, z_10_5,
         r_10_10, z_10_10)
table1 <- matrix(vec, nrow = 12, ncol = 6, byrow = TRUE)
table1 <- cbind(c(2,3,5,10), table1)
table1 <- round(table1, digits = 2)
table1 <- as.data.frame(table1)
colnames(table1) <- c( "rwname", "3RES", "3CC", "5RES", "5CC", "10RES", "10CC")
row.names(table1) <- c("2_10", "2_100","2_1000","3_10", "3_100","3_1000","5_10","5_100","5_1000","10_10", "10_100", "10_1000")
#making the table itself
# gt_table1 <- gt(table1, rowname_col = "rwname") %>%
#   tab_stubhead(label = "Number of Blocks (k)") %>%
#   tab_spanner(label = "k x 3 nodes", columns = c("3RES", "3CC")) %>%
#   tab_spanner(label = "k x 5 nodes", columns = c("5RES", "5CC")) %>%
#   tab_spanner(label = "k x 10 nodes", columns = c("10RES", "10CC")) %>%
#   cols_label("3RES" = "Residual",
#              "3CC" = "Correctly classified",
#              "5RES" = "Residual",
#              "5CC" = "Correctly classified",
#              "10RES" = "Residual",
#              "10CC" = "Correctly classified") %>%
#   tab_footnote(footnote = "The network with 3 clusters, 9 nodes and 1000 observation
#                was discarded because the random network generation process created only 2 blocks instead of 3.",
#                locations = cells_body(
#                  columns = c("3RES", "3CC"),
#                  rows = "rwname" == 3)
#                )

library(stargazer)
sg_table1 <- stargazer(table1, summary = FALSE, rownames = FALSE,
                        title = "Table 1: Results of the simulation study for the ordinary SBM"
                        )

##########################################################################################################
#MCMC diagnostics with CODA

library(coda)

#making a list of the theta MCMC samples across all simulation
MCMC_list <- list()
for(l_e in 1:36){

  theta_MCMC_m <- matrix(0, nrow = 10000, (sim_results[[l_e]]$network$K * (sim_results[[l_e]]$network$K + 1))/2)
  for(iteration in 1:length(sim_results[[l_e]]$estimate$theta)){

    theta_MCMC <- sim_results[[l_e]]$estimate$theta[[iteration]]
    theta_MCMC <- theta_MCMC[lower.tri(theta_MCMC, diag = TRUE)]
    theta_MCMC_m[iteration, ] <- theta_MCMC
  }

  MCMC_list[[l_e]] <- theta_MCMC_m
}

geweke_list <- list()
heidelberg_list <- list()
autocor_list <- list()
effective_list <- list()

for(l_e in 1:36){

   #geweke_list[[l_e]] <- geweke.diag(as.mcmc(MCMC_list[[l_e]][3000:10000]))
  #heidelberg_list[[l_e]] <- heidel.diag(as.mcmc(MCMC_list[[l_e]]))
  autocor_list[[l_e]] <- autocorr.diag(as.mcmc(MCMC_list[[l_e]]))
  #effective_list[[l_e]] <- effectiveSize(as.mcmc(MCMC_list[[l_e]]))


}

conv_logic <- logical(length = 36)
for(l_e in 1:36){
  conv_logic[l_e] <- any(abs(geweke_list[[l_e]]$z) > 2.58)
}

sum(conv_logic) / 36

#number of simulation where the heidelberg test passed for all of the theta values
heidelberg_list[[29]]

station_logic <- logical(length = 36)
for(l_e in 1:36){
  station_logic[l_e] <- any(heidelberg_list[[l_e]][ ,1] == 0)
}

sum(station_logic)


autocorr.plot(MCMC_list[[18]][3000:10000,1])

  per_obs <- list(ten = numeric(length = 12), hundred = numeric(length = 12), thousand = numeric(length = 12))
  counter_1 <- 0
  counter_2 <- 0
  counter_3 <- 0
  for(l_e in 1:36){
    if(l_e %% 3 == 1){
      counter_1 <- counter_1 + 1
      per_obs[[1]][counter_1] <- median(residual_list[[l_e]])
    } else if(l_e %% 3 == 2){
      counter_2 <- counter_2 + 1
      per_obs[[2]][counter_2] <- median(residual_list[[l_e]])
    } else if(l_e %% 3 == 0){
      counter_3 <- counter_3 + 1
      per_obs[[3]][counter_3] <- median(residual_list[[l_e]])
    }
  }

mad(per_obs[[1]])
mad(per_obs[[2]])
mad(per_obs[[3]])

MCMC_EAP <- list()
for(l_e in 1:length(MCMC_list)){

  MCMC_EAP[[l_e]] <- colMeans(MCMC_list[[l_e]][3000:10000, ])

}

sim_EAP <- list()
for(l_e in 1:length(sim_results)){

  temp_d <- sim_results[[l_e]]$network$theta
  sim_EAP[[l_e]] <- temp_d[lower.tri(temp_d, diag = TRUE)]

}

bias_df <- cbind(unlist(MCMC_EAP), unlist(sim_EAP))


abline(0,1)
cor(bias_df[,1], bias_df[, 2])

relab_eap <- cbind(unlist(relab_eap_list_est), unlist(relab_eap_list_true))

layout(1:2)
plot(bias_df[,1], bias_df[, 2], ylim = c(0,1), xlim = c(0,1), ylab = "Simulated Theta", xlab = "Estimated Theta")
abline(0,1)
plot(relab_eap[,1], relab_eap[, 2], ylim = c(0,1), xlim = c(0,1), ylab = "Simulated Theta", xlab = "Estimated Theta")
abline(0,1)


cor(relab_eap[, 1], relab_eap[, 2])

layout(1:2)
plot(MCMC_list[[18]][1:3000,1], type = "l", ylim = c(0.280, 0.305), xlim = c(1, 3000), ylab = "Theta", xlab = "Number of iterations")
plot(MCMC_list[[18]][3000:10000,1], type = "l", ylim = c(0.280, 0.305), xlim = c(1, 7000), ylab = "Theta", xlab = "Number of iterations")


layout(1:2)
acf(MCMC_list[[18]][1:3000 ,1], main = "", ylim = c(0,1), xlim = c(0, 35))
acf(MCMC_list[[18]][3000:10000,1], main = "")


###################################################################################
#Creating the new table

table_m <- numeric(length = 36)

for(i in 1:36){

  table_m[i] <- mean(residual_list[[i]])

}
vec <- numeric()
for(i in 1:36){

  vec <- c(vec, table_m[i], poccc_list[i])

}

table1 <- matrix(vec, nrow = 12, ncol = 6, byrow = TRUE)
table1 <- cbind(c(2,3,5,10), table1)
table1 <- round(table1, digits = 2)
table1 <- as.data.frame(table1)
colnames(table1) <- c( "rwname", "3RES", "3CC", "5RES", "5CC", "10RES", "10CC")
row.names(table1) <- c("2_10", "2_100","2_1000","3_10", "3_100","3_1000","5_10","5_100","5_1000","10_10", "10_100", "10_1000")

library(stargazer)
sg_table1 <- stargazer(table1, summary = FALSE, rownames = FALSE,
                       title = "Table 1: Results of the simulation study for the ordinary SBM"
)

mean(table_m)
