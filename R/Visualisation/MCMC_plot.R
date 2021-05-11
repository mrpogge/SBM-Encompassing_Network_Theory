#first run the algorithm and save the mcmc output

plot.MCMC_SBM <- function(obj, type = c("line_theta", "density_theta", "density_z")){

type <- match.arg(type, choices = c("line_theta", "density_theta", "density_z"))

#selecting theta from the object
  theta_c <- obj$theta
  m <- matrix(0, 1000, (obj$N_Blocks * (obj$N_Blocks+1))/2)

  for(i in 1:1000){
    m[i, ] <- theta_c[[i]][lower.tri(theta_c[[i]], diag = T)]
  }

#auto.layout function from ENMwizard
#auto.layout = function(n, layout=T){
#    ### figure out how many rows
#    sq = sqrt(n)
#    rws = round(sq)
#
#    #### if it's a perfect square, fill the matrix
#    if (sqrt(n) == round(sqrt(n))){
#      m = matrix(1:4, nrow=sq, byrow=T)
#    } else {
#
#      #### repeat twice the numbers that fit nicely
#      topNum = trunc(n/rws)*rws
#      numbs = 1:topNum
#      if (topNum==n){
#        m = matrix(numbs, nrow=rws, byrow=T)
#      } else {
#        #### get the rest figured out
#        rem = n-topNum  ### remaining numbers
#        rest = (topNum+1):n
#        cols = topNum/rws
#        rest = c(rep(0, times= cols-length(rest)), rest, rep(0, times= cols-length(rest)))
#        m = matrix(c(numbs, rest), nrow=rws+1, byrow=T)
#      }
#    }
#
#    if (layout){
#      layout(m)
#    } else {
#      m
#    }
#  }
#plotting according to the type

  if(type == "line_theta"){

  #creating theta line graph from MCMC output
    #layout(auto.layout(ncol(m), layout = FALSE))
    par(mfrow = c(3,5))
    for(i in 1:ncol(m)){
      plot(m[,i], type = "l")
    }

  } else if(type == "density_theta"){
  #creating theta histogram from MCMC output
    #auto.layout(ncol(m))
    par(mfrow = c(3,5))
    for(i in 1:ncol(m)){
      hist(m[,i])
    }
  } else if(type == "density_z"){
  #creating z histogram from MCMC output
    #auto.layout(obj$N_nodes)
    par(mfrow = c(3,5))
    for(i in 1:obj$N_nodes){
      hist(obj$z[,i])
    }
  }


}

############################# from ENMwizard
auto.layout = function(n, layout=T){
  ### figure out how many rows
  sq = sqrt(n)
  rws = round(sq)

  #### if it's a perfect square, fill the matrix
  if (sqrt(n) == round(sqrt(n))){
    numbs = sort(rep(1:n, times=2))
    m = matrix(numbs, nrow=sq, byrow=T)
  } else {

    #### repeat twice the numbers that fit nicely
    topNum = trunc(n/rws)*rws
    numbs = sort(rep(1:topNum, times=2))
    if (topNum==n){
      m = matrix(numbs, nrow=rws, byrow=T)
    } else {
      #### get the rest figured out
      rem = n-topNum  ### remaining numbers
      rest = sort(rep((topNum+1):n, times=2))
      cols = (topNum/rws)*2
      rest = c(rep(0, times=(cols-length(rest))/2), rest, rep(0, times=(cols-length(rest))/2))
      m = matrix(c(numbs, rest), nrow=rws+1, byrow=T)
    }
  }

  if (layout){
    layout(m)
  } else {
    m
  }
}



theta_c <- res$theta
m <- matrix(0, 1000, (K * (K+1))/2)

for(i in 1:1000){
 m[i, ] <- theta_c[[i]][lower.tri(theta_c[[i]], diag = T)]
}



par(mfrow = c(3,5))
for(i in 1:ncol(m)){
  plot(m[,i], type = "l")
}




par(mfrow = c(3,5))
for(i in 1:ncol(m)){

  hist(m[,i])
}



par(mfrow = c(4,5))
for(i in 1:p){

  hist(z_c[,i])
}


colMeans(m[400:1000,])
Theta_true[lower.tri(Theta_true, diag = TRUE)]

m_hat <- matrix(0, K, K)
m_hat[lower.tri(m_hat, diag = TRUE)] <- colMeans(m[400:1000,])
m_hat
Theta_true

colMeans(m_atheta[200:1000,])
Theta_true[lower.tri(Theta_true, diag = TRUE)]

m_atheta_hat <- matrix(0, K, K)
m_hat[lower.tri(m_atheta_hat, diag = TRUE)] <- colMeans(m_atheta[400:1000,])
m_hat
Theta_true
