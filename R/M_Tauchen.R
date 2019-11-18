
#Discretizes a vector VAR(1) process of size M according to Tauchen (1986)

#The package includes the functions TgridVAR, RtauchenVAR 

#TgridVAR generates the grid for the aproximating discrete-valued 
#vector Markov Chain of the VAR(1) process

#Inputs: 
#Ni: Vector of size M with the number of grid points for each component of the
#autoregresive vector y_t
#SSigma_eps: Diagonal matrix of variances of the white noise
#A: Matrix of persistences of the process
#m: Number of standard deviations around the mean that the grid measures

#Output:
#List of M vectors with the sizes specified in Ni. Each vector contains 
#the posible states of each component of the autoregresive vector y_t

TgridVAR=function(Ni, SSigma_eps, A, m){
#Defines the length of the random vector of the process
  M = length(Ni)
#Initializes values for the grid, variances matrix of the vector, and lenght of step for each dimension
  grid = vector(mode = "list", length = M)
  SSigma_y = matrix(0,M,M) 
  w = rep(0,M)
#Vectorization of white noise variances matrix
  vec_SSigma_eps = matrix(c(SSigma_eps))
#Computes the vectorized unconditional variance of the random vector 
  if(M==1){
    vec_SSigma_y = solve(1 - kronecker(A,A)) %*% vec_SSigma_eps
  }
  else{
    vec_SSigma_y = solve(diag(M^2) - kronecker(A,A)) %*% vec_SSigma_eps
  }
#Converts the result to a matrix of variances
  SSigma_y = matrix(vec_SSigma_y,M,M)
#Obtains the standard deviations asociated with the matrix
  ssigma_y = sqrt(diag(SSigma_y))
#Sets the values for the lenght of steps of the gird accordingly for each component
  w = 2*m*ssigma_y/(Ni-1)
  
#Creates a list in which each element is the grid for each component of y_t
  for(j in 1:M){
    aux_vec_yi = -m*ssigma_y[j] + c(0:(Ni[j]-1))*w[j]
    grid[[j]]= aux_vec_yi
  }
#Returns the list with the grid for each dimension
  return(grid)
}

Tgridstates=function(Ni, SSigma_eps, A, m){
  
}

#RTauchenVAR generates the transition matrix for all posible states. The states 
#are in a vector of the lenght specified in Tauchen(1986).

Rtauchen_vec=function(Ni, SSigma_eps, A, m){
  
  grid = TgridVAR(Ni, SSigma_eps, A, m)
  N_star = prod(Ni)
  M = length(Ni)

  grid_index = vector(mode = "list", length = M)
  for(j in 1:M){
    aux_vec_l = c(1:(Ni[j]))
    grid_index[[j]]= aux_vec_l
  }
  
  states_grid = expand.grid(grid)
  rownames(states_grid) = paste0("state", seq_len(nrow(states_grid)))
  colnames(states_grid) = paste0("y", seq_len(ncol(states_grid)))
  states_grid = data.matrix(states_grid)
  
  states_grid_index = expand.grid(grid_index)
  rownames(states_grid_index) = paste0("state", seq_len(nrow(states_grid_index)))
  colnames(states_grid_index) = paste0("y", seq_len(ncol(states_grid_index)))
  states_grid_index = data.matrix(states_grid_index)
  
  
  
  if(M==1){
    vec_SSigma_eps = matrix(c(SSigma_eps))
    vec_SSigma_y = solve(1 - kronecker(A,A)) %*% vec_SSigma_eps
  }else{
    vec_SSigma_eps = matrix(c(SSigma_eps))
    vec_SSigma_y = solve(diag(M^2) - kronecker(A,A)) %*% vec_SSigma_eps
  }
  SSigma_y = matrix(vec_SSigma_y,M,M)
  SSigma_eps = matrix(vec_SSigma_eps,M,M)
  ssigma_y = sqrt(diag(SSigma_y))
  ssigma_eps = sqrt(diag(SSigma_eps))
  w = 2*m*ssigma_y/(Ni-1)
  
  h = vector(mode = "list", length = N_star)
  
  for(j in 1:N_star){
    h[[j]]= grid
  }
  
  for(j in 1:N_star){
    y=states_grid[j,]
    mu=A%*%y
    for(i in 1:M){
      if(Ni[i] > 1){
        for(l in 1:Ni[i]){
          if(l == 1){
            h[[j]][[i]][l] = pnorm((grid[[i]][l] - mu[i] +(w[i]/2))/ssigma_eps[i])
          }else if(l == Ni[i]){
            h[[j]][[i]][l] = 1- pnorm((grid[[i]][l] - mu[i]-(w[i]/2))/ssigma_eps[i])
          }else{
            h[[j]][[i]][l] = pnorm((grid[[i]][l] - mu[i] +(w[i]/2))/ssigma_eps[i])-pnorm((grid[[i]][l]-mu[i]-(w[i]/2))/ssigma_eps[i])
          }
        }
      }
      else{
        h[[j]][[i]][i] = 1
      }
    }
  }
  
  P = array(0, dim = c(N_star,N_star))
  aux_prod = rep(0,M)
  
  for(j in 1:N_star){
    aux_matrix = data.matrix(expand.grid(h[[j]]))
    for(k in 1:N_star){
      P[j,k]=prod(aux_matrix[k,])
      }
    }
  return(P)
}

#Ejemplo Multivariado:
# A=matrix(c(0.7,0.3,0.2,0.5),ncol=2)
# Ni=c(5,5)
# SSigma_eps=diag(2)*(0.02^2)
# m=3

#Ejemplo Univariado:
# A=0.98
# Ni=5
# SSigma_eps=0.02^2
# m=3
