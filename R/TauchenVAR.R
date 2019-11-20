
#Discretizes a vector VAR(1) process of size M according to Tauchen (1986)

#The package includes the functions TgridVAR, Tstates_grid, Tstates_grid_index y RtauchenVAR 

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
#the posible states of each dimension of the autoregresive vector y_t

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

#Tstates_grid returns an N_star x M matrix, in which each row is a posible state
#of the MArkov chain, and each row is the value that the dimension i of the vector 
#takes in that state

Tstates_grid=function(Ni, SSigma_eps, A, m){
  grid = TgridVAR(Ni, SSigma_eps, A, m)
  
  states_grid = expand.grid(grid)
  rownames(states_grid) = paste0("state", seq_len(nrow(states_grid)))
  colnames(states_grid) = paste0("y", seq_len(ncol(states_grid)))
  states_grid = data.matrix(states_grid)
  
  return(states_grid)
}

#Tstates_grid_index returns an N_star x M matrix, in which each row is a posible state
#of the MArkov chain, and each row is the index of the grid for state i that
#the dimension i of the vector takes in that state

Tstates_grid_index=function(Ni){
  M = length(Ni)
  grid_index = vector(mode = "list", length = M)
  for(j in 1:M){
    aux_vec_l = c(1:(Ni[j]))
    grid_index[[j]]= aux_vec_l
  }
  
  states_grid_index = expand.grid(grid_index)
  rownames(states_grid_index) = paste0("state", seq_len(nrow(states_grid_index)))
  colnames(states_grid_index) = paste0("y", seq_len(ncol(states_grid_index)))
  states_grid_index = data.matrix(states_grid_index)
  
  return(states_grid_index)
}

#RTauchenVAR generates the transition matrix for all posible states. The states 
#are in a vector of the lenght N_star, as specified in Tauchen(1986).

RtauchenVAR=function(Ni, SSigma_eps, A, m){
#Takes the lenght of the vector
  M = length(Ni)
#Takes the number of states of the MArkov chain
  N_star = prod(Ni)
#Obtains the grid for each dimension
  grid = TgridVAR(Ni, SSigma_eps, A, m)
  states_grid = Tstates_grid(Ni, SSigma_eps, A, m)
  states_grid_index = Tstates_grid_index(Ni)
#obtains the variances and places them accordingly
  if(M==1){
    vec_SSigma_eps = matrix(c(SSigma_eps))
    vec_SSigma_y = solve(1 - kronecker(A,A)) %*% vec_SSigma_eps
  }else{
    vec_SSigma_eps = matrix(c(SSigma_eps))
    vec_SSigma_y = solve(diag(M^2) - kronecker(A,A)) %*% vec_SSigma_eps
  }
  
  SSigma_y = matrix(vec_SSigma_y,M,M)
  ssigma_y = sqrt(diag(SSigma_y))
  
  SSigma_eps = matrix(vec_SSigma_eps,M,M)
  ssigma_eps = sqrt(diag(SSigma_eps))

#Sets the length of the step for each dimension
  w = 2*m*ssigma_y/(Ni-1)
#Initializes the list to storage the transition probabilities for each diemnsion
  h = vector(mode = "list", length = N_star)
  for(j in 1:N_star){
    h[[j]]= grid
  }
#Computes the transition probabilities for each dimension of each state, as in Tauchen(1986).
#h[[j]][[i]][l] is the probability that, starting in state j, the dimension i of the vector 
#takes the value l of the grid for that dimension, on the next state.
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

#Initializes the transition matrix for the Markov chain  
  P = array(0, dim = c(N_star,N_star))

#Generates the transition probabilities form state j to state k as in Tauchen(1986).
#For each state, all the probabilities of transition for each dimension are stored 
#in  data.matrix(expand.grid(h[[j]])). Then, byu independence, the probability of changing
#state is the product of the probability of changing value in each dimension.
  
  for(j in 1:N_star){
    aux_matrix = data.matrix(expand.grid(h[[j]]))
    for(k in 1:N_star){
      P[j,k]=prod(aux_matrix[k,])
      }
    }
  return(P)
}


