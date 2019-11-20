library(Rtauchen)
# This example computes the transition probability matrix of the finite-state Markov chain approximation of an AR(1) process with:
# n = 5 points in the Markov chain
# ssigma = 0.02
# lambda = 0.95
# m = 3
results1 = Rtauchen(5, 0.02, 0.98, 3)
results2 = RtauchenVAR(5, 0.02^2, 0.98, 3)
print(results1-results2)

#Ejemplo Multivariado:
A=matrix(c(0.7,0.3,0.2,0.5),ncol=2)
Ni=c(3,5)
SSigma_eps=diag(2)*(0.02^2)
m=3

#Muestra la lista con la grid para cada dimensión
grid = TgridVAR(Ni, SSigma_eps, A, m)
print(grid)
for (i in 1:length(Ni)) {
  print(grid[[i]])
}
#Muestra la matriz de transión de los n1*n2*...*nI estados y muestra la suma horizontal
results3 = RtauchenVAR(Ni, SSigma_eps, A, m)
print(results3)
for (i in 1:prod(Ni)) {
  print(sum(results3[i,]))
}

#Muestra los valores de cada estado y cada dimensión
statesgrid = Tstates_grid(Ni, SSigma_eps, A, m)
print(statesgrid)

##Muestra los valores de cada estado y cada dimensión con respecto al indice que ocupan en la
#grid generada por TgridVAR
statesgridindex = Tstates_grid_index(Ni)
print(statesgridindex)


