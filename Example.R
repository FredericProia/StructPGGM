source("AlgoEstim.R")
library(mvtnorm)

# Example with :
# - q = 2 outcomes
# - p = 100 predictors
# - n = 500
# We draw w = +1/2 or -1/2 and a randomly selected row of Oyx is filled with w while the other is identically zero 
q = 2
p = 100
n = 500
w = (2*rbinom(1, 1, 1/2)-1)/2
Oyx = matrix(0, nrow=q, ncol=p)
l = rbinom(1, 1, 1/2)+1
Oyx[l,] = w

# First finite difference operator (with a small additional term for invertibility)
L = matrix(0, nrow=p, ncol=p)
for (i in 1:(p-1)){
  L[i, i] = 1
  L[i, i+1] = -1
}
L = (t(L)%*%L)/2 + 10^(-6)*diag(p)

# Residual covariance of the form (rho^|i-j|) with rho = 1/2
R = matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2)

# Simulation
Oyy = solve(R)
B = -t(Oyx)%*%R
E = rmvnorm(n, rep(0, q), sigma = R)
X = matrix(0, nrow=n, ncol=p)
Y = matrix(0, nrow=n, ncol=q)
for (i in 1:n){
  X[i,] = rnorm(p)
  Y[i,] = t(B)%*%X[i,] + E[i,]
}

# Example of PGGM estimation with penalties and no Oracle
Est = Estim(X, Y, 0.0001, 0.15, 0.01, 1.1, L, 1000, NULL, FALSE, TRUE)
Est$EstOyx
Est$EstOyy

# Same with Oracle
Est = Estim(X, Y, 0.0001, 0.15, 0.01, 1.1, L, 1000, Oyy, FALSE, TRUE)
Est$EstOyx
Est$EstOyy
