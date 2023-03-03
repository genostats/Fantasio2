require(Fantasio2)

N <- 100
P <- 1e4

Y <- sample(0:1, N, TRUE)
covar <- matrix(runif(N), N)

H <- matrix( runif(N*P), nrow = N)

st1 <- system.time(A <- Fantasio2:::glm.HBD.0(Y, cbind(1,covar), H))

st2 <- system.time( {
B <- matrix(0, nrow = P, ncol = 2)
for(i in 1:P) {
  fit <- summary(glm(Y ~ covar + H[,i], family = binomial()))
  B[i,] <- fit$coefficients[3,1:2]
}
})

plot( B[,1], A$beta )
