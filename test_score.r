set.seed(1);
expit <- function(x) 1/(1 + exp(-x))
source("~/COURS/SDS/stats/cours/GLMs/score_test_regression_logistique.r")

X <- cbind(1,runif(1000))
G <- runif(1000)
L <- cbind(X, G) %*% c(0, -.5, .4)
Y <- rbinom(1000, 1, expit(L))
summary(glm(Y ~ X + G - 1, family = binomial()))

two_steps2(Y, X, G)

fit <- glm(Y ~ X - 1, family = binomial)
pi <- fit$fitted.values
W <- pi*(1-pi)
Y1 <- Y - pi
XWX <- t(X) %*% (W*X)
ei <- eigen(XWX)
ei$values <- 1/sqrt(ei$values)
A <- (ei$values * t(ei$vectors)) %*% t(W * X)

Fantasio2:::logitModelScore(Y1, W, A, matrix(G, ncol = 1), 0, 0)

# -----------------------
# sans covariables

Y1 <- Y - mean(Y)
w <- mean(Y)*(1 - mean(Y))

Fantasio2:::logitModelScore_nocovar(Y1, w, matrix(G, ncol = 1), 0, 0)

