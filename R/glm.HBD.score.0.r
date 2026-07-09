# Y = vecteur de 0 et 1
# covar.matrix = matrice de covariables [doit contenir un intercept]
# H = matrice des pHBD ou des FLOD...
# variance = ne sera pris en compte que si covar.matrix est réduit à un intercept
#            c'est le vecteur des variances des scores si elles ont déjà été calculées
#            -> utile pour des permutations car la variance est toujours la même
glm.HBD.score.0 <- function(Y, covar.matrix = matrix(1, length(Y)), H, variance, pval=TRUE) {

  use.houba <- is(H, "mmatrix")

  # in case there are some missing phenotypes, deal with it by 
  # extracting the corresponding lines
  w <- which(is.na(Y))
  if( length(w) > 0 ) {
    Y <- Y[-w]
    covar.matrix <- covar.matrix[-w, ]
    H <- H[-w,]
  }

  X <- as.matrix(covar.matrix)

  # if there's only one intercept, we use the dedicated function
  if(ncol(X) == 1 && all(X == 1)) {
    pi <- mean(Y)
    w <- pi*(1-pi)      
    Y1 <- Y - pi
    compute.var <- missing(variance)
    if(!compute.var) {
      if(length(variance) != ncol(H)) 
        stop("length(variance) should be equal to ncol(H)")
    }
    if(use.houba) {
      R <- as.data.frame(logitModelScore_nocovar_mmatrix(Y1, w, H, 0, ncol(H) - 1L, compute.var))
    } else {
      R <- as.data.frame(logitModelScore_nocovar_matrix(Y1, w, H, 0, ncol(H) - 1L, compute.var))
    }

    if(!compute.var) R$variance <- as.vector(variance)
  } else {
    # do we need to add an intercept
    # mean of (X beta - 1)^2 with beta = (X'X)⁻¹ X' 1 : residues of 1 ~ X ...
    if(mean( (X %*% solve( crossprod(X), colSums(X) ) - 1)**2 ) > 1e-8) {
      X <- cbind(1, X)
      warning("An intercept column was added to the covariate matrix")
    }
  
    # fit model under H0 and prepare second step
    fit <- glm(Y ~ X - 1, family = binomial())
    pi <- fit$fitted.values
    W <- pi*(1-pi)
    Y1 <- Y - pi
    WX <- W*X
    XWX <- crossprod(X, WX)
    ei <- eigen(XWX)
    ei$values <- 1/sqrt(ei$values)
    A <- tcrossprod(ei$values * t(ei$vectors), WX)
  
    # final step
    if(use.houba) {
      R <- as.data.frame(logitModelScore_mmatrix(Y1, W, A, H, 0, ncol(H) - 1L))
    } else {
      R <- as.data.frame(logitModelScore_matrix(Y1, W, A, H, 0, ncol(H) - 1L))
    }
  }
  R$z.value <- R$score / sqrt(R$variance)
  
  if(pval) { 
    R$p.left <- pnorm(R$z.value, lower.tail = TRUE)
    R$p.bilateral <- pchisq(R$z.value**2, df = 1, lower.tail = FALSE) 
    R$p.right <- pnorm(R$z.value, lower.tail = FALSE) 
  }
  
  R
}

