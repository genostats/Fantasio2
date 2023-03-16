# Y = vecteur de 0 et 1
# covar.matrix = matrice de covariables [doit contenir un intercept]
# H = matrice des pHBD ou des FLOD...
 
glm.HBD.0 <- function(Y, covar.matrix = matrix(1, length(Y)), H, test = match.arg("bilateral", "right", "left")) {

  test <- match.arg(test)

  # déc QR de X pour améliorer la stabilité de l'algo
  if(all(covar.matrix[,1] == 1)) {
    qr.X <- qr(covar.matrix)
  } else {
    X1 <- cbind(1, covar.matrix);
    qr.X1 <- qr(X1);
    if(qr.X1$rank == ncol(X1)) {
      warning("An intercept column was added to the covariate matrix")
      qr.X <- qr.X1
    } else {
      qr.X <- qr(covar.matrix)
    }
  }
  X <- qr.Q(qr.X)

  # ajout d'une colonne vide !
  X <- cbind(X, 0)
  R <- as.data.frame(logitModel(Y, X, H, 0, ncol(H)-1))
  R$z.value <- R$beta/R$sd.beta
   
  if(test == "bilateral") 
    R$p <- pchisq(R$z.value**2, df = 1, lower.tail = FALSE) 
  else if(test == "right")
    R$p <- pnorm(R$z.value, lower.tail = FALSE) 
  else 
    R$p <- pnorm(R$z.value, lower.tail = TRUE)

  R
} 
