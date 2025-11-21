#options(echo=FALSE);

#param       <- commandArgs(trailingOnly=T)
#filein      <- eval(paste(text=param[1])) 
#filout_pop  <- eval(paste(text=param[2])) 
#filout_ind  <- eval(paste(text=param[3])) 

#data=read.table(filein,h=T)
#data=subset(data,is.na(data$LIKELIHOOD)==F)


T <- function(alpha,L) {   ##T
	alpha[1]*10^L[1] + alpha[2]*10^L[2] + alpha[3]*10^L[3] + alpha[4]*10^L[4] + (1-sum(alpha))
}

P <- function(alpha,L) {   ##P
	alpha*10^L / sum(alpha*10^L) 
}

LOD <- function(alpha,L) {   ##fonction LOD
	# LOD=0
	# for (i in 1:nbind) {
		# LOD=LOD+log10(T(alpha,L[i,]))
	# }
	# -(LOD)
	-sum(log10(T(alpha,L)))
}

matingtype <- function(atlas){
  data <- atlas@estimations
  nbsubs <- ncol(data$l1)
  nbind <- nrow(data$l1)
  nameind <- atlas@submap_summary$id

  id <- diag(4)
  alpha <- matrix(rep(0,nbsubs*6),nbsubs,6)
  colnames(alpha)=c("SUB","1C","2C","2x1C","AV","OUT")
  PPost=NULL
  
  for (i in 1:nbsubs) {

  	L <- data.frame(matrix(c(data$l1[,i], data$l0[,i]),nbind,6))
  	L <- (L[,3:6]-L[,2])/log(10)
	  colnames(L) <- c("1C","2C","2x1C","AV")
	
	  res1 <- constrOptim(rep(0.20,4), LOD, NULL, ui=rbind(-id,id,-c(1,1,1,1)), ci=c(rep(-1,4),rep(0,4),-1), L=L)
	  res2 <- constrOptim(rep(0.01,4), LOD, NULL, ui=rbind(-id,id,-c(1,1,1,1)), ci=c(rep(-1,4),rep(0,4),-1), L=L)
	
  	res <- ifelse(LOD(res1$par,L) > LOD(res2$par,L), res2, res1)
	  alpha[i,] <- c(i,res$par,1-sum(res$par))
	
  	PPost <- rbind(PPost,t(apply(cbind(L,rep(log10(1),nbind)),1,P,alpha=alpha[i,2:6])))

  }	
  
  alpha_median <- ifelse (nbsubs>1, apply(alpha[,2:6],2,median), alpha[,2:6])
  alpha_median <- round(alpha_median/sum(alpha_median),4)
  alpha_median[5] <- round(1-sum(alpha_median[1:4]),4)
  
  if (nbsubs>1) {
	  PPost_median <- NULL
	  for (i in 1:nbind) {
		  vect <- apply(PPost[i+(0:(nbsubs-1))*nbind,],2,median)
		  vect <- round(vect/sum(vect),4)
		  vect[5] <- round(1-sum(vect[1:4]),4)
		  PPost_median <- rbind(PPost_median,vect) 
	  }
  } else {
	  PPost_median <- PPost
  }
  rownames(PPost_median)=nameind
  
  return(list(pop_mt = alpha_median, ind_par_mt = PPost_median))
}