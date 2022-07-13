#' Logistic regression on HBD probability or FLOD score
#' 
#' @param x an atlas object
#' @param expl_var the explanatory variable 'FLOD' or 'HBD_prob'
#' @param covar_df a dataframe containing covariates
#' @param covar covariates of interest such as 'age', 'sex' , ...
#' if missing, all covariates of the dataframe are considered
#' @param n.cores number of cores for parallelization calculation (default = 1)
#' @param run whether the fonction is called or not (default = FALSE)
#' @param phen.code phenotype coding :
#'        - 'R' : 0:control ; 1:case ; NA:unknown (default)
#'        - 'plink' : 1:control ; 2:case ; 0/-9/NA:unknown
#' if 'plink' the function automatically convert it to 'R' to run logistic regression
#' 
#' @export

glm.HBD <- function( x, expl_var, covar_df, covar, n.cores = 1, run = FALSE, phen.code) {
  
  if(class(x)[1] != "atlas")
    stop("Need an atlas")
  
  if (run) { 
    if (expl_var == 'HBD_prob') {
      # Recovery HBD_prob
      hbd <- as.data.frame(x@HBD_recap)
    } else if (expl_var == 'FLOD') {
      # Recovery FLOD
      hbd <- as.data.frame(x@FLOD_recap)
    } else {
      stop("Explanatory variable must be 'HBD_prob' or 'FLOD'")
    }
    
    # Recovery phenotype
    id <- sub(".*:", "" , row.names(hbd))
    id.index <- match ( id, x@bedmatrix@ped$id )
    pheno <- x@bedmatrix@ped$pheno [id.index]
    if (phen.code == 'plink') {
      pheno <- ifelse(pheno == 1, 0, ifelse(pheno == 2, 1, NA))# Translate phenotype
    }
    
    # Recovery chr, snps, pos_cM and pos_Bp 
    final <- get.positions(x)	
    
    if(n.cores == 1 ) {
      # unadjusted 
      if (missing(covar_df)) {
        message("No covariates given for the analysis = unadjusted data. To use covariates import a dataframe.")
        message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i])"))
        for(i in 1:ncol(hbd)){
          model <- glm( pheno ~ hbd[,i] , family = binomial)
          final[i,'estimate'] 	<- summary(model)$coef[2,1]
          final[i,'std_error'] 	<- summary(model)$coef[2,2]
          final[i,'z_value'] 	<- summary(model)$coef[2,3]
          final[i,'p_value'] 	<- summary(model)$coef[2,4]
        }
        x@logisticRegression$unadj <- final
        message("-----------> GLM on UNADJUSTED data Done \n")
      }
      
      # adjusted 
      else {	
        if(missing(covar)) {
          message(paste0("No covariates specified - All covariates of the dataframe will be used : " , gsub(",", " +", toString(colnames(covar_df)))))
          message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i] + ", gsub(",", " +", toString(colnames(covar_df))) ,")" ))
          df <- na.omit(covar_df[id,])				 # take all covar given in the dataframe
        } else {
          message(paste0("Covariates = ", gsub(",", " +", toString(covar))))
          message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i] + ", gsub(",", " +", toString(covar)) ,")"))
          df <- na.omit(as.data.frame(covar_df[ id , covar])) #rownames covar_df  = individual id 	
        }
        
        for(i in 1:ncol(hbd)){
          model <- glm( pheno ~ hbd[,i] + . , data = df,  family = binomial)
          final[i,'estimate'] 	<- summary(model)$coef[2,1]
          final[i,'std_error'] 	<- summary(model)$coef[2,2]
          final[i,'z_value'] 	<- summary(model)$coef[2,3]
          final[i,'p_value'] 	<- summary(model)$coef[2,4]
        }
        x@logisticRegression$adj <- final 
        message("-----------> GLM on ADJUSTED data Done \n")
      }
    }
    
    # Parallelization - n.cores > 1
    else {
      ## Com Margot : à bouger :
      library(foreach)
      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
      
      if(missing(covar_df)) {
        message("No covariates given for the analysis = unadjusted data. To use covariates import a dataframe.")
        message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i])"))
        
        res <- foreach (i = 1:ncol(hbd), .combine = rbind) %dopar% {
          model <- glm( pheno ~ hbd[,i] ,  family = binomial)
          estimate <- summary(model)$coef[2,1]
          std_error <- summary(model)$coef[2,2]
          z_value <- summary(model)$coef[2,3]
          p_value <- summary(model)$coef[2,4]
          data.frame( estimate, std_error, z_value, p_value)
        }
        x@logisticRegression$unadj <- cbind(final, res) 
        message("-----------> GLM on UNADJUSTED data Done \n")
      }
      
      else {
        
        if(missing(covar)) {
          message(paste0("No covariates specified - All covariates of the dataframe will be used : " , gsub(",", " +", toString(colnames(covar_df)))))
          message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i] + ", gsub(",", " +", toString(colnames(covar_df))) ,")"))
          df <- na.omit(covar_df[id,])				 # take all covar given in the dataframe
        } else {
          message( paste0("Covariates = ", gsub(",", " +", toString(covar))))
          message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i] + ", gsub(",", " +", toString(covar)),")"))
          df <- na.omit(as.data.frame(covar_df[ id , covar]))  	
        } 
        
        res <- foreach (i = 1:ncol(hbd), .combine = rbind) %dopar% {
          model <- glm( pheno ~ hbd[,i] + . , data = df, family = binomial)
          estimate <- summary(model)$coef[2,1]
          std_error <- summary(model)$coef[2,2]
          z_value <- summary(model)$coef[2,3]
          p_value <- summary(model)$coef[2,4]
          data.frame( estimate, std_error, z_value, p_value)
        }
        x@logisticRegression$adj <- cbind(final, res)
        message("-----------> GLM on ADJUSTED data Done \n")
      }
      parallel::stopCluster(cl)
    }
  }
  x
} 
